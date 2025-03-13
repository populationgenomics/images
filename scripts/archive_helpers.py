import logging
import re
import subprocess
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from enum import Enum

import duckdb
import pandas as pd
import polars as pl
from google.api_core import exceptions as gcp_exceptions
from google.cloud import artifactregistry_v1

from scripts.Image import Image

logging.getLogger().setLevel(logging.INFO)

# Initialize the Artifact Registry client
artifact_client = artifactregistry_v1.ArtifactRegistryClient()


class Status(Enum):
    SUCCESS = 'success'
    FAILURE = 'failure'


def copy_image(source_image: Image, target_image: Image) -> tuple[Status, str]:
    """
    Copies a Docker image from one Google Cloud repository to another using skopeo.
    If the image already exists in the target repository, the function skips the copy and returns success.

    Parameters:
        source_image (str): Full path to the source image (including digest) like 'gcr.io/source-project/repo/image@sha256:xyz'.
        target_image (str): Full path to the target image like 'gcr.io/target-project/repo/image@sha256:xyz'.
    """
    try:
        # Check if the image already exists in the target repository
        inspect_command = ['skopeo', 'inspect', f'docker://{target_image}']
        logging.info(f'Checking if image {target_image} already exists...')
        result = subprocess.run(
            inspect_command, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
        )

        if result.returncode == 0:
            logging.info(f'Image {target_image} already exists. Skipping copy.')
            return Status.SUCCESS, ''

        # Copy the image if it doesn't exist
        copy_command = [
            'skopeo',
            'copy',
            '--all',
            f'docker://{source_image}',
            f'docker://{target_image}',
        ]
        logging.info(
            f'Copying image from {source_image} to {target_image} using skopeo...'
        )
        subprocess.run(copy_command, check=True)

        logging.info(
            f'Image successfully copied from {source_image} to {target_image}.'
        )
        return Status.SUCCESS, ''

    except subprocess.CalledProcessError as e:
        logging.error(f'Error occurred while copying image with skopeo: {e}')
        return Status.FAILURE, str(e)


def migrate_images(
    image_data,
    new_repo_parent,
    archive_prefix='ARCHIVED',
    repo_suffix='-archive',
    max_workers=10,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Archives images from the provided DuckDB relation by copying them to a new repository
    and adding a archive_prefix tag to the old images.

    Parameters:
        image_data (DuckDBPyRelation): The relation containing image data.
        new_repo_parent (str): The new repository parent, e.g., 'gcr.io/new-project'.
    """
    # Convert the DuckDB relation to a pandas DataFrame for easier iteration
    df = image_data.to_df()

    success = pd.DataFrame(columns=df.columns + ['message'])
    failure = pd.DataFrame(columns=df.columns + ['message'])

    # 'projects/sabrina-dev-337923/locations/australia-southeast1/repositories/'
    _, new_project, _, new_location, _ = new_repo_parent.split('/')

    # Function to handle image processing
    def process_image(row):
        source_image = Image(row)
        target_image = source_image.create_archival_image(new_repo_parent, repo_suffix)

        # Copy the image to the new repository
        status, message = copy_image(source_image, target_image)

        result = {'row': row, 'status': status}
        result['row']['message'] = message

        if status == Status.FAILURE:
            result['success'] = False
        else:
            result['success'] = True

            # Add the archive_prefix tag to the old image
            source_image.add_tag(f'{archive_prefix}-{source_image.image_digest_path}')

        return result

    # Create a ThreadPoolExecutor to parallelize the processing
    with ThreadPoolExecutor(
        max_workers=max_workers
    ) as executor:  # Adjust max_workers as necessary
        futures = {executor.submit(process_image, row): row for _, row in df.iterrows()}

        for future in as_completed(futures):
            result = future.result()
            row = result['row']
            if result['success']:
                success.loc[len(success)] = row
            else:
                failure.loc[len(failure)] = row

    return success, failure


def unmark_images(
    image_data, archive_prefix='ARCHIVED', max_workers=10
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Reverses the archiving process by removing any archive_prefix tags from images in the provided DuckDB relation.

    Parameters:
        image_data (DuckDBPyRelation): The relation containing image data.

    Returns:
        pd.DataFrame: A DataFrame containing rows where the archive_prefix tag was successfully removed.
    """
    # First, find all unique project, repo, name combos
    unique_images = duckdb.query(
        image_data,
        """
        SELECT DISTINCT image_location, image_project, image_repository, image_name
        FROM image_data
        """,
    )

    # Convert the DuckDB relation to a pandas DataFrame for easier iteration
    df = unique_images.to_df()

    success = pd.DataFrame(columns=df.columns + ['message'])
    failure = pd.DataFrame(columns=df.columns + ['message'])

    # Function to handle unarchiving process for each image
    def process_image(row, archive_prefix='ARCHIVED'):
        # Create the image reference without the digest to search for all versions of the image
        source_image = Image(row)

        result = {'row': row, 'success': False}

        logging.info(f'Processing image: {source_image}')

        try:
            # Get all tags for the given image name
            tags = source_image.search_base_image_tags(
                lambda t: t.startswith(archive_prefix)
            )

            for tag in tags:
                # Remove the archive_prefix tag
                source_image.remove_tag(tag)
                logging.info(f"Successfully removed tag '{tag}' from {source_image}.")

            result['success'] = True
        except subprocess.CalledProcessError as e:
            logging.error(
                f"Failed to remove '{archive_prefix}-*' tag(s) from {source_image}: {e}"
            )
            result['row']['message'] = str(e)
            result['success'] = False

        return result

    # Create a ThreadPoolExecutor to parallelize the processing
    with ThreadPoolExecutor(
        max_workers=max_workers
    ) as executor:  # Adjust max_workers as necessary
        futures = {
            executor.submit(process_image, row, archive_prefix): row
            for _, row in df.iterrows()
        }

        for future in as_completed(futures):
            result = future.result()
            row = result['row']
            if result['success']:
                success.loc[len(success)] = row
            else:
                failure.loc[len(failure)] = row

    return success, failure


def get_images_for_registry(
    registry: str,
    default_registry: str = 'projects/cpg-common/locations/australia-southeast1/repositories/',
):
    ar_client = artifactregistry_v1.ArtifactRegistryClient()

    if '/' not in registry:
        request = artifactregistry_v1.ListDockerImagesRequest(
            parent=f'{default_registry}{registry}'
        )
    else:
        request = artifactregistry_v1.ListDockerImagesRequest(parent=registry)

    # Make the request
    try:
        page_result = ar_client.list_docker_images(request=request)
    except Exception:
        page_result = []

    responses = []
    for response in page_result:
        name_match = re.search(
            'projects/([^/]+)/locations/([^/]+)/repositories/([^/]+)/dockerImages/([^@]+)@sha256:(.+)$',
            response.name,
        )

        if isinstance(name_match, re.Match):
            image_details = {
                'build_time': response.build_time,
                'image_size_bytes': response.image_size_bytes,
                'tags': response.tags,
                'update_time': response.update_time,
                'upload_time': response.upload_time,
                'name': response.name,
                'image_project': name_match[1],
                'image_location': name_match[2],
                'image_repository': name_match[3],
                'image_name': name_match[4],
                'image_digest': name_match[5],
            }

        responses.append(image_details)
    return responses


def get_image_data(image_data: dict[str, tuple[str, list[str]]]) -> dict:
    all_images = defaultdict(list)
    for key, data in image_data.items():
        parent, repos = data
        for repo in repos:
            images = get_images_for_registry(registry=repo, default_registry=parent)
            all_images[key].extend(images)

        all_images[key] = pl.DataFrame(
            [
                {
                    'build_time': i['build_time'],
                    'image_size_bytes': i['image_size_bytes'],
                    'tags': list(i['tags']),
                    'update_time': i['update_time'],
                    'upload_time': i['upload_time'],
                    'name': i['name'],
                    'image_project': i['image_project'],
                    'image_location': i['image_location'],
                    'image_repository': i['image_repository'],
                    'image_name': i['image_name'],
                    'image_digest': i['image_digest'],
                }
                for i in all_images[key]
            ]
        )

    return all_images


def delete_images(
    delete_images,
    new_repo_parent=None,
    archive_prefix='ARCHIVED',
    repo_suffix='-archive',
    dry_run=True,
    max_workers=10,
):
    """
    Deletes images from GCP Artifact Registry based on the given DataFrame.

    Args:
        delete_images (pd.DataFrame): DataFrame containing image details to delete.
                                      Expected columns: ['image_project', 'image_location', 'image_repository', 'image_name', 'image_digest']
        max_workers (int): Maximum number of threads to use for parallel processing.
    """

    # Convert the DuckDB relation to a pandas DataFrame for easier iteration
    df = delete_images.to_df()

    success = pd.DataFrame(columns=delete_images.columns + ['message'])
    failure = pd.DataFrame(columns=delete_images.columns + ['message'])

    def process_image(row):
        image = Image(row)

        # Delete tag
        delete_tag = f'{archive_prefix}-{image.digest}'
        existing_tags = image.tags.copy()
        if delete_tag not in existing_tags:
            existing_tags.append(delete_tag)

        # Check if the image exists
        if not image.exists():
            logging.info(
                f'Image {image.image_digest_path} does not exist. Skipping deletion.'
            )
            return {'row': row, 'success': True}

        result = {'row': row, 'success': False}

        if dry_run:
            logging.info(f'DRY RUN: Deleting image tags: {image.image_digest_path}...')
            result['row']['message'] = 'DRY RUN: Image not deleted.'
            result['success'] = True
            return result

        try:
            logging.info(f'Deleting image tags: {image.image_digest_path}...')
            image.remove_tags(existing_tags)
            message = image.delete_image()
            result['row']['message'] = message
            result['success'] = True
        except gcp_exceptions.NotFound as e:
            logging.warning(f'Image {image.image_digest_path} not found.')
            result['row']['message'] = str(e)
            result['success'] = True
        except gcp_exceptions.FailedPrecondition as e:
            logging.warning(
                f'Image {image.image_digest_path} failed precondition check: {e}.'
            )

            # We do not re-add the archival tag if a precondition fails
            existing_tags.remove(delete_tag)
            image.add_tags(existing_tags)

            # Now, we remove the archival (copy) image
            if new_repo_parent:
                archival_image = image.create_archival_image(
                    new_repo_parent, repo_suffix
                )
                message = archival_image.delete_image()

            result['row']['message'] = str(e) + f'\n\n{message}'
            result['success'] = False
        except subprocess.CalledProcessError as e:
            logging.error(f'Failed to delete {image.image_digest_path}. Error: {e}')
            logging.error(f'Readding back the image Tags: {existing_tags}')
            image.add_tags(existing_tags)
            result['row']['message'] = str(e)
            result['success'] = False

        return result

    # Create a ThreadPoolExecutor to parallelize the processing
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(process_image, row): row for _, row in df.iterrows()}

        for future in as_completed(futures):
            result = future.result()
            row = result['row']
            if result['success']:
                success.loc[len(success)] = row
            else:
                failure.loc[len(failure)] = row

    return success, failure


def list_duplicates(
    images,
    new_repo_parent,
    repo_suffix='-archive',
    max_workers=10,
):
    # Convert the DuckDB relation to a pandas DataFrame for easier iteration
    df = images.to_df()

    dupes = pd.DataFrame(columns=images.columns + ['message'])
    no_dupes = pd.DataFrame(columns=images.columns + ['message'])

    def process_images(row):
        image = Image(row)
        archival_image = image.create_archival_image(new_repo_parent, repo_suffix)
        result = {'row': row, 'dupe': False}

        # Check if the image exists
        if not image.exists():
            logging.info(f'Image {image.image_digest_path} does not exist. ')
            return result

        # Check if the archival image exists
        if archival_image.exists():
            result['dupe'] = True
            result['row']['message'] = (
                f'Archival image {archival_image.image_digest_path} exists.'
            )
            return result

        result['row']['message'] = (
            f'Archival image {archival_image.image_digest_path} does not exist.'
        )
        return result

    # Create a ThreadPoolExecutor to parallelize the processing
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(process_images, row): row for _, row in df.iterrows()
        }

        for future in as_completed(futures):
            result = future.result()
            row = result['row']
            if result['dupe']:
                dupes.loc[len(dupes)] = row
            else:
                no_dupes.loc[len(no_dupes)] = row

    return dupes, no_dupes
