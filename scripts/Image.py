import logging
from dataclasses import dataclass
from typing import Callable

from google.api_core import exceptions as gcp_exceptions
from google.cloud import artifactregistry_v1

logging.getLogger().setLevel(logging.INFO)

# Initialize the Artifact Registry client
artifact_client = artifactregistry_v1.ArtifactRegistryClient()


@dataclass
class Image:
    name: str
    project: str
    location: str
    repository: str
    digest: str
    tags: list[str]

    # Optional fields
    build_time: str | None
    size_bytes: str | None
    update_time: str | None
    upload_time: str | None

    def __init__(self, image_data):
        self.name = image_data.get('image_name')
        self.project = image_data.get('image_project')
        self.location = image_data.get('image_location')
        self.repository = image_data.get('image_repository')
        self.digest = image_data.get('image_digest')
        self.tags = list(image_data.get('tags'))

        self.build_time = image_data.get('build_time')
        self.size_bytes = image_data.get('image_size_bytes')
        self.update_time = image_data.get('update_time')
        self.upload_time = image_data.get('upload_time')

    def clone(self):
        return Image(
            {
                'image_name': self.name,
                'image_project': self.project,
                'image_location': self.location,
                'image_repository': self.repository,
                'image_digest': self.digest,
                'tags': self.tags,
                'build_time': self.build_time,
                'image_size_bytes': self.size_bytes,
                'update_time': self.update_time,
                'upload_time': self.upload_time,
            }
        )

    def __str__(self):
        return self.image_digest_path

    @property
    def image_path(self):
        return f'{self.location}-docker.pkg.dev/{self.project}/{self.repository}/{self.name}'

    @property
    def image_digest_path(self):
        # Create the image reference without the digest to search for all versions of the image
        image_path = f'{self.image_path}@sha256:{self.digest}'
        image_path = image_path.replace('%2F', '/')

        return image_path

    @property
    def gcp_repository_resource_name(self):
        return f'projects/{self.project}/locations/{self.location}/repositories/{self.repository}'

    @property
    def gcp_package_resource_name(self):
        return f'{self.gcp_repository_resource_name}/packages/{self.name}'

    @property
    def gcp_version_resource_name(self):
        return f'{self.gcp_package_resource_name}/versions/sha256:{self.digest}'

    def get_gcp_tag_resource_name(self, tag):
        return f'projects/{self.project}/locations/{self.location}/entrygroups/{self.repository}/entries/{self.digest}/tags/{tag}'

    def get_gcp_tag_resource_name_for_delete(self, tag):
        return f'{self.gcp_package_resource_name}/tags/{tag}'

    def get_image_with_tag(self, tag):
        return f'{self.image_path}:{tag}'

    def list_tags(self):
        "List all tags for the given image name."
        request = artifactregistry_v1.ListTagsRequest(
            parent=self.gcp_package_resource_name,
        )

        tags = []
        for tag in artifact_client.list_tags(request=request):
            if self.digest in tag.version:
                tags.append(tag.name.split('/')[-1])

        return tags

    def search_base_image_tags(self, filter_fxn: Callable[[str], bool] | None = None):
        "Get all tags for the given image name."

        # List tags for the specified image
        request = artifactregistry_v1.ListTagsRequest(
            parent=self.gcp_package_resource_name,
        )

        tags = []
        for tag in artifact_client.list_tags(request=request):
            tags.append(tag.name.split('/')[-1])  # Extract the tag name

        if not filter_fxn:
            return tags

        return list(filter(filter_fxn, tags))

    def exists(self):
        try:
            request = artifactregistry_v1.GetVersionRequest(
                name=self.gcp_version_resource_name,
            )
            artifact_client.get_version(request=request)
            return True
        except gcp_exceptions.NotFound:
            logging.info(f'Image {self.image_digest_path} does not exist.')
            return False
        except Exception as e:
            logging.error(f'Error checking if image exists: {e}')
            return False

    def image_tag_exists(self, tag):
        "Check if this tag exists for the image in the regsitry already."
        raise NotImplementedError

    def add_tag(self, tag):
        "Add a tag to a container image using Google Cloud Artifact Registry client."
        try:
            # Create the request to add the tag
            request = artifactregistry_v1.CreateTagRequest(
                parent=self.gcp_package_resource_name,
                tag_id=tag,
                tag=artifactregistry_v1.Tag(
                    name=self.get_gcp_tag_resource_name(tag),
                    version=self.gcp_version_resource_name,
                ),
            )

            # Call the API to add the tag
            artifact_client.create_tag(request=request)

            logging.info(
                f'Successfully added tag {tag} to image {self.image_digest_path}'
            )
        except gcp_exceptions.AlreadyExists as already_exists_err:
            logging.warning(
                f'Tag {tag} already exists for image {self.image_digest_path}: {already_exists_err}'
            )
        except Exception as e:
            logging.error(
                f'Failed to add tag {tag} to image {self.image_digest_path}: {e}'
            )
            raise

    def add_tags(self, tags: list[str]):
        "Add multiple tags to a container image using gcloud CLI."
        for tag in tags:
            self.add_tag(tag)

    def remove_tag(self, tag):
        "Remove a tag from the image using Google Cloud Artifact Registry client."
        try:
            logging.info(self.get_gcp_tag_resource_name_for_delete(tag))
            request = artifactregistry_v1.DeleteTagRequest(
                name=self.get_gcp_tag_resource_name_for_delete(tag),
            )
            artifact_client.delete_tag(request=request)
            logging.info(
                f'Successfully removed tag {tag} from image {self.image_digest_path}'
            )
        except gcp_exceptions.NotFound:
            logging.warning(f'Tag {tag} not found for image {self.image_digest_path}')
        except Exception as e:
            logging.error(
                f'Failed to remove tag {tag} from image {self.image_digest_path}: {e}'
            )
            raise

    def remove_tags(self, tags):
        "Remove multiple tags from the image using gcloud CLI."
        for tag in tags:
            self.remove_tag(tag)

    def delete_image(self):
        "Delete a container image using Google Cloud Artifact Registry client."

        assert self.exists(), f'Image {self.image_digest_path} does not exist.'

        try:
            request = artifactregistry_v1.DeleteVersionRequest(
                name=self.gcp_version_resource_name,
            )
            operation = artifact_client.delete_version(request=request)

            # Wait for the operation to complete
            operation.result()

            logging.info(f'Successfully deleted image: {self.image_digest_path}')
        except gcp_exceptions.FailedPrecondition as failed_precondition_err:
            # https://cloud.google.com/sdk/gcloud/reference/artifacts/versions/delete
            # FailedPrecondition possible reasons:
            # - Image is referenced by another image manifest or index
            # - Image is tagged (only deletes untagged images by default)
            # - Image needs to be kept as it meets the retention policy
            # - Retention policy restricts deletion
            # - Oranization policy restricts deletion
            # - Image does not exist
            # - You do not have the relevent permission to delete this image version
            # The relevant error message will be raised to be handled by the caller
            logging.warning(
                f'Image {self.image_digest_path} failed precondition check: {failed_precondition_err}'
            )
            raise failed_precondition_err
        except Exception as e:
            logging.error(f'Failed to delete image {self.image_digest_path}: {e}')
            raise e

    def create_archival_image(
        self, new_repo_parent: str, repo_suffix: str = '-archive'
    ) -> 'Image':
        "Create a new image reference with the archive_prefix tag."
        new_image: Image = self.clone()

        # 'projects/sabrina-dev-337923/locations/australia-southeast1/repositories/'
        _, new_project, _, new_location, _ = new_repo_parent.removesuffix('/').split(
            '/'
        )

        new_image.location = new_location
        new_image.project = new_project
        new_image.repository = self.repository + repo_suffix
        return new_image
