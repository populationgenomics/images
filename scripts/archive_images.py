import logging
import re
from pathlib import Path

from common.image_repository_helpers import (
    DeleteVersionStatus,
    Image,
    Repository,
    add_tags,
    copy_image,
    delete_version,
    list_images_in_repository,
    version_exists,
)

logging.getLogger().setLevel(logging.INFO)

SUPPORTED_REPOSITORIES = ['images']


def get_archive_list():
    archive_list_file = Path(__file__).parent.parent / 'archived_images.txt'
    archive_list: list[str] = []
    with open(archive_list_file) as f:
        for line in f:
            archive_list.append(line.strip())
    # remove duplicates from the list and return
    return list(set(archive_list))


def validate_archive_list(archive_list: list[str]):
    for image in archive_list:
        match = re.search(
            '([^/]+)/([^@]+)@sha256:(.+)$',
            image,
        )
        if not match:
            raise ValueError(f'Invalid image specified: {image}')

        repo = match[1]
        name = match[2]
        digest = match[3]

        if not repo or not name or not digest:
            raise ValueError(f'Invalid image specified: {image}')

        if repo not in SUPPORTED_REPOSITORIES:
            raise ValueError(f'Unsupported repository: {repo}')


def move_image(source_image: Image, dest_image: Image, dest_repository: Repository):
    # Check if the tags for the source image already exist in the destination repository
    # this shouldn't happen but is worth checking
    conflicting_tags = dest_repository.find_conflicting_tags(dest_image)
    copy_created = False

    if conflicting_tags:
        raise Exception(
            f'Tag(s) {",".join(conflicting_tags)} for image {dest_image.name} already exist in the destination repository.'
        )

    # Copy the image to the destination repository
    if dest_repository.includes_image_version(dest_image):
        logging.info(
            f'Image {dest_image.version_id} already exists in the destination repository. Skipping copy.'
        )
    else:
        copy_image(source_image.docker_name, dest_image.docker_name)
        copy_created = True

    # add the tags to the copied image
    add_tags(dest_image, dest_image.tags)

    # double check that the new version exists before deleting the old one
    exists = version_exists(dest_image)
    if not exists:
        raise Exception(
            f'Image {dest_image.version_id} does not exist in the destination repository. Something has gone very wrong'
        )
    delete_op_status = delete_version(source_image)
    # If the delete operation fails with a failed precondition error then we want to
    # delete the copy that was just created, but only if it was actually just created
    # if it already existed then it is a bit too risky to delete it.
    if delete_op_status == DeleteVersionStatus.FAILED_PRECONDITION and copy_created:
        logging.warning(
            f'Deletion newly created copy of image at {dest_image.version_id} as the source image failed to delete'
        )
        exists = version_exists(source_image)
        if not exists:
            raise Exception(
                f'Image {source_image.version_id} does not exist in the source repository. Deletion failed so it should be there, something has gone wrong so bailing here to avoid deleting the destination image'
            )
        delete_version(dest_image)


def archive_images_in_repository(repository: str, archive_list: list[str]):
    logging.info('Getting images from repositories, this takes a while...')
    active_images = list_images_in_repository(repository)
    archived_images = list_images_in_repository(f'{repository}-archive')

    to_archive = [
        image for image in active_images if image.active_version_id in archive_list
    ]
    to_unarchive = [
        image
        for image in archived_images
        if image.active_version_id not in archive_list
    ]

    logging.info(f'Found {len(to_archive)} images to archive.')
    logging.info(f'Found {len(to_unarchive)} images to unarchive.')

    for image in to_archive:
        logging.info(f'Archiving image: {image.active_version_id}')
        move_image(
            image, image.convert_to_archived(), Repository(images=archived_images)
        )

    for image in to_unarchive:
        logging.info(f'Unarchiving image: {image.active_version_id}')
        move_image(image, image.convert_to_active(), Repository(images=active_images))

    logging.info('Done!')


def archive_images():
    archive_list = get_archive_list()
    validate_archive_list(archive_list)

    for repository in SUPPORTED_REPOSITORIES:
        to_archive = [img for img in archive_list if img.startswith(f'{repository}/')]
        archive_images_in_repository(repository, to_archive)


if __name__ == '__main__':
    archive_images()
