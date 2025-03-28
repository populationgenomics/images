import datetime
import re
from dataclasses import dataclass
from typing import Literal

from google.api_core.datetime_helpers import DatetimeWithNanoseconds
from google.cloud import artifactregistry_v1
from google.cloud.artifactregistry_v1.types import DockerImage
from google.protobuf.timestamp_pb2 import Timestamp


# google's python types autogenerated from protobufs have incorrect types
# for the timestamp fields, it says they are a "Timestamp" but they are actually
# a DatetimeWithNanoseconds | None, so need a function to handle the conversion
def image_timestamp_to_datetime(timestamp: Timestamp | DatetimeWithNanoseconds | None):
    if timestamp is None:
        return None
    return datetime.datetime.fromisoformat(str(timestamp))


@dataclass
class Image:
    full_path: str
    build_time: datetime.datetime | None
    update_time: datetime.datetime | None
    upload_time: datetime.datetime | None
    size_bytes: int
    tags: list[str]
    digest: str
    project: str
    location: str
    repository: str
    name: str

    @property
    def status(self) -> Literal['active', 'archived']:
        return 'archived' if self.repository.endswith('-archive') else 'active'

    @property
    def short_path(self):
        return f'{self.repository}/{self.name}'

    # Compute what some paths would have been before image was archived
    # This is useful for comparing image versions across archived and non-archived
    # repositories.
    @property
    def active_full_path(self):
        return self.full_path.replace('-archive/dockerImages', '/dockerImages')

    @property
    def active_repository(self):
        return self.repository.removesuffix('-archive')

    @property
    def active_short_path(self):
        return f'{self.active_repository}/{self.name}'

    @property
    def active_version_id(self):
        return f'{self.active_short_path}@sha256:{self.digest}'

    @property
    def active_docker_name(self):
        return f'{self.location}-docker.pkg.dev/{self.project}/{self.active_repository}/{self.name}@sha256:{self.digest}'

    # Compute what some paths would look like after the image is archived
    # This is useful when moving images back and forth between archived and active
    @property
    def archived_full_path(self):
        return (
            self.full_path.replace('/dockerImages', '-archive/dockerImages')
            if '-archive/dockerImages' not in self.full_path
            else self.full_path
        )

    @property
    def archived_repository(self):
        return (
            f'{self.repository}-archive'
            if not self.repository.endswith('-archive')
            else self.repository
        )

    @property
    def archived_short_path(self):
        return f'{self.archived_repository}/{self.name}'

    @property
    def archived_version_id(self):
        return f'{self.archived_short_path}@sha256:{self.digest}'

    @property
    def archived_docker_name(self):
        repository = self.archived_repository
        return f'{self.location}-docker.pkg.dev/{self.project}/{repository}/{self.name}@sha256:{self.digest}'

    def convert_to_active(self):
        "Get a version of the image tas if it were active"
        return Image(
            full_path=self.active_full_path,
            build_time=self.build_time,
            update_time=self.update_time,
            upload_time=self.upload_time,
            size_bytes=self.size_bytes,
            tags=self.tags,
            digest=self.digest,
            project=self.project,
            location=self.location,
            repository=self.active_repository,
            name=self.name,
        )

    def convert_to_archived(self):
        "Get a version of the image as if it were archived"
        return Image(
            full_path=self.archived_full_path,
            build_time=self.build_time,
            update_time=self.update_time,
            upload_time=self.upload_time,
            size_bytes=self.size_bytes,
            tags=self.tags,
            digest=self.digest,
            project=self.project,
            location=self.location,
            repository=self.archived_repository,
            name=self.name,
        )

    @staticmethod
    def from_artifact_repository_image(image_data: DockerImage) -> 'Image':
        """
        Given a DockerImage object, create a new Image object.
        """
        image_path = image_data.name.replace('%2F', '/')

        name_match = re.search(
            'projects/([^/]+)/locations/([^/]+)/repositories/([^/]+)/dockerImages/([^@]+)@sha256:(.+)$',
            image_path,
        )

        if not name_match:
            raise ValueError(f'Invalid image path: {image_path}')

        project = name_match[1]
        location = name_match[2]
        repository = name_match[3]
        name = name_match[4]
        digest = name_match[5]

        # Validate that info was pulled from the image path successfully
        if not (
            type(project) is str
            and type(location) is str
            and type(repository) is str
            and type(name) is str
            and type(digest) is str
        ):
            raise ValueError(f'Invalid image path: {image_path}')

        return Image(
            full_path=image_path,
            build_time=image_timestamp_to_datetime(image_data.build_time),
            update_time=image_timestamp_to_datetime(image_data.update_time),
            upload_time=image_timestamp_to_datetime(image_data.upload_time),
            size_bytes=image_data.image_size_bytes,
            # Filter out digest from tags, no need to double up on it
            tags=[l for l in list(image_data.tags) if l != digest],
            digest=digest,
            project=project,
            location=location,
            repository=repository,
            name=name,
        )


def list_images_in_repository(repository: str):
    """
    Get a list of all images in the specifiec repository.
    returns a list of Image dataclass instances
    """
    ar_client = artifactregistry_v1.ArtifactRegistryClient()

    request = artifactregistry_v1.ListDockerImagesRequest(
        parent=f'projects/cpg-common/locations/australia-southeast1/repositories/{repository}',
    )

    page_result = ar_client.list_docker_images(request=request)
    return [Image.from_artifact_repository_image(response) for response in page_result]
