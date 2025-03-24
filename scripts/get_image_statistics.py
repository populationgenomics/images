from pathlib import Path

import polars as pl
from common.image_logs_helpers import get_image_logs
from common.image_repository_helpers import Image, list_images_in_repository

repositories = ['images', 'images-dev', 'images-tmp', 'images-archive']


def get_image_stats():
    # get the logs first as it is much faster than getting the images
    logs_df = get_image_logs()

    all_images: list[Image] = []
    for repo in repositories:
        all_images.extend(list_images_in_repository(repo))

    # get the image data, use "active" fields so it can be best compared to the logs
    images_df = pl.DataFrame(
        [
            {
                'full_path': image.active_full_path,
                'build_time': image.build_time,
                'update_time': image.update_time,
                'upload_time': image.upload_time,
                'size_bytes': image.size_bytes,
                'tags': image.tags,
                'digest': image.digest,
                'project': image.project,
                'location': image.location,
                'repository': image.active_repository,
                'name': image.name,
                'status': image.status,
                'short_path': image.active_short_path,
            }
            for image in all_images
        ]
    )
    logs_df.write_parquet(
        Path(__file__).parent / 'image_statistics/src/data/logs.parquet'
    )
    images_df.write_parquet(
        Path(__file__).parent / 'image_statistics/src/data/images.parquet'
    )


if __name__ == '__main__':
    get_image_stats()
