#!/usr/bin/env python3
import json
import subprocess
import re
import os


def extract_version_from_file(file_path):
    """
    Extract the version from a Dockerfile by searching for a line like:
      ARG VERSION=1.0.0
    """
    try:
        with open(file_path, 'r') as f:
            content = f.read()
        pattern = re.compile(r'^\s*ENV\s+VERSION\s*=\s*([^\s]+)', re.MULTILINE)
        match = pattern.search(content)
        return match.group(1) if match else None
    except Exception:
        return None


def get_next_version_tag(folder, version):
    """
    Query GCP to list tags for the given image and determine the next available
    version suffix for the extracted version.

    This example uses the 'gcloud container images list-tags' command.
    Adjust it if you are using Artifact Registry.
    """
    # Construct the fully qualified image name.
    BASE_IMAGE_PATH_PROD = os.environ.get(
        'GCP_BASE_IMAGE',
        'australia-southeast1-docker.pkg.dev/cpg-common/images',
    )
    BASE_IMAGE_PATH_ARCHIVE = os.environ.get(
        'GCP_BASE_IMAGE',
        'australia-southeast1-docker.pkg.dev/cpg-common/images-archive',
    )
    full_image_name_prod = f'{BASE_IMAGE_PATH_PROD}/{folder}'
    full_image_name_archive = f'{BASE_IMAGE_PATH_ARCHIVE}/{folder}'

    tags_list = []

    for full_image_name in [full_image_name_prod, full_image_name_archive]:
        cmd = [
            'gcloud',
            'container',
            'images',
            'list-tags',
            full_image_name,
            '--format=json',
        ]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            return f'{version}-1'
        try:
            tags_list += json.loads(result.stdout)
        except Exception:
            pass

    max_suffix = 0
    pattern = re.compile(rf'^{re.escape(version)}-(\d+)$')
    for entry in tags_list:
        # Each entry should have a list of tags.
        tags = entry.get('tags', [])
        for tag in tags:
            match = pattern.match(tag)
            if match:
                num = int(match.group(1))
                if num > max_suffix:
                    max_suffix = num
    new_suffix = max_suffix + 1
    return f'{version}-{new_suffix}'


def main():
    # Path to the checkout of the commit before the push.
    before_dir = os.environ.get('BEFORE_DIR', 'before')
    include_entries = []

    # List Dockerfiles tracked in the current commit.
    result = subprocess.run(
        ['git', 'ls-files', '*Dockerfile'], capture_output=True, text=True
    )
    dockerfiles = result.stdout.splitlines()

    for file in dockerfiles:
        current_version = extract_version_from_file(file)
        if current_version is None:
            continue

        before_file = os.path.join(before_dir, file)
        before_version = None
        if os.path.exists(before_file):
            before_version = extract_version_from_file(before_file)

        # We add the entry if the version is new or has changed.
        if before_version != current_version:
            folder_path = os.path.dirname(file)
            folder = os.path.basename(folder_path) if folder_path else 'root'

            # Determine the next available tag based on current_version.
            new_tag = get_next_version_tag(folder, current_version)

            include_entries.append({'name': folder, 'tag': new_tag})

    # Build the final matrix structure.
    matrix = {'include': include_entries}

    print(str(matrix).replace(' ', ''), end='')


main()
