#!/usr/bin/env python3
#  ruff: noqa: S603

import json
import subprocess
import re
import os
import sys


def extract_version_from_file(file_path: str) -> str | None:
    """
    Extract the version from a Dockerfile by searching for a line like:
    ENV VERSION=1.0.0
    """
    with open(file_path) as f:
        content = f.read()
    pattern = re.compile(r'^\s*ENV\s+VERSION\s*=\s*([^\s]+)', re.MULTILINE)
    match = pattern.search(content)
    return match.group(1) if match else None


def get_next_version_tag(folder: str, version: str) -> str:
    """
    Query GCP to list tags for the given image and determine the next available
    version suffix for the extracted version.
    """
    base_image_path_prod = os.environ.get(
        'GCP_BASE_IMAGE',
        'australia-southeast1-docker.pkg.dev/cpg-common/images',
    )
    base_image_path_archive = os.environ.get(
        'GCP_BASE_ARCHIVE_IMAGE',
        'australia-southeast1-docker.pkg.dev/cpg-common/images-archive',
    )
    full_image_name_prod = f'{base_image_path_prod}/{folder}'
    full_image_name_archive = f'{base_image_path_archive}/{folder}'

    tags_list = []
    import logging

    logging.basicConfig(level=logging.ERROR)
    for full_image_name in [full_image_name_prod, full_image_name_archive]:
        cmd = [
            'gcloud',
            'container',
            'images',
            'list-tags',
            full_image_name,
            '--format=json',
        ]
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)

        # If existing tags are found, proceed to determine the next version.
        tags_list += json.loads(result.stdout)

    max_suffix = 0
    pattern = re.compile(rf'^{re.escape(version)}-(\d+)$')

    # If no tags are found, it returns the next version as -1.
    for entry in tags_list:
        tags = entry.get('tags', [])
        for tag in tags:
            match = pattern.match(tag)
            if match:
                num = int(match.group(1))
                max_suffix = max(max_suffix, num)
    new_suffix = max_suffix + 1
    return f'{version}-{new_suffix}'


def get_before_commit():
    """
    Determines the correct 'before' commit:
    - If on a feature branch: Finds the commit where the branch
    diverged from `origin/main`
    - If on `main` with a merge: Finds the last two merge
    commits and compares them.
    """
    current_branch = subprocess.run(
        ['git', 'rev-parse', '--abbrev-ref', 'HEAD'],
        capture_output=True,
        text=True,
        check=True,
    ).stdout.strip()

    if current_branch == 'main':
        # Get the last two merge commits
        merge_commits = subprocess.run(
            ['git', 'log', '--merges', '--format=%H', '-n', '2'],
            capture_output=True,
            text=True,
            check=True,
        ).stdout.splitlines()
        merge_commit_threshold = 2
        if len(merge_commits) >= merge_commit_threshold:
            return merge_commits[1]  # Compare latest merge with the one before it
        return 'HEAD~1'  # Default to previous commit if no merges
    # Find the commit where this branch diverged from main
    return subprocess.run(
        ['git', 'merge-base', 'HEAD', 'origin/main'],
        capture_output=True,
        text=True,
        check=True,
    ).stdout.strip()


def main():
    before_commit = get_before_commit()

    # Get changed Dockerfiles
    result = subprocess.run(
        ['git', 'diff', '--name-only', before_commit, 'HEAD', '--', '*Dockerfile'],
        capture_output=True,
        text=True,
        check=True,
    )
    dockerfiles = result.stdout.splitlines()

    include_entries = []

    for file in dockerfiles:
        current_version = extract_version_from_file(file)
        if current_version is None:
            continue

        # Get only the last folder name
        folder_path = os.path.dirname(file)
        folder = os.path.basename(folder_path) if folder_path else 'root'

        # Determine the next available tag based on current_version.
        new_tag = get_next_version_tag(folder, current_version)

        include_entries.append({'name': folder, 'tag': new_tag})

    # Build the final matrix structure.
    matrix = {'include': include_entries}
    print(json.dumps(matrix, separators=(',', ':')), file=sys.stderr)
    print(json.dumps(matrix, separators=(',', ':')), end='')


main()
