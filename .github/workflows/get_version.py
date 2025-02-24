#!/usr/bin/env python3
import json
import subprocess
import re
import os
import sys


def get_changed_files():
    # Get the list of files changed between HEAD~1 and HEAD
    result = subprocess.run(
        ["git", "diff", "HEAD~1", "HEAD", "--name-only"], capture_output=True, text=True
    )
    return result.stdout.splitlines()


def get_file_diff(file_path):
    # Get the diff for the given file
    result = subprocess.run(
        ["git", "diff", "HEAD~1", "HEAD", "--", file_path],
        capture_output=True,
        text=True,
    )
    return result.stdout


def extract_version_from_diff(diff_text):
    # Look for an added line containing "ARG VERSION=" and extract the version value.
    pattern = re.compile(r"^\+\s*ENV\s+VERSION\s*=\s*([^\s]+)", re.MULTILINE)
    match = pattern.search(diff_text)
    return match.group(1) if match else None


def get_next_version_tag(full_image_name, version):
    """
    Call GCP to list all tags for the given image and determine the next
    available version suffix.
    """
    # Command to list tags (using Container Registry command; adjust if needed)
    cmd = [
        "gcloud",
        "container",
        "images",
        "list-tags",
        full_image_name,
        "--format=json",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        # Command failed (image might not exist yet): assume no tags exist.
        return f"{version}-1"
    try:
        tags_list = json.loads(result.stdout)
    except Exception:
        tags_list = []

    max_suffix = 0
    # Look for tags that exactly match: version-N (e.g. 1.0.0-1, 1.0.0-2)
    pattern = re.compile(rf"^{re.escape(version)}-(\d+)$")
    for entry in tags_list:
        # Each entry is expected to have a "tags" key that is a list of tags.
        tags = entry.get("tags", [])
        for tag in tags:
            match = pattern.match(tag)
            if match:
                num = int(match.group(1))
                if num > max_suffix:
                    max_suffix = num
    new_suffix = max_suffix + 1
    return f"{version}-{new_suffix}"


def main():
    changed_files = get_changed_files()
    include_entries = []
    BASE_IMAGE_PATH = os.environ.get(
        "GCP_BASE_IMAGE", "australia-southeast1-docker.pkg.dev/cpg-common/images"
    )

    for file in changed_files:
        if file.endswith("Dockerfile"):
            diff_text = get_file_diff(file)
            version = extract_version_from_diff(diff_text)
            if version:
                # The folder name is extracted from the Dockerfile's directory.
                folder_path = os.path.dirname(file)
                folder = os.path.basename(folder_path) if folder_path else "root"
                # Build the fully qualified image name.
                full_image_name = f"{BASE_IMAGE_PATH}/{folder}"
                # Get next available version tag from GCP.
                new_tag = get_next_version_tag(full_image_name, version)
                include_entries.append({"name": folder, "tag": new_tag})

    # Build the final matrix structure.
    matrix = {"include": include_entries}

    print(str(matrix).replace(" ", ""), end="")


main()
