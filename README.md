# Images

This is a repository of single-software images to build and push in the `cpg-common` artifact registry.
This might also be a good place to put scripts for moving containers with Skopeo.

The directory should be named the same as the resulting image on `cpg-common`.

## Image moving scripts

Please put your scripts that move images to our cpg-common repository in the `/scripts` directory.

They should follow the format: `move-{IMAGE_NAME}:{IMAGE_VERSION}.sh`:

```shell
#!/usr/bin/env bash
set -ex
SOURCE_IMAGE=""
IMAGE_NAME=""
IMAGE_TAG=""

gcloud auth configure-docker australia-southeast1-docker.pkg.dev
skopeo copy ${SOURCE_IMAGE} docker://australia-southeast1-docker.pkg.dev/cpg-common/images/${IMAGE_NAME}:${IMAGE_TAG}
```

Then if you could follow the convention for the analysis-runner:

```shell
IMAGE_NAME="TRTools"
IMAGE_TAG="v4.0.2"
analysis-runner \
    --dataset fewgenomes --access-level standard \
    --description "Images: Move ${IMAGE_NAME}:${IMAGE_TAG}" \
    --output-dir images/mv-${IMAGE_NAME}-${IMAGE_TAG} \
    scripts/move-${IMAGE_NAME}:${IMAGE_TAG}.sh
```


## Single tool images

You should create a folder called `images/${IMAGE_NAME}`, with any relevant files. You must name the docker file: `Dockerfile`.
Once the image is ready to be built, you should create a pull-request, this will automatically assign code-reviewers.

Once your branch has been merged, you can then deploy your image using the GitHub Action, providing the name and tag of the image.

This will automatically build, and push the image to the `cpg-common` artifact registry.
