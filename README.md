# Images

This is a repository of:

- Single-software images to build and push in the `cpg-common` artifact registry.
- Scripts that move external images into our artifact registry using `skopeo`.

This is a BAD place for:

- multi-software images
- analysis-specific images, e.g., an image that contains an in-progress script for processing data

## Single tool images

Assuming you want to create an image for a tool called `mytool` of version `1.0.1`, add a folder called `images/mytool` with a `Dockerfile` file. Parametrise the `Dockerfile` by `${VERSION}` that will correspond to the version of your tool and the tag of your future image. For example, inside the `Dockerfile` you would have:

```Dockerfile
ARG VERSION=1.0.1
RUN pip install mytool==${VERSION}
```

Finally, create a pull-request with your changes and assign someone to review it. Once merged, the GitHub CI will automatically build your image and push it into the `cpg-common` artifact registry.
The image will be given a tag like `1.0.1-1`, `1.0.1-2`, etc, made up from your upstream tool version number with an automatically incremented sequence number suffix.
This means that you will have a permanent tag to use to refer to each build of your image even when the upstream tool version has not changed.

To update the tool version, modify the corresponding entry in the `ARG VERSION` line and create a pull request, CI will automatically rebuild the image.
If the upstream tool is unchanged but you are modifying the packaging in the _Dockerfile_, just make those modifications and create a pull request. CI will increment the sequence number suffix automatically.

## Image moving scripts

Whenever possible, avoid moving images with `scopeo`, and build them with a `Dockefile` instead following the section above. Moving is not recommended for two reasons:

- It hides how the image was build originally (and sometimes the original `Dockefile` is not even shared),
- The source image can be removed by the author in the future, making your analysis no longer reproducible.

If you still need to move, create a script named `move-{IMAGE_NAME}:{IMAGE_VERSION}.sh` following the template below:

```shell
#!/usr/bin/env bash
set -ex
SOURCE_IMAGE=""
IMAGE_NAME=""
IMAGE_TAG=""

gcloud auth configure-docker australia-southeast1-docker.pkg.dev
skopeo copy ${SOURCE_IMAGE} docker://australia-southeast1-docker.pkg.dev/cpg-common/images/${IMAGE_NAME}:${IMAGE_TAG}
```

And submit it with the analysis runner:

```shell
IMAGE_NAME="mytool"
IMAGE_TAG="1.0.1"
analysis-runner \
    --dataset fewgenomes --access-level standard \
    --description "Images: Move ${IMAGE_NAME}:${IMAGE_TAG}" \
    --output-dir images/mv-${IMAGE_NAME}-${IMAGE_TAG} \
    scripts/move-${IMAGE_NAME}:${IMAGE_TAG}.sh
```
