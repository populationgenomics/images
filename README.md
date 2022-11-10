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
ARG VERSION=${VERSION}
RUN pip install mytool==${VERSION}
```

Then, add an entry into [`images.toml`](images.toml) with the name of your image matching the desired tag/version:

```toml
mytool = '1.0.1'
```

Finally, create a pull-request with your changes and assign someone to review it. Once merged, the GitHub CI will automatically build your image and push it into the `cpg-common` artifact registry.

To update the tool version, modify the corresponding entry in TOML and create a pull request, CI will automatically rebuild the image.

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
