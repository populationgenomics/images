# Images

This is a repository of:

- Single-software images to build and push in the `cpg-common` artifact registry.
- Scripts that move external images into our artifact registry using Skopeo.

This is a BAD place for:

- multi-software images
- analysis-specific images, eg: an image that contains in-progress software for processing data

## Single tool images

Assuming you want to create an image for a tool called `mytool`, version `1.0.1`, add a folder called `images/mytool` with a `Dockerfile`. Parametrise the `Dockerfile` by `${VERSION}`, that will correspond to the version of your tool and the tag of your future image. For example, inside the `Dockerfile` you would have:

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
