"""
Prepare ready contig TOML
"""

import os
import sys

import toml

IMAGES_PREFIX = os.environ['IMAGES_PREFIX']

d = toml.load('images.toml')
# Manually adding cpg_workflows image, as we don't maintain it with
# images.toml, instead we check out the latest version from the
# production-pipelines repository.
d['cpg_workflows'] = sys.argv[1]
d = {k: f'{IMAGES_PREFIX}/{k}:{v}' for k, v in d.items()}
print(toml.dumps({'images': d}))
print(toml.dumps({'images': d}))
