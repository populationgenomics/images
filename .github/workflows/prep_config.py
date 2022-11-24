"""
Prepare ready contig TOML
"""

import os
import toml

IMAGES_PREFIX = os.environ['IMAGES_PREFIX']

d = toml.load('images.toml')
d = {k: f'{IMAGES_PREFIX}/{k}:{v}' for k, v in d.items()}
print(toml.dumps({'images': d}))
