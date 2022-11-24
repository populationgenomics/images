"""
Prepare ready contig TOML
"""

import os
from pathlib import Path
import tomllib

IMAGES_PREFIX = os.environ['IMAGES_PREFIX']

with Path('images.toml').open('tb') as f:
    d = tomllib.load(f)
    d = {k: f'{IMAGES_PREFIX}/k:{v}' for k, v in d.items()}

print(tomllib.dumps({'images': d}))
