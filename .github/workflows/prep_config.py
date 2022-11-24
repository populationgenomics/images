"""
Prepare ready contig TOML
"""

import os
from pathlib import Path
import toml

IMAGES_PREFIX = os.environ['IMAGES_PREFIX']

with Path('images.toml').open('rb') as f:
    d = toml.load(f)
    d = {k: f'{IMAGES_PREFIX}/k:{v}' for k, v in d.items()}

print(toml.dumps({'images': d}))
