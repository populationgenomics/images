"""
Prepare test matrix (to build images in parallel)
"""

import sys
from pathlib import Path
import toml

d = toml.load('images.toml')
before_d = {}
if Path(path := 'before/images.toml').exists():
    before_d = toml.load(path)

# Changed only
d = {
    name: tag
    for name, tag in d.items()
    if not before_d.get(name) or before_d[name] != tag
}
if d:
    d = {'include': [{'name': name, 'tag': tag} for name, tag in d.items()]}
else:
    d = {}

print(str(d).replace(' ', ''), end='', file=sys.stderr)
print(str(d).replace(' ', ''), end='')
