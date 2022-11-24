"""
Prepare test matrix (to build images in parallel)
"""

import sys
from pathlib import Path
import toml

matrix = toml.load('images.toml')
before_d = {}
if Path(path := 'before/images.toml').exists():
    before_d = toml.load(path)

# Changed only
matrix = {
    name: tag
    for name, tag in matrix.items()
    if not before_d.get(name) or before_d[name] != tag
}
if matrix:
    matrix = {'include': [{'name': name, 'tag': tag} for name, tag in matrix.items()]}
else:
    matrix = ''

print(str(matrix).replace(' ', ''), end='', file=sys.stderr)
print(str(matrix).replace(' ', ''), end='')
