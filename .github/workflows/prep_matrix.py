"""
Prepare test matrix (to build images in parallel)
"""

from pathlib import Path
import tomllib

with Path('images.toml').open('rb') as f:
    d = tomllib.load(f)
before_d = {}
if (path := Path('before/images.toml')).exists():
    before_d = tomllib.load(path.open('rb'))
d = {
    'include': [
        {'name': name, 'tag': tag}
        for name, tag in d.items()
        if not before_d.get(name) or before_d[name] != tag
    ]
}
print(str(d).replace(' ', ''), end='')
