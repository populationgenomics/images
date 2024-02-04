# Hail GCP Notebook


HGVS has this akward psycopg2 dependency which we won't use in the image, and the binary equivalent (`psycopg2-binary`) is already installed.

The `requirements.in` specifies the high level requirements we want to install - `pip-compile` resolves this into a complete list (effectively a `pip freeze`), we remove `psycopg2` specifically, and use that to install with `--no-deps`.

```bash
pip install pip-tools

# update / freeze deps
pip-compile | grep -v "psycopg2==" > requirements.txt
```
