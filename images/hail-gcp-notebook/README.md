# Hail GCP Notebook


HGVS has this akward psycopg2 dependency which we won't use in the image, and the binary equivalent (`psycopg2-binary`) is already installed. 

```bash
pip-compile | grep -v "psycopg2==" > requirements.txt 
```
