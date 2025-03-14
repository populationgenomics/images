.PHONY: compile lint install install-dev install-all

compile:
	pip-compile --output-file=requirements.txt pyproject.toml
	pip-compile --extra=dev --output-file=requirements-dev.txt pyproject.toml

lint:
	ruff lint
	ruff format

install:
	pip install -r requirements.txt

install-dev:
	pip install -r requirements-dev.txt

install-all: install install-dev
