.PHONY: compile lint install install-dev install-all

compile-requirements:
	docker run --platform linux/amd64 -v $$(pwd):/opt/deps python:3.11 /bin/bash -c '\
		cd /opt/deps; \
		pip install pip-tools; \
		pip-compile requirements.in;\
		pip-compile requirements-dev.in;\
	'

lint:
	ruff check
	ruff format

install:
	pip install -r requirements.txt

install-dev:
	pip install -r requirements-dev.txt

