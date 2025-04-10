.PHONY: compile-requirements lint install install-dev

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
	pip install --no-deps -r requirements.txt

install-dev:
	pip install --no-deps -r requirements-dev.txt
