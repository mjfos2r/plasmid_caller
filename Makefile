IMAGE_NAME = plasmid_caller
VERSION := $(shell python -c 'from $(IMAGE_NAME).__about__ import __version__; print(__version__)')

TAG1 = mjfos2r/$(IMAGE_NAME):$(VERSION)
TAG2 = mjfos2r/$(IMAGE_NAME):latest

all: | refresh build push

refresh:
	uv lock
	uv sync --refresh

build:
	docker build -t $(TAG1) -t $(TAG2) .

build_no_cache:
	docker build --no-cache -t $(TAG1) -t $(TAG2) .

tag:
	git tag -s v$(VERSION) -m "Release version $(VERSION)"

push:
	docker push $(TAG1)
	docker push $(TAG2)
