IMAGE_NAME := freesurfer/petsurfer-bids
TAG := $(shell python3 -c "import tomllib; print(tomllib.load(open('pyproject.toml','rb'))['project']['version'])")

.PHONY: docker

docker:
	docker build -t $(IMAGE_NAME):$(TAG) .
