IMAGE_NAME := freesurfer/petsurfer-bids
TAG := latest

.PHONY: docker

docker:
	docker build -t $(IMAGE_NAME):$(TAG) .
