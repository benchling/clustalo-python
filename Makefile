ROOT_DIR:=$(shell dirname $(realpath $(firstword $(MAKEFILE_LIST))))

build:
	docker build -t clustalo-build .
	docker run --rm -v ${ROOT_DIR}:/mnt clustalo-build cp -r /src/dist/ /mnt/

upload: build
	python -m twine upload dist/*
