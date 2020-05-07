.PHONY: build
ROOT_DIR:=$(shell dirname $(realpath $(firstword $(MAKEFILE_LIST))))

build:
	# Build bdist_wheel for local MacOS development - turn off OPENMP on Mac
	OPENMP_DISABLED=true python3 setup.py bdist_wheel
	OPENMP_DISABLED=true python2 setup.py bdist_wheel
	python3 setup.py sdist

upload: build
	python -m twine upload dist/*
