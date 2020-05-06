.PHONY: build
ROOT_DIR:=$(shell dirname $(realpath $(firstword $(MAKEFILE_LIST))))

build:
	python3 setup.py sdist
	# Build bdist_wheel for local MacOS development - turn off OPENMP on Mac (but leave it alone for sdist)
	OPENMP_DISABLED=true python3 setup.py sdist bdist_wheel
	OPENMP_DISABLED=true python2 setup.py bdist_wheel

upload: build
	python -m twine upload dist/*
