from setuptools import setup, Extension

import os

OPENMP_DISABLED = os.environ.get('OPENMP_DISABLED', False)
libraries=['clustalo', 'stdc++']
extra_compile_args = []

if not OPENMP_DISABLED:
    libraries.append('gomp')
    extra_compile_args.append('-fopenmp')

module = Extension('clustalo',
                   sources = ['clustalo.c'],
                   include_dirs=['/usr/local/sci/clustalo/current/include/clustalo'],
                   library_dirs=['/usr/local/sci/clustalo/current/lib'],
                   libraries=libraries,
                   extra_compile_args=extra_compile_args)

setup(name='clustalo',
      version='0.1.3',
      description='Python wrapper around libclustalo',
      author='Benchling Engineering',
      author_email='eng@benchling.com',
      url='https://github.com/benchling/clustalo-python',
      ext_modules=[module])
