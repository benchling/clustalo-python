from distutils.core import setup, Extension

module = Extension('_clustalo',
                   sources = ['clustalo/clustalo.c'],
                   include_dirs=['/usr/include/clustalo', '/usr/local/include/clustalo'],
                   libraries=['clustalo', 'stdc++', 'gomp'],
                   extra_compile_args=['-fopenmp'])

setup(name='clustalo',
      version='0.1',
      description='Python wrapper around libclustalo',
      packages=['clustalo'],
      ext_modules=[module])
