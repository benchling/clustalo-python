from distutils.core import setup, Extension

module = Extension('clustalo',
                   sources = ['clustalo.c'],
                   include_dirs=['/usr/include/clustalo', '/usr/local/include/clustalo'],
                   libraries=['clustalo', 'stdc++', 'gomp'],
                   extra_compile_args=['-fopenmp'])

setup(name = 'Clustalo',
      version = '0.1',
      description = 'Python wrapper around Clustal Omega',
      ext_modules = [module])
