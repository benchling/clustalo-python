import logging
import os
import platform
import shutil
import subprocess

from setuptools import setup, Extension

OPENMP_DISABLED = os.environ.get('OPENMP_DISABLED', False)
libraries = ['clustalo', 'stdc++']


class BuildFlags:
    _tools_formulae = {
        'Darwin': ('brew', 'libomp'),
    }
    _args = {
        'Darwin': {
            'compiler': ['-Xpreprocessor', '-fopenmp'],
            'linker': ['-lomp']
        }
    }
    _default_args = {
        'compiler': ['-fopenmp'],
        'linker': ['-fopenmp']
    }
    compiler: list
    linker: list

    def __init__(self):
        self._system = platform.system()
        tool, formula = self._tools_formulae.get(self._system, ('apt', 'libomp-dev'))
        args = self._args.get(self._system, self._default_args)

        not_found = self._libomp_check(tool, formula)
        if not_found is not None:
            logging.warning(f'{repr(not_found)} not found -- cannot compile parallelized code')
            for key in args:
                args[key] = []

        for key, val in args.items():
            self.__setattr__(key, val)

    @staticmethod
    def _libomp_check(tool, formula):
        if shutil.which(tool) is None:
            return tool

        formulae = subprocess.check_output([tool, 'list']).decode()
        if formula not in formulae:
            return formula

        return None


flags = BuildFlags()
module = Extension('clustalo',
                   sources=['clustalo.c'],
                   include_dirs=['/usr/include/clustalo', '/usr/local/include/clustalo'],
                   library_dirs=['/usr/local/lib'],
                   libraries=libraries,
                   extra_compile_args=flags.compiler,
                   extra_link_args=flags.linker)

setup(name='clustalo',
      version='0.1.2',
      description='Python wrapper around libclustalo',
      author='Benchling Engineering',
      author_email='eng@benchling.com',
      url='https://github.com/benchling/clustalo-python',
      ext_modules=[module])
