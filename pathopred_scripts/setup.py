from setuptools import setup, Extension

from Cython.Build import cythonize

import numpy



extension = Extension('_load_data', ['_load_data.pyx'],
                      extra_compile_args="-O2 -march=native -pipe -mtune=native".split(),
                      extra_link_args="-O2 -march=native -pipe -mtune=native".split())
setup(name='ss_load_data', ext_modules=cythonize(extension), include_dirs=[numpy.get_include()])

