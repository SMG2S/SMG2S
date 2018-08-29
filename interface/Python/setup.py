#!/usr/bin/env python

from distutils.core import setup, Extension
import mpi4py
from glob import glob

mpi_incdir = mpi4py.get_include()
print(mpi_incdir)

with open('README.txt') as file:
    long_description = file.read()

module = Extension('_smg2s',
                           sources=['smg2s/smg2s_wrap.cxx'],
                           include_dirs=[mpi_incdir, './smg2s/include'],
                           swig_opts=['-I./smg2s/include'],
                           extra_compile_args=['-stdlib=libc++', '-std=c++0x'],
                           )

setup (name = 'smg2s',
        version = '1.0.1',
       author      = "Xinzhe Wu",
       author_email = 'xinzhe.wu1990@gmail.com',
       description = 'SMG2S: Scalable Matrix Generator with Given Spectrum',
       long_description = long_description,
       ext_modules = [module],
       py_modules = ["smg2s/smg2s"],
       url = 'http://smg2s.github.io',
       license = 'GNU Lesser General Public License v3.0',
)

