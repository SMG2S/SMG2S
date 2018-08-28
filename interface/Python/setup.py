#!/usr/bin/env python

from distutils.core import setup, Extension
import mpi4py

mpi_incdir = mpi4py.get_include()


module = Extension('_smg2s',
                           sources=['smg2s_wrap.cxx'],
                           include_dirs=[mpi_incdir],
                           )

setup (name = 'smg2s',
       version = '1.0',
       author      = "Xinzhe Wu",
       description = """SMG2S classes""",
       ext_modules = [module],
       py_modules = ["smg2s"],
)

