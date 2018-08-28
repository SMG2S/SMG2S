#!/usr/bin/env python

from distutils.core import setup, Extension


module = Extension('_smg2s',
                           sources=['smg2s_wrap.cxx'],
                           )

setup (name = 'smg2s',
       version = '1.0',
       author      = "Xinzhe Wu",
       description = """SMG2S classes""",
       ext_modules = [module],
       py_modules = ["smg2s"],
)

