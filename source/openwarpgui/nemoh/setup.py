#!/usr/bin/env python
"""
This is the standard python setup.py file to install this module
"""

import sys
import subprocess
import os
import glob
import shutil

from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import numpy as np

__author__ = "yedtoss"
__copyright__ = "Copyright (C) 2014 TopCoder Inc. All rights reserved."
__version__ = "1.0"

args = sys.argv[1:]

# Make a `cleanall` rule to get rid of intermediate and library files
if "cleanall" in args:
    print "Deleting cython files..."
    # Just in case the build directory was created by accident,
    for name in glob.glob('*.so'):
      os.remove(name)
    for name in glob.glob('*.pyd'):
      os.remove(name)
    for name in glob.glob('*.pyc'):
      os.remove(name)
    shutil.rmtree('build', True)
    for name in glob.glob('*.c'):
      os.remove(name)

    # Now do a normal clean
    sys.argv[1] = "clean"

setup(
  name = 'Nemoh Python wrapper',
  ext_modules = cythonize([
    Extension("solver_fortran", ["solver_fortran.pyx"],
              libraries=["nemoh"],
              include_dirs=[np.get_include()],
              # Disabled as do not work in all os
              #library_dirs=[os.path.join(os.getcwd(), 'bundled', 'simulation', 'libs')],
              #runtime_library_dirs=[os.path.join(os.getcwd(), 'bundled', 'simulation', 'libs')]
    ),

    ]),
)