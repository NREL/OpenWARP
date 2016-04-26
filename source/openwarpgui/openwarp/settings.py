# -*- coding: utf-8 -*-
"""
This Python module contains configurations(aka settings) of the application.

Updated since version 1.1:
    1. Added path to paraview script template file.

Updated since version 1.2 (OpenWarp - Add Logging Functionality):
    1. Added support for logging.

"""

__author__ = "caoweiquan322, yedtoss"
__copyright__ = "Copyright (C) 2014-2016 TopCoder Inc. All rights reserved."
__version__ = "1.2"

import sys
import os

# Represents the suffix of executables in different platforms.
bin_suffix = sys.platform == 'win32' and '.exe' or ''

# Represents the directory of this application.
app_dir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))

# Represents the application home in the user home
user_home = os.path.join(os.path.expanduser("~"), "OpenWarpFiles")

# Represents the full path of the "nglib-mesh" executable binary.
# Required, non-empty.
if sys.platform == 'win32':
    MESH_GENERATOR_BIN = os.path.join(app_dir, 'bundled', 'mesh-generator', 'bin', 'nglib-mesh' + bin_suffix)
elif sys.platform == "darwin":
    MESH_GENERATOR_BIN = os.path.join(app_dir, 'bundled', 'mesh-generator', 'build', 'nglib-mesh' + bin_suffix)
else:
	MESH_GENERATOR_BIN = os.path.join(app_dir, 'bundled', 'mesh-generator', 'bin', 'nglib_mesh' + bin_suffix)

# Represents the full path of the "preProcessor" executable binary.
# Required, non-empty.
PREPROCESSOR_BIN = os.path.join(app_dir, 'bundled', 'simulation', 'bin', 'preProcessor' + bin_suffix)

# Represents the full path of the "Solver" executable binary.
# Required, non-empty.
SOLVER_BIN = os.path.join(app_dir, 'bundled', 'simulation', 'bin', 'Solver' + bin_suffix)

# Represents the full path of the "postProcessor" executable binary.
# Required, non-empty.
POSTPROCESSOR_BIN = os.path.join(app_dir, 'bundled', 'simulation', 'bin', 'postProcessor' + bin_suffix)

# Represents the full path of the ParaView executable binary.
# Required, non-empty.
if sys.platform == 'win32':
    PARAVIEW_BIN = os.path.join(app_dir, 'bundled', 'paraview', 'bin', 'paraview' + bin_suffix)
elif sys.platform == "darwin":
    PARAVIEW_BIN = os.path.join(app_dir, 'bundled', 'paraview.app/Contents/MacOS/paraview')
else:
	PARAVIEW_BIN = os.path.join(app_dir, 'bundled', 'paraview_linux/bin/paraview')


# Represents the full path of the ParaView script template for loading data files.
# Required, non-empty.
PARAVIEW_SCRIPT_TEMPLATE = os.path.join(app_dir, 'bundled', 'paraview', 'paraview_script_template.py')

# Represents the root directory for user data.
# Required, non-empty.
USER_DATA_DIRECTORY = os.path.join(user_home, 'user_data')

# Represents the root directory for temporary files
# Required, non-empty.
TEMP_DATA_DIRECTORY = os.path.join(user_home, 'temp')

# Represents the accepted file extensions for paraview to visualize.
VISUALIZATION_FILE_EXTENSIONS = ['stl', 'vtk', 'tec']

# Represents the web server port.
# Required. It must be a positive integer in range of [0, 65535].
WEB_SERVER_PORT = 8386

LOG_DIR = os.path.join(user_home, 'logs')

LOG_FILE = os.path.join(LOG_DIR, 'logs.log')
