# -*- coding: utf-8 -*-
"""
Copyright (C) 2014 TopCoder Inc., All Rights Reserved.

This Python module provides a template of script for ParaView to load multiple data files.
This script could NOT run directly, but should be treated as a template and run after edited some fields.

@author:  caoweiquan322
@version: 1.0
"""
import os.path
# Check if paraview.simple module has been imported.
try:
    paraview.simple
except:
    from paraview.simple import *

# Reset cameras
paraview.simple._DisableFirstRenderCameraReset()

def render_data_file(filepath):
    try:
        filename = os.path.basename(filepath)
        if filename.endswith('vtk'):
            reader = LegacyVTKReader(FileNames = [filepath], guiName = filename)
            Show()
        elif filename.endswith('stl'):
            reader = STLReader(FileNames = [filepath], guiName = filename)
            Show()
        elif filename.endswith('tec'):
            reader = TecplotReader(FileNames = [filepath], guiName = filename)
            # Tec files are generally not shown.
        else:
            pass
    except:
        # Yield
        pass

# The whole process to launch
# parameter_files is meant to be replaced by a real list of filepath.
files = <parameter_files>
for filepath in files:
    render_data_file(filepath)
Render()

