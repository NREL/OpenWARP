#!/usr/bin/env python
"""
This module defines configuration variables used by the application

Changes in version 1.1 (Code Acceleration of the Calculation of Influence Coefficients of Nemoh):
    Added USE_ODE_INFLUENCE_COEFFICIENTS settings

Changes in version 1.2 (Implementation of Higher Order Panel Methods):
    Added USE_HIGHER_ORDER, NUM_PANEL_HIGHER_ORDER, B_SPLINE_ORDER settings

Changes in version 1.3 (Dipoles Implementation in NEMOH):
    Added USE_DIPOLES_IMPLEMENTATION and THIN_PANELS settings

Changes in version 1.4 (Hydrodynamic Data Exporter Assembly v1.0)
       Added parameters controlling wether or not to compute the drift forces or yaw moment

Changes in version 1.5 (Irregular Frequencies Assembly)
       Added parameters to remove irregular frequencies or not.

Changes in version 1.6 (OpenWarp - Add Logging Functionality)
       Added support for logging.
"""


__author__ = "yedtoss"
__copyright__ = "Copyright (C) 2014-2016 TopCoder Inc. All rights reserved."
__version__ = "1.6"

import os

import numpy as np
np.set_printoptions(threshold=10)


BASE_TEST_DIR = os.path.join('..', 'test_files', 'remove_irregular_frequencies')

HDF5_FILE = os.path.join(BASE_TEST_DIR, 'db.hdf5')
"""
The path to the hdf5 file.
String, non-empty, not null
Required.
"""

NEMOH_CALCULATIONS_FILE = os.path.join(BASE_TEST_DIR, 'Nemoh.cal')
"""
The path to the nemoh calculations.
String.
Not Required.
If present a conversion from this nemoh file to the hdf5 file is performed
"""

NEMOH_INPUT_FILE = os.path.join(BASE_TEST_DIR, 'input.txt')
"""
The path to the nemoh input.
String.
Not Required.
If present a conversion from this nemoh file to the hdf5 file is performed
"""


MESH_TEC_FILE = os.path.join(BASE_TEST_DIR, 'mesh', 'mesh.tec')
"""
The path to the mesh data in tec format
String.
Not Required
If not present no export is performed
"""

FK_FORCE_TEC_FILE = os.path.join(BASE_TEST_DIR, 'results', 'fkforce.tec')
"""
The path to the froudkrylov force data in tec format
String.
Not Required
If not present no export is performed
"""

RADIATION_COEFFICIENTS_TEC_FILE = os.path.join(BASE_TEST_DIR, 'results', 'radiationcoefficients.tec')
"""
The path to the file where to save the added mass and damping forces for the radiation problems
in tec format.
String.
Not Required
If not present no export is performed
"""

DIFFRACTION_FORCE_TEC_FILE = os.path.join(BASE_TEST_DIR, 'results', 'diffractionforce.tec')
"""
The path to the file where to save the diffraction force for the diffraction problems
in tec format.
String.
Not Required
If not present no export is performed
"""

EXCITATION_FORCE_TEC_FILE = os.path.join(BASE_TEST_DIR, 'results', 'excitationforce.tec')
"""
The path to the file where to save the the excitation force for the diffraction problems
in tec format.
String.
Not Required
If not present no export is performed
"""

IRF_TEC_FILE = os.path.join(BASE_TEST_DIR, 'results', 'irf.tec')
"""
The path to the file where to save the the irf in tec format.
String.
Not Required
If not present no export is performed
"""


WAVE_FIELD_TEC_FILE = os.path.join(BASE_TEST_DIR, 'results', 'WaveField.tec')
"""
The path to the file where to save the the wave field in tec format.
String.
Not Required
If not present no export is performed
"""


NEMOH_INT = 'i'
"""
Represents the integer type to use when performing computations.
It should be a valid numpy integer type.
See http://docs.scipy.org/doc/numpy/reference/arrays.scalars.html#arrays-scalars-built-in
for possible values
"""

NEMOH_FLOAT = 'f'
"""
Represents the float type to use when performing computations.
It should be a valid numpy float type.
See http://docs.scipy.org/doc/numpy/reference/arrays.scalars.html#arrays-scalars-built-in
for possible values
"""


NEMOH_COMPLEX = 'F'
"""
Represents the complex type to use when performing computations.
It should be a valid numpy complex type.
See http://docs.scipy.org/doc/numpy/reference/arrays.scalars.html#arrays-scalars-built-in
for possible values
"""

GREEN_TABULATION_NUMX = 328
"""
Represents Number of points in x direction of tabulated data
Set to negative value of None to not use
Default is 328
"""


GREEN_TABULATION_NUMZ = 46
"""
Represents Number of points in z direction of tabulated data
Set to negative value of None to not use
Default 46
"""

GREEN_TABULATION_SIMPSON_NPOINTS = 251
"""
Represents Number of sub intervals used to approximate the green function integral using simpson rule
Set to negative value of None to not use
"""

USE_ODE_INFLUENCE_COEFFICIENTS = False
"""
Indicate whether or not to use the ode method to compute the influence coefficients
Default False
Set to None to not use.
Corresponding setting in hdf5 file is required
"""

USE_HIGHER_ORDER = False
"""
Whether or not to use higher order panel method.
Default False
Set to None to not use.
Corresponding setting in hdf5 file is required
"""

NUM_PANEL_HIGHER_ORDER = 1
"""
The number of panel per patch in the higher order method
Default 1
Set to None to not use.
Corresponding setting in hdf5 file is required
"""

B_SPLINE_ORDER = 1
"""
The order of the B-Spline for the potential in the higher order
panel method
Default 1
Set to None to not use.
Corresponding setting in hdf5 file is required
"""

USE_DIPOLES_IMPLEMENTATION = False
"""
Wether or not to use the dipoles Implementation
Default False
Set to None to not use it.
Corresponding setting in hdf5 file is required
"""

THIN_PANELS = [-1]
"""
A list containing the index of the panel which are thin dipoles.
The index should be 0-based.
Set the list to [-1] to indicate that all panels are thin dipoles.
Set to None to not use it.
Corresponding setting in hdf5 file is required
"""

COMPUTE_DRIFT_FORCES = False
"""
Wether or not to compute the drift forces.
Set to None not to use it
Corresponding setting in hdf5 file is required
"""

COMPUTE_YAW_MOMENT = False
"""
Wether or not to compute the yaw moment.
Set to None not to use it
Corresponding setting in hdf5 file is required
"""

REMOVE_IRREGULAR_FREQUENCIES = True
"""
Wether or not to remove the irregular frequencies.
Set to None not to use it
Corresponding setting in hdf5 file is required
"""

LOGGING_CONFIGURATION_FILE = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'logging.json')

LOG_FILE = os.path.join(BASE_TEST_DIR, 'logs', 'logs.log')
