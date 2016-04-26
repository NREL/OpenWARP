#!/usr/bin/env python
"""
This module defines variables exposing the structure of the hdf5 file used by the application.
Group names will and should always be terminated by /
Datasets name should not be terminated by /

Changes in Version 1.1 (Drift forces and QTF Implementation of Nemoh)
                  Added structure path to store drift forces and yaw moment

Changes in version 1.2 (Code Acceleration of the Calculation of Influence Coefficients of Nemoh):
    Added structure to store whether to run ode influence coefficients or not

Changes in version 1.3 (Implementation of Higher Order Panel Methods):
    Added structure to store USE_HIGHER_ORDER, NUM_PANEL_HIGHER_ORDER, B_SPLINE_ORDER settings

Changes in version 1.4 (Dipoles Implementation in NEMOH):
    Added structure to store USE_DIPOLES_IMPLEMENTATION and THIN_PANELS settings

Changes in version 1.5 (Hydrodynamic Data Exporter Assembly v1.0)
       Added parameters to store excitation forces, radiation damping coefficients, added mass,
       as well as mesh properties like center of buoyancy, volume displacement.

Changes in version 1.6 (Irregular Frequencies Assembly)
       Added parameters to remove irregular frequencies or not: H5_SOLVER_REMOVE_IRREGULAR_FREQUENCIES
       and another controlling if a panel is in the body domain or interior free surface domain:
       H5_SOLVER_IS_INTERIOR_DOMAIN.

Changes in version 1.7 (OpenWarp - Add Logging Functionality)
       Added support for logging.
"""

import numpy as np

__author__ = "yedtoss"
__copyright__ = "Copyright (C) 2014-2016 TopCoder Inc. All rights reserved."
__version__ = "1.7"


#
#
#
#   INPUTS
#
#
#

H5_INPUT_GROUP = 'input/'
"""
The hdf5 top level input group
"""

H5_CALCULATIONS = H5_INPUT_GROUP + 'calculations/'
"""
The hdf5 input calculations group. Contains mostly all the values previously in nemoh.cal file
"""

H5_INPUT_SOLVER = H5_INPUT_GROUP + 'solver/'
"""
The hdf5 input solver group. Contains all values previously in input.txt
"""

H5_SOLVER_USE_HIGHER_ORDER = H5_INPUT_SOLVER + 'USE_HIGHER_ORDER'
"""
Whether or not to use higher order panel method. 1 for using it 0 for no
"""

H5_SOLVER_REMOVE_IRREGULAR_FREQUENCIES = H5_INPUT_SOLVER + 'REMOVE_IRREGULAR_FREQUENCIES'
"""
Whether or not to remove irregular frequencies. 1 for using it 0 for no
"""

H5_SOLVER_USE_DIPOLES_IMPLEMENTATION = H5_INPUT_SOLVER + 'USE_DIPOLES_IMPLEMENTATION'
"""
Whether or not to use dipoles Implementation. 1 for using it 0 for no.
"""

H5_SOLVER_COMPUTE_YAW_MOMENT = H5_INPUT_SOLVER + 'COMPUTE_YAW_MOMENT'
"""
Whether or not to compute yaw moment. 1 for using it 0 for no.
"""

H5_SOLVER_COMPUTE_DRIFT_FORCES = H5_INPUT_SOLVER + 'COMPUTE_DRIFT_FORCES'
"""
Whether or not to compute drift forces. 1 for using it 0 for no.
"""


H5_SOLVER_THIN_PANELS = H5_INPUT_SOLVER + 'THIN_PANELS'
"""
Array containing whether a given panel is a conventional or a dipole one's.
1 for dipoles panel, 0 for conventional
"""

H5_SOLVER_IS_INTERIOR_DOMAIN = H5_INPUT_SOLVER + 'IS_INTERIOR_DOMAIN'
"""
Array indicating whether or not a panel is in the body or the interior free surface domain. 1 to indicate the interior, 0 for the body.
"""


H5_SOLVER_NUM_PANEL_HIGHER_ORDER = H5_INPUT_SOLVER + 'NUM_PANEL_HIGHER_ORDER'
"""
The number of panel per patch in the higher order method
"""

H5_SOLVER_B_SPLINE_ORDER = H5_INPUT_SOLVER + 'B_SPLINE_ORDER'
"""
The order of the B-Spline for the potential in the higher order
"""

H5_SOLVER_SWITCH_ODE_INFLUENCE = H5_INPUT_SOLVER + 'SWITCH_ODE_INFLUENCE'
"""
Indicate whether or not to use the ode method to compute the influence coefficients
"""

H5_SOLVER_GREEN_TABULATION_NUMX = H5_INPUT_SOLVER + 'GREEN_TABULATION_NUMX'
"""
Number of points in x direction of tabulated data
"""

H5_SOLVER_GREEN_TABULATION_NUMZ = H5_INPUT_SOLVER + 'GREEN_TABULATION_NUMZ'
"""
Number of points in z direction of tabulated data
"""

H5_SOLVER_GREEN_TABULATION_SIMPSON_NPOINTS = H5_INPUT_SOLVER + 'GREEN_TABULATION_SIMPSON_NPOINTS'
"""
Number of sub intervals used to approximate the green function integral
using simpson rule
"""

H5_SOLVER_TYPE = H5_INPUT_SOLVER + 'type'
"""
The solver type. (0) for Direct Gauss (1) for GMRES (2) GMRES with FMM acceleration (2 not implemented yet)
Previouly line 2 of input.txt
"""

H5_SOLVER_GMRES_RESTART = H5_INPUT_SOLVER + 'gmres_restart'
"""
The Restart parameter for GMRES.
Previously line 3 of input.txt
"""

H5_SOLVER_GMRES_STOPPING = H5_INPUT_SOLVER + 'gmres_stopping'
"""
Stopping criterion for GMRES
Previously line 4 of input.txt
"""
H5_SOLVER_GMRES_MAX_ITERATIONS = H5_INPUT_SOLVER + 'gmres_max_iterations'
"""
Maximum iterations for GMRES
Previously line 5 of input.txt
"""

H5_ENVIRONMENT = H5_CALCULATIONS + 'environment/'
"""
Group that contains the environment of the calculations.
Lines 2 to 5 of previous nemoh.cal
"""
H5_NUM_WAVE_FREQUENCIES = H5_ENVIRONMENT + 'num_wave_frequencies'
"""
The number of wave frequencies
Line 27, 1st number of previous nemoh.cal with only 1 body
"""
H5_MIN_WAVE_FREQUENCIES = H5_ENVIRONMENT + 'min_wave_frequencies'
"""
The minimum wave frequency, rad/s
Line 27, 2nd number of previous nemoh.cal with only 1 body
"""
H5_MAX_WAVE_FREQUENCIES = H5_ENVIRONMENT + 'max_wave_frequencies'

"""
The maximum wave frequency, rad/s
Line 27, 3rd number of previous nemoh.cal with only 1 body
"""

H5_NUM_WAVE_DIRECTIONS = H5_ENVIRONMENT + 'num_wave_directions'
"""
The number of wave directions
Line 28, 1st number of previous nemoh.cal with only 1 body
"""
H5_MIN_WAVE_DIRECTIONS = H5_ENVIRONMENT + 'min_wave_directions'
"""
The minimum wave direction, degree
Line 28, 2nd number of previous nemoh.cal with only 1 body
"""
H5_MAX_WAVE_DIRECTIONS = H5_ENVIRONMENT + 'max_wave_directions'
"""
The maximum wave direction, degree
Line 28, 3rd number of previous nemoh.cal with only 1 body
"""

H5_SHOW_PRESSURE = H5_ENVIRONMENT + 'show_pressure'
"""
Flag controlling whether or not to show pressure
Line 31 of previous nemoh.cal with only 1 body
"""

H5_KOCHIN_NUMBER = H5_ENVIRONMENT + 'kochin_number'
"""
Kochin, Number of directions of calculation (0 for no calculations)
Line 32, 1st number of previous nemoh.cal with only 1 body
"""
H5_KOCHIN_MIN = H5_ENVIRONMENT + 'kochin_min'
"""
Kochin, Minimum directions of calculation
Line 32, 2nd number of previous nemoh.cal with only 1 body
"""
H5_KOCHIN_MAX = H5_ENVIRONMENT + 'kochin_max'
"""
Kochin, Maximum directions of calculation
Line 32, 3rd number of previous nemoh.cal with only 1 body
"""

H5_FREE_SURFACE_POINTS_X = H5_ENVIRONMENT + 'free_surface_points_x'
"""
Free surface elevation, Number of points in x direction (0 for no calcutions)
Line 33, 1st number of previous nemoh.cal with only 1 body
"""
H5_FREE_SURFACE_POINTS_Y = H5_ENVIRONMENT + 'free_surface_points_y'
"""
Free surface elevation, Number of points in y direction (0 for no calcutions)
Line 33, 2nd number of previous nemoh.cal with only 1 body
"""
H5_FREE_SURFACE_DIMENSION_X = H5_ENVIRONMENT + 'free_surface_dimension_x'
"""
Free surface elevation, dimensions of domain in x direction
Line 33, 3rd number of previous nemoh.cal with only 1 body
"""
H5_FREE_SURFACE_DIMENSION_Y = H5_ENVIRONMENT + 'free_surface_dimension_y'
"""
Free surface elevation, dimensions of domain in y direction
Line 33, 4th number of previous nemoh.cal with only 1 body
"""

H5_ENV_VOLUME = H5_ENVIRONMENT + 'volume'
"""
Fluid specific volume (KG/M**3)
Line 2 of previous nemoh.cal
"""
H5_ENV_GRAVITY = H5_ENVIRONMENT + 'gravity'
"""
Gravity  (M/S**2)
Line 3 of previous nemoh.cal
"""
H5_ENV_DEPTH = H5_ENVIRONMENT + 'depth'
"""
Water depth (M)
Line 4 of previous nemoh.cal
"""


H5_ENV_WAVE_POINT = H5_ENVIRONMENT + 'wave_point'
"""
Wave Point
Line 5 of previous nemoh.cal
"""

H5_COMPUTE_IRF = H5_ENVIRONMENT + 'compute_irf'
"""
Flag controlling the irf computation. (0 for no calculation)
Line 30 1st number of previous nemoh.cal with only 1 body
"""
H5_IRF_TIME_STEP = H5_ENVIRONMENT + 'irf_time_step'
"""
IRF time step. (0 for no calculation)
Line 30 2nd number of previous nemoh.cal with only 1 body
"""
H5_IRF_DURATION = H5_ENVIRONMENT + 'irf_duration'
"""
IRF duration. (0 for no calculation)
Line 30 3rd number of previous nemoh.cal with only 1 body
"""

H5_BODY_MESH = 'mesh'
"""
Contains the mesh array of a body.
Content of the file reference in Line 9 of previous nemoh.cal.
Example Cylinder.dat
"""
H5_BODY_NUM_POINTS = 'env_num_points'
"""
The number of points of a body
Line 10 1st number of previous nemoh.cal
"""
H5_BODY_NUM_PANELS = 'env_num_panels'
"""
The number of panels of a body
Line 10 2nd number of previous nemoh.cal
"""


H5_BODIES = H5_CALCULATIONS + 'bodies/'
"""
Group to contain all the bodies
"""

body1 = H5_BODIES + 'body1/'
H5_BODY_BASE = 'body'
"""
The base group name to use for body inside the H5_BODIES group
"""

H5_FREEDOM_DEGREE = 'freedom_degree'
"""
Freedom degree of a body
"""

H5_GENERALISED_FORCES = 'generalised_forces'
"""
Generalised forces of a body
"""


#
#
#
#   Intermediate  OUTPUTS
#
#
#
H5_OUTPUT = 'output/'
"""
The hdf5 intermediate output group
"""
H5_MESH = H5_OUTPUT + 'mesh/'
"""
The group for mesh.
Contains values previously stored in Mesh
"""

H5_L12 = H5_MESH + 'l12/'
"""
The L12 group name
Contains values previously stored in Mesh/L12.dat
"""
H5_L12_X = H5_L12 + 'x'
"""
Nodes coordinates
"""
H5_L12_P = H5_L12 + 'p'
"""
Connectivities of the mesh
"""
H5_L12_COUNT = H5_L12 + 'params'
"""
The parameters for the mesh like the Symmetry about the xOz plane (1 for yes) (i_sym) variable
"""
H5_L10 = H5_MESH + 'l10/'
"""
The L10 group name
Contains values previously stored in Mesh/L10.dat
"""
H5_L10_CPANEL = H5_L10 + 'c_panel'
"""
To which body belongs the panel
"""
H5_L10_XM = H5_L10 + 'xm'
"""
Centre of panels
"""
H5_L10_N = H5_L10 + 'n'
"""
Normal vectors
"""
H5_L10_A = H5_L10 + 'a'
"""
Area of panel
"""
H5_L10_COUNT = H5_L10 + 'param'
"""
The parameters for the mesh like the number of points, of panels
"""
hdf5_tec = H5_MESH + 'mesh_tec'
H5_MESH_INTEGRATION = H5_MESH + 'integration'
"""
The integration results
"""
H5_MESH_FREE_SURFACE = H5_MESH + 'free_surface/'
"""
Free surface group name
"""
H5_MESH_FREE_SURFACE_VECTORS = H5_MESH_FREE_SURFACE + 'vectors'
"""
Free surface vectors
"""
H5_MESH_FREE_SURFACE_INDEX = H5_MESH_FREE_SURFACE + 'index'
"""
Free surface indices
"""
H5_MESH_KOCHIN = H5_MESH + 'kochin'
"""
Kochin values
"""


H5_NORMAL_VELOCITY = H5_OUTPUT + 'normal_velocity/'
"""
Group name for the normal velocity.
Contains values previously in NormalVelocities.dat
"""
H5_NORMAL_VELOCITY_W = H5_NORMAL_VELOCITY + "w"
"""
The wave frequency
"""
H5_NORMAL_VELOCITY_BETA = H5_NORMAL_VELOCITY + "beta"
"""
The wave directions
"""

H5_NORMAL_VELOCITY_BETA_RAW = H5_NORMAL_VELOCITY + "beta_raw"
"""
The wave directions
"""

H5_NORMAL_VELOCITY_SWITCH_POTENTIAL = H5_NORMAL_VELOCITY + "switch_potential"
"""
Array of flags controlling whether to show potential for each problem
"""
H5_NORMAL_VELOCITY_SWITCH_FREE_SURFACE = H5_NORMAL_VELOCITY + "switch_free_surface"
"""
Array of flags controlling whether to show free surface for each problem
"""
H5_NORMAL_VELOCITY_SWITCH_KOCHIN = H5_NORMAL_VELOCITY + "switch_kochin"
"""
Array of flags controlling whether to show kochin for each problem
"""
H5_NORMAL_VELOCITY_VELOCITIES = H5_NORMAL_VELOCITY + "velocities"
"""
The velocities array
"""


#
#
#
#   RESULTS
#
#
#
H5_RESULTS = 'results/'
"""
The group containing results.
Previously values inside Results/ directory
"""
H5_RESULTS_FK_FORCES = H5_RESULTS + 'fk_forces'
"""
The froude krylov forces
Mostly the value in previous results/FKForces.dat
"""
H5_RESULTS_CASE = H5_RESULTS + 'case/'
"""
The result case group
"""
H5_RESULTS_CASE_FORCE = H5_RESULTS_CASE + 'force'
"""
The case force
"""
H5_RESULTS_CASE_MOTION = H5_RESULTS_CASE + 'motion'
"""
The case motion
"""

H5_RESULTS_CASE_RADIATION = H5_RESULTS_CASE + 'radiation'
"""
The case radiation
"""

H5_RESULTS_CASE_BETA = H5_RESULTS_CASE + 'beta'
"""
The case wave directions
"""
H5_RESULTS_CASE_W = H5_RESULTS_CASE + 'w'
"""
The case wave frequencies
"""
H5_RESULTS_CASE_THETA = H5_RESULTS_CASE + 'theta'
"""
The case angle
"""
H5_RESULTS_FORCES = H5_RESULTS + 'forces'
"""
The forces.
Mostly the values in the previous results/Forces.dat
"""
H5_RESULTS_FK_FORCES_RAW = H5_RESULTS + 'fk_forces_raw'
"""
The raw froude krylov forces in complex number
"""
H5_RESULTS_KOCHIN = H5_RESULTS + 'kochin'
"""
The kochin number
Previously Kochin.*.dat
"""
H5_RESULTS_FREE_SURFACE_PANEL = H5_RESULTS + 'free_surface_panel'
"""
The free surface panel
Previously  results/free_surface*.dat  1st part
"""
H5_RESULTS_FREE_SURFACE_POINTS = H5_RESULTS + 'free_surface_points'
"""
The free surface points.
Previously  results/freesurface*.dat 2nd part
"""
H5_RESULTS_POTENTIAL = H5_RESULTS + 'potential'
"""
The potential
Previous results/potential*.dat
"""

H5_RESULTS_DRIFT_FORCES = H5_RESULTS + 'drift_forces'
"""
The mean drift forces
"""

H5_RESULTS_YAW_MOMENT = H5_RESULTS + 'yaw_moment'
"""
The mean yaw moment
"""

H5_RESULTS_ADDED_MASS = H5_RESULTS + 'added_mass'
"""
The added mass coefficients per wave frequency
"""

H5_RESULTS_ADDED_MASS_INFINITE = H5_RESULTS + 'added_mass_infinite'
"""
The infinite frequency added mass coefficients
"""

H5_RESULTS_ADDED_MASS_ZERO = H5_RESULTS + 'added_mass_zero'
"""
The zero frequency added mass coefficients
"""

H5_RESULTS_RADIATION_DAMPING = H5_RESULTS + 'radiation_damping'
"""
The radiation damping coefficients
"""

H5_RESULTS_EXCITATION_FORCES = H5_RESULTS + 'excitation_forces'
"""
The excitation forces coefficients
"""

H5_RESULTS_VOLUME_DISPLACEMENT = H5_RESULTS + 'volume_displacement'
"""
The volume displacement per body
"""

H5_RESULTS_CENTER_BUOYANCY = H5_RESULTS + 'center_buoyancy'
"""
The center of buoyancy per body
"""

H5_RESULTS_WATER_PLANE_AREA = H5_RESULTS + 'water_plane_area'
"""
The water plane area per body
"""

H5_RESULTS_STIFNESS = H5_RESULTS + 'stifness'
"""
The hydrostatic stifness matrix per body
"""






#
#
#
#   ATTRIBUTES CONFIGURATION
#   attributes should be a simple dictionary. the key must be a string
#   The value should of a type recognised by h5py: int types, float types or numpy string
#
#
#
#
#
#
#
#
#
#
#
#   INPUTS
#
#
#

H5_SOLVER_SWITCH_ODE_INFLUENCE_ATTR = {
    "description": np.string_("Indicate whether or not to use the ode method to compute the influence coefficients")
}
"""
Attribute for the H5_SOLVER_SWITCH_ODE_INFLUENCE path
"""

H5_SOLVER_USE_DIPOLES_IMPLEMENTATION_ATTR = {
    "description": np.string_("Whether or not to use dipoles Implementation. 1 for using it 0 for no.")
}
"""
Whether or not to use dipoles Implementation. 1 for using it 0 for no.
"""

H5_SOLVER_REMOVE_IRREGULAR_FREQUENCIES_ATTR = {
    "description": np.string_("Whether or not to remove irregular frequencies. 1 for using it 0 for no.")
}
"""
Whether or not to remove irregular frequencies. 1 for using it 0 for no
"""

H5_SOLVER_IS_INTERIOR_DOMAIN_ATTR = {
    "description": np.string_("Array indicating whether or not a panel is in the body or the interior free surface domain. 1 to indicate the interior, 0 for the body.")
}
"""
Array indicating whether or not a panel is in the body or the interior free surface domain. 1 to indicate the interior, 0 for the body.
"""



H5_SOLVER_COMPUTE_YAW_MOMENT_ATTR = {
    "description": np.string_("Whether or not to compute yaw moment. 1 for using it 0 for no.")
}
"""
Attribute for the H5_SOLVER_COMPUTE_YAW_MOMENT path
"""


H5_SOLVER_COMPUTE_DRIFT_FORCES_ATTR = {
    "description": np.string_("Whether or not to compute drift forces. 1 for using it 0 for no.")
}
"""
Attribute for the H5_SOLVER_COMPUTE_DRIFT_FORCES path
"""


H5_SOLVER_THIN_PANELS_ATTR = {
    "description": np.string_("Array containing whether a given panel is a conventional or a dipole one's." + 
        "1 for dipoles panel, 0 for conventional")
}
"""
Array containing whether a given panel is a conventional or a dipole one's.
1 for dipoles panel, 0 for conventional
"""

H5_SOLVER_USE_HIGHER_ORDER_ATTR = {
    "description": np.string_("Whether or not to use higher order panel method. 1 for using it 0 for no")
}
"""
Attribute for the H5_SOLVER_USE_HIGHER_ORDER path
"""

H5_SOLVER_NUM_PANEL_HIGHER_ORDER_ATTR = {
    "description": np.string_("The number of panel per patch in the higher order method")
}
"""
The number of panel per patch in the higher order method
"""

H5_SOLVER_B_SPLINE_ORDER_ATTR = {
    "description": np.string_("The order of the B-Spline for the potential in the higher order")
}
"""
The order of the B-Spline for the potential in the higher order
"""

H5_SOLVER_GREEN_TABULATION_NUMX_ATTR = {
    "description": np.string_("Number of points in x direction of tabulated data")
}

H5_SOLVER_GREEN_TABULATION_NUMZ_ATTR = {
    "description": np.string_("Number of points in z direction of tabulated data")
}

H5_SOLVER_GREEN_TABULATION_SIMPSON_NPOINTS_ATTR = {
    "description": np.string_("Number of sub intervals used to approximate the green function integral using simpson rule")
}


H5_SOLVER_TYPE_ATTR = {
    "description": np.string_("The solver type. (0) for Direct Gauss (1) for GMRES (2) GMRES with FMM acceleration (2 not implemented yet)")
}
"""
The solver type attributes. (0) for Direct Gauss (1) for GMRES (2) GMRES with FMM acceleration (2 not implemented yet)
Previously line 2 of input.txt
"""

H5_SOLVER_GMRES_RESTART_ATTR = {
    "description": np.string_("The Restart parameter for GMRES.")
}
"""
Attributes for The Restart parameter for GMRES.
Previously line 3 of input.txt
"""

H5_SOLVER_GMRES_STOPPING_ATTR = {
    "description": np.string_("Stopping criterion for GMRES")
}
"""
Attributes for the Stopping criterion for GMRES
Previously line 4 of input.txt
"""
H5_SOLVER_GMRES_MAX_ITERATIONS_ATTR = {
    "description": np.string_("Maximum iterations for GMRES")
}
"""
Attributes for the Maximum iterations for GMRES
Previously line 5 of input.txt
"""

H5_ENVIRONMENT_ATTR = {
    "description": np.string_("Group that contains the environment of the calculations.")
}
"""
Attributes for the Group that contains the environment of the calculations.
Lines 2 to 5 of previous nemoh.cal
"""
H5_NUM_WAVE_FREQUENCIES_ATTR = {
    "description": np.string_("The number of wave frequencies")
}
"""
The number of wave frequencies
Line 27, 1st number of previous nemoh.cal with only 1 body
"""
H5_MIN_WAVE_FREQUENCIES_ATTR = {
    "description": np.string_("The minimum wave frequency, rad/s")
}
"""
The minimum wave frequency, rad/s
Line 27, 2nd number of previous nemoh.cal with only 1 body
"""
H5_MAX_WAVE_FREQUENCIES_ATTR = {
    "description": np.string_("The maximum wave frequency, rad/s")
}


"""
The maximum wave frequency, rad/s
Line 27, 3rd number of previous nemoh.cal with only 1 body
"""

H5_NUM_WAVE_DIRECTIONS_ATTR = {
    "description": np.string_("The number of wave directions")
}
"""
The number of wave directions
Line 28, 1st number of previous nemoh.cal with only 1 body
"""
H5_MIN_WAVE_DIRECTIONS_ATTR = {
    "description": np.string_("The minimum wave direction, degree")
}
"""
The minimum wave direction, degree
Line 28, 2nd number of previous nemoh.cal with only 1 body
"""
H5_MAX_WAVE_DIRECTIONS_ATTR = {
    "description": np.string_("The maximum wave direction, degree")
}
"""
The maximum wave direction, degree
Line 28, 3rd number of previous nemoh.cal with only 1 body
"""

H5_SHOW_PRESSURE_ATTR = {
    "description": np.string_("Flag controlling whether or not to show pressure")
}
"""
Flag controlling whether or not to show pressure
Line 31 of previous nemoh.cal with only 1 body
"""

H5_KOCHIN_NUMBER_ATTR = {
    "description": np.string_("Kochin, Number of directions of calculation (0 for no calculations)")
}
"""
Kochin, Number of directions of calculation (0 for no calculations)
Line 32, 1st number of previous nemoh.cal with only 1 body
"""
H5_KOCHIN_MIN_ATTR = {
    "description": np.string_("Kochin, Minimum directions of calculation")
}
"""
Kochin, Minimum directions of calculation
Line 32, 2nd number of previous nemoh.cal with only 1 body
"""
H5_KOCHIN_MAX_ATTR = {
    "description": np.string_("Kochin, Maximum directions of calculation")
}
"""
Kochin, Maximum directions of calculation
Line 32, 3rd number of previous nemoh.cal with only 1 body
"""

H5_FREE_SURFACE_POINTS_X_ATTR = {
    "description": np.string_("The free surface elevation, Number of points in x direction (0 for no calcutions)")
}
"""
Free surface elevation, Number of points in x direction (0 for no calcutions)
Line 33, 1st number of previous nemoh.cal with only 1 body
"""
H5_FREE_SURFACE_POINTS_Y_ATTR = {
    "description": np.string_("The free surface elevation, Number of points in y direction (0 for no calcutions)")
}
"""
Free surface elevation, Number of points in y direction (0 for no calcutions)
Line 33, 2nd number of previous nemoh.cal with only 1 body
"""
H5_FREE_SURFACE_DIMENSION_X_ATTR = {
    "description": np.string_("The free surface elevation, dimensions of domain in x direction")
}
"""
Free surface elevation, dimensions of domain in x direction
Line 33, 3rd number of previous nemoh.cal with only 1 body
"""
H5_FREE_SURFACE_DIMENSION_Y_ATTR = {
    "description": np.string_("The free surface elevation, dimensions of domain in y direction")
}
"""
Free surface elevation, dimensions of domain in y direction
Line 33, 4th number of previous nemoh.cal with only 1 body
"""

H5_ENV_VOLUME_ATTR = {
    "description": np.string_("Fluid specific volume (KG/M**3)")
}
"""
Fluid specific volume (KG/M**3)
Line 2 of previous nemoh.cal
"""
H5_ENV_GRAVITY_ATTR = {
    "description": np.string_("Gravity  (M/S**2)")
}
"""
Gravity  (M/S**2)
Line 3 of previous nemoh.cal
"""
H5_ENV_DEPTH_ATTR = {
    "description": np.string_("Water depth (M)")
}
"""
Water depth (M)
Line 4 of previous nemoh.cal
"""


H5_ENV_WAVE_POINT_ATTR = {
    "description": np.string_("Wave Point")
}
"""
Wave Point
Line 5 of previous nemoh.cal
"""

H5_COMPUTE_IRF_ATTR = {
    "description": np.string_("Flag controlling the irf computation. (0 for no calculation)")
}
"""
Flag controlling the irf computation. (0 for no calculation)
Line 30 1st number of previous nemoh.cal with only 1 body
"""
H5_IRF_TIME_STEP_ATTR = {
    "description": np.string_("IRF time step. (0 for no calculation)")
}
"""
IRF time step. (0 for no calculation)
Line 30 2nd number of previous nemoh.cal with only 1 body
"""
H5_IRF_DURATION_ATTR = {
    "description": np.string_("IRF duration. (0 for no calculation)")
}
"""
IRF duration. (0 for no calculation)
Line 30 3rd number of previous nemoh.cal with only 1 body
"""

H5_BODY_MESH_ATTR = {
    "description": np.string_("Contains the mesh array of a body.")
}
"""
Contains the mesh array of a body.
Content of the file reference in Line 9 of previous nemoh.cal.
Example Cylinder.dat
"""
H5_BODY_NUM_POINTS_ATTR = {
    "description": np.string_("The number of points of a body")
}
"""
The number of points of a body
Line 10 1st number of previous nemoh.cal
"""
H5_BODY_NUM_PANELS_ATTR = {
    "description": np.string_("The number of panels of a body")
}
"""
The number of panels of a body
Line 10 2nd number of previous nemoh.cal
"""


H5_BODIES_ATTR = {
    "description": np.string_("Group to contain all the bodies")
}
"""
Group to contain all the bodies
"""


H5_BODY_BASE_ATTR = {
    "description": np.string_("The base group name to use for body inside the H5_BODIES group")
}
"""
The base group name to use for body inside the H5_BODIES group
"""

H5_FREEDOM_DEGREE_ATTR = {
    "description": np.string_("Freedom degree of a body")
}
"""
Freedom degree of a body
"""

H5_GENERALISED_FORCES_ATTR = {
    "description": np.string_("Generalised forces of a body")
}
"""
Generalised forces of a body
"""


#
#
#
#   Intermediate  OUTPUTS
#
#
#
H5_OUTPUT_ATTR = {
    "description": np.string_("The hdf5 intermediate output group")
}
"""
The hdf5 intermediate output group
"""
H5_MESH_ATTR = {
    "description": np.string_("The group for mesh.")
}
"""
The group for mesh.
Contains values previously stored in Mesh
"""

H5_L12_ATTR = {
    "description": np.string_("The L12 group name")
}
"""
The L12 group name
Contains values previously stored in Mesh/L12.dat
"""
H5_L12_X_ATTR = {
    "description": np.string_("Nodes coordinates")
}
"""
Nodes coordinates
"""
H5_L12_P_ATTR = {
    "description": np.string_("Connectivities of the mesh")
}
"""
Connectivities of the mesh
"""
H5_L12_COUNT_ATTR = {
    "description": np.string_("The parameters for the mesh like the Symmetry about the xOz plane (1 for yes) (i_sym) "
                              "variable")
}
"""
The parameters for the mesh like the Symmetry about the xOz plane (1 for yes) (i_sym) variable
"""
H5_L10_ATTR = {
    "description": np.string_("The L10 group name")
}
"""
The L10 group name
Contains values previously stored in Mesh/L10.dat
"""
H5_L10_CPANEL_ATTR = {
    "description": np.string_("To which body belongs the panel")
}
"""
To which body belongs the panel
"""
H5_L10_XM_ATTR = {
    "description": np.string_("Centre of panels")
}
"""
Centre of panels
"""
H5_L10_N_ATTR = {
    "description": np.string_("Normal vectors")
}
"""
Normal vectors
"""
H5_L10_A_ATTR = {
    "description": np.string_("Area of panel")
}
"""
Area of panel
"""
H5_L10_COUNT_ATTR = {
    "description": np.string_("The parameters for the mesh like the number of points, of panels")
}
"""
The parameters for the mesh like the number of points, of panels
"""
H5_MESH_INTEGRATION_ATTR = {
    "description": np.string_("The integration results")
}
"""
The integration results
"""
H5_MESH_FREE_SURFACE_ATTR = {
    "description": np.string_("Free surface group name")
}
"""
Free surface group name
"""
H5_MESH_FREE_SURFACE_VECTORS_ATTR = {
    "description": np.string_("Free surface vectors")
}
"""
Free surface vectors
"""
H5_MESH_FREE_SURFACE_INDEX_ATTR = {
    "description": np.string_("Free surface indices")
}
"""
Free surface indices
"""
H5_MESH_KOCHIN_ATTR = {
    "description": np.string_("Kochin values")
}
"""
Kochin values
"""


H5_NORMAL_VELOCITY_ATTR = {
    "description": np.string_("Group name for the normal velocity.")
}
"""
Group name for the normal velocity.
Contains values previously in NormalVelocities.dat
"""
H5_NORMAL_VELOCITY_W_ATTR = {
    "description": np.string_("The wave frequency")
}
"""
The wave frequency
"""
H5_NORMAL_VELOCITY_BETA_ATTR = {
    "description": np.string_("The wave directions")
}
"""
The wave directions
"""

H5_NORMAL_VELOCITY_BETA_RAW_ATTR = {
    "description": np.string_("The raw wave directions")
}
"""
The raw wave directions
"""

H5_NORMAL_VELOCITY_SWITCH_POTENTIAL_ATTR = {
    "description": np.string_("Array of flags controlling whether to show potential for each problem")
}
"""
Array of flags controlling whether to show potential for each problem
"""
H5_NORMAL_VELOCITY_SWITCH_FREE_SURFACE_ATTR = {
    "description": np.string_("Array of flags controlling whether to show free surface for each problem")
}
"""
Array of flags controlling whether to show free surface for each problem
"""
H5_NORMAL_VELOCITY_SWITCH_KOCHIN_ATTR = {
    "description": np.string_("Array of flags controlling whether to show kochin for each problem")
}
"""
Array of flags controlling whether to show kochin for each problem
"""
H5_NORMAL_VELOCITY_VELOCITIES_ATTR = {
    "description": np.string_("The velocities array")
}
"""
The velocities array
"""


#
#
#
#   RESULTS
#
#
#
H5_RESULTS_ATTR = {
    "description": np.string_("The group containing results.")
}
"""
The group containing results.
Previously values inside Results/ directory
"""
H5_RESULTS_FK_FORCES_ATTR = {
    "description": np.string_("The froude krylov forces")
}
"""
The froude krylov forces
Mostly the value in previous results/FKForces.dat
"""
H5_RESULTS_CASE_ATTR = {
    "description": np.string_("The result case group")
}
"""
The result case group
"""
H5_RESULTS_CASE_FORCE_ATTR = {
    "description": np.string_("The case force")
}
"""
The case force
"""
H5_RESULTS_CASE_MOTION_ATTR = {
    "description": np.string_("The case motion")
}
"""
The case motion
"""

H5_RESULTS_CASE_RADIATION_ATTR = {
    "description": np.string_("The case radiation")
}
"""
The case motion
"""

H5_RESULTS_CASE_BETA_ATTR = {
    "description": np.string_("The case wave directions")
}

"""
The case wave directions
"""
H5_RESULTS_CASE_W_ATTR = {
    "description": np.string_("The case wave frequencies")
}
"""
The case wave frequencies
"""
H5_RESULTS_CASE_THETA_ATTR = {
    "description": np.string_("The case angle")
}
"""
The case angle
"""
H5_RESULTS_FORCES_ATTR = {
    "description": np.string_("The forces.")
}
"""
The forces.
Mostly the values in the previous results/Forces.dat
"""
H5_RESULTS_FK_FORCES_RAW_ATTR = {
    "description": np.string_("The raw froude krylov forces in complex number")
}
"""
The raw froude krylov forces in complex number
"""
H5_RESULTS_KOCHIN_ATTR = {
    "description": np.string_("The kochin number")
}
"""
The kochin number
Previously Kochin.*.dat
"""
H5_RESULTS_FREE_SURFACE_PANEL_ATTR = {
    "description": np.string_("The free surface panel")
}
"""
The free surface panel
Previously  results/free_surface*.dat  1st part
"""
H5_RESULTS_FREE_SURFACE_POINTS_ATTR = {
    "description": np.string_("The free surface points.")
}
"""
The free surface points.
Previously  results/freesurface*.dat 2nd part
"""
H5_RESULTS_POTENTIAL_ATTR = {
    "description": np.string_("The potential")
}
"""
The potential
Previous results/potential*.dat
"""

H5_RESULTS_DRIFT_FORCES_ATTR = {
    "description": np.string_("The mean drift forces")
}
"""
The mean drift forces attributes
"""

H5_RESULTS_YAW_MOMENT_ATTR = {
    "description": np.string_("The mean yaw moment")
}
"""
The mean yaw moment
"""


H5_RESULTS_ADDED_MASS_ATTR = {
    "description": np.string_("The added mass coefficients per wave frequency")
}
"""
The added mass coefficients per wave frequency
"""

H5_RESULTS_ADDED_MASS_INFINITE_ATTR = {
    "description": np.string_("The infinite frequency added mass coefficients")
}
"""
The infinite frequency added mass coefficients
"""

H5_RESULTS_ADDED_MASS_ZERO_ATTR = {
    "description": np.string_("The zero frequency added mass coefficients")
}
"""
The zero frequency added mass coefficients
"""

H5_RESULTS_RADIATION_DAMPING_ATTR = {
    "description": np.string_("The radiation damping coefficients")
}
"""
The radiation damping coefficients
"""

H5_RESULTS_EXCITATION_FORCES_ATTR = {
    "description": np.string_("The excitation forces coefficients")
}
"""
The excitation forces coefficients
"""

H5_RESULTS_VOLUME_DISPLACEMENT_ATTR = {
    "description": np.string_("The volume displacement per body")
}
"""
The volume displacement per body
"""

H5_RESULTS_CENTER_BUOYANCY_ATTR = {
    "description": np.string_("The center of buoyancy per body")
}
"""
The center of buoyancy per body
"""

H5_RESULTS_WATER_PLANE_AREA_ATTR = {
    "description": np.string_("The water plane area per body")
}
"""
The water plane area per body
"""

H5_RESULTS_STIFNESS_ATTR = {
    "description": np.string_("The hydrostatic stifness matrix per body")
}
"""
The hydrostatic stifness matrix per body
"""

