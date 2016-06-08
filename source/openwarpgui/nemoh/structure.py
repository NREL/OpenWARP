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
       Added parameters to remove irregular frequencies or not: self.H5_SOLVER_REMOVE_IRREGULAR_FREQUENCIES
       and another controlling if a panel is in the body domain or interior free surface domain:
       self.H5_SOLVER_IS_INTERIOR_DOMAIN.

Changes in version 1.7 (OpenWarp - Add Logging Functionality)
       Added support for logging.

Changes in version 1.8 (OPENWARP - PROVIDE A COMMAND LINE INTERFACE USING PYTHON):
    Refactor this module to re-use the same structure as the hdf5 file.
"""

import numpy as np

__author__ = "yedtoss, TCSASSEMBLER"
__copyright__ = "Copyright (C) 2014-2016 TopCoder Inc. All rights reserved."
__version__ = "1.8"



class BaseStructure(object):
    """
    Base structure
    """


    def top_level_input_group(self):
        """
        Function to get the top level input group
        :return the top level input group
        """
        return ''

    def separator(self):
        """
        Function to get the path separator
        :return the path separator
        """
        return '/'

    def get_separator(self, include_separator):
        """
        Function to get the path separator depending on the flag
        :include_separator flag to tell whether or not to return the separator
        :return the path separator or empty string
        """
        if include_separator:
            return self.separator()
        else:
            return ''

    def environment(self):
        """
        Function to get environment name
        :return the environment name
        """
        return self.top_level_input_group() + 'environment' + self.separator()

    def floating_bodies(self, include_separator=True):
        """
        Function to get the floating_bodies name
        :return the floating_bodies name
        """
        return self.top_level_input_group() + 'floating_bodies' + self.get_separator(include_separator)

    def load_cases(self):
        """
        Function to get the load_cases name
        :return the load_cases name
        """
        return self.top_level_input_group() + 'load_cases' + self.separator()

    def post_processing(self):
        """
        Function to get the post_processing name
        :return the post_processing name
        """
        return self.top_level_input_group() + 'postprocessing' + self.separator()

    def calculation(self):
        """
        Function to get the calculation name
        :return the calculation name
        """
        return self.top_level_input_group() + 'calculation' + self.separator()

    def meshing(self):
        """
        Function to get the meshing name
        :return the path meshing name
        """
        return 'meshing' + self.separator()



    def __init__(self):
        """
        Initialize the members
        """

        #
        #
        #
        #   INPUTS
        #
        #
        #

        self.H5_INPUT_GROUP = self.top_level_input_group()
        """
        The hdf5 top level input group
        """

        self.LOAD_CASES_GROUP = self.load_cases()

        self.POST_PROCESSING_GROUP = self.post_processing()

        self.H5_CALCULATIONS = self.H5_INPUT_GROUP + 'calculations' + self.separator()
        """
        The hdf5 input calculations group. Contains mostly all the values previously in nemoh.cal file
        """

        self.H5_INPUT_SOLVER = self.H5_INPUT_GROUP + 'solver' + self.separator()
        """
        The hdf5 input solver group. Contains all values previously in input.txt
        """

        self.H5_SOLVER_USE_HIGHER_ORDER = self.calculation() + 'use_higher_order'
        """
        Whether or not to use higher order panel method. 1 for using it 0 for no
        """

        self.H5_SOLVER_REMOVE_IRREGULAR_FREQUENCIES = self.calculation() + 'remove_irregular_frequencies'
        """
        Whether or not to remove irregular frequencies. 1 for using it 0 for no
        """

        self.H5_SOLVER_USE_DIPOLES_IMPLEMENTATION = self.calculation() + 'use_dipoles_implementation'
        """
        Whether or not to use dipoles Implementation. 1 for using it 0 for no.
        """

        self.H5_SOLVER_COMPUTE_YAW_MOMENT = self.calculation() + 'compute_yaw_moment'
        """
        Whether or not to compute yaw moment. 1 for using it 0 for no.
        """

        self.H5_SOLVER_COMPUTE_DRIFT_FORCES = self.calculation() + 'compute_drift_forces'
        """
        Whether or not to compute drift forces. 1 for using it 0 for no.
        """


        self.H5_SOLVER_THIN_PANELS = self.calculation() + 'thin_panels'
        """
        Array containing whether a given panel is a conventional or a dipole one's.
        1 for dipoles panel, 0 for conventional
        """

        self.H5_SOLVER_IS_INTERIOR_DOMAIN = self.calculation() + 'IS_INTERIOR_DOMAIN'
        """
        Array indicating whether or not a panel is in the body or the interior free surface domain. 1 to indicate the interior, 0 for the body.
        """


        self.H5_SOLVER_NUM_PANEL_HIGHER_ORDER = self.calculation() + 'num_panel_higher_order'
        """
        The number of panel per patch in the higher order method
        """

        self.H5_SOLVER_B_SPLINE_ORDER = self.calculation() + 'b_spline_order'
        """
        The order of the B-Spline for the potential in the higher order
        """

        self.H5_SOLVER_SWITCH_ODE_INFLUENCE = self.calculation() + 'use_ode_influence_coefficients'
        """
        Indicate whether or not to use the ode method to compute the influence coefficients
        """

        self.H5_SOLVER_GREEN_TABULATION_NUMX = self.calculation() + 'green_tabulation_numx'
        """
        Number of points in x direction of tabulated data
        """

        self.H5_SOLVER_GREEN_TABULATION_NUMZ = self.calculation() + 'green_tabulation_numz'
        """
        Number of points in z direction of tabulated data
        """

        self.H5_SOLVER_GREEN_TABULATION_SIMPSON_NPOINTS = self.calculation() + 'green_tabulation_simpson_npoints'
        """
        Number of sub intervals used to approximate the green function integral
        using simpson rule
        """

        self.H5_SOLVER_TYPE = self.calculation() + 'solver_type'
        """
        The solver type. (0) for Direct Gauss (1) for GMRES (2) GMRES with FMM acceleration (2 not implemented yet)
        Previouly line 2 of input.txt
        """

        self.H5_SOLVER_GMRES_RESTART = self.calculation() + 'gmres_restart'
        """
        The Restart parameter for GMRES.
        Previously line 3 of input.txt
        """

        self.H5_SOLVER_GMRES_STOPPING = self.calculation() + 'gmres_stopping'
        """
        Stopping criterion for GMRES
        Previously line 4 of input.txt
        """
        self.H5_SOLVER_GMRES_MAX_ITERATIONS = self.calculation() + 'gmres_max_iterations'
        """
        Maximum iterations for GMRES
        Previously line 5 of input.txt
        """

        
        """
        Group that contains the environment of the calculations.
        Lines 2 to 5 of previous nemoh.cal
        """
        self.H5_NUM_WAVE_FREQUENCIES = self.load_cases() + 'num_wave_frequencies'
        """
        The number of wave frequencies
        Line 27, 1st number of previous nemoh.cal with only 1 body
        """
        self.H5_MIN_WAVE_FREQUENCIES = self.load_cases() + 'min_wave_frequencies'
        """
        The minimum wave frequency, rad/s
        Line 27, 2nd number of previous nemoh.cal with only 1 body
        """
        self.H5_MAX_WAVE_FREQUENCIES = self.load_cases() + 'max_wave_frequencies'

        """
        The maximum wave frequency, rad/s
        Line 27, 3rd number of previous nemoh.cal with only 1 body
        """

        self.H5_NUM_WAVE_DIRECTIONS = self.load_cases() + 'num_wave_directions'
        """
        The number of wave directions
        Line 28, 1st number of previous nemoh.cal with only 1 body
        """
        self.H5_MIN_WAVE_DIRECTIONS = self.load_cases() + 'min_wave_directions'
        """
        The minimum wave direction, degree
        Line 28, 2nd number of previous nemoh.cal with only 1 body
        """
        self.H5_MAX_WAVE_DIRECTIONS = self.load_cases() + 'max_wave_directions'
        """
        The maximum wave direction, degree
        Line 28, 3rd number of previous nemoh.cal with only 1 body
        """

        self.H5_SHOW_PRESSURE = self.post_processing() + 'show_pressure'
        """
        Flag controlling whether or not to show pressure
        Line 31 of previous nemoh.cal with only 1 body
        """

        self.H5_KOCHIN_NUMBER = self.post_processing() + 'kochin_num_directions'
        """
        Kochin, Number of directions of calculation (0 for no calculations)
        Line 32, 1st number of previous nemoh.cal with only 1 body
        """
        self.H5_KOCHIN_MIN = self.post_processing() + 'kochin_min_directions'
        """
        Kochin, Minimum directions of calculation
        Line 32, 2nd number of previous nemoh.cal with only 1 body
        """
        self.H5_KOCHIN_MAX = self.post_processing() + 'kochin_max_directions'
        """
        Kochin, Maximum directions of calculation
        Line 32, 3rd number of previous nemoh.cal with only 1 body
        """

        self.H5_FREE_SURFACE_POINTS_X = self.post_processing() + 'free_surface_num_points_x'
        """
        Free surface elevation, Number of points in x direction (0 for no calcutions)
        Line 33, 1st number of previous nemoh.cal with only 1 body
        """
        self.H5_FREE_SURFACE_POINTS_Y = self.post_processing() + 'free_surface_num_points_y'
        """
        Free surface elevation, Number of points in y direction (0 for no calcutions)
        Line 33, 2nd number of previous nemoh.cal with only 1 body
        """
        self.H5_FREE_SURFACE_DIMENSION_X = self.post_processing() + 'free_surface_dimension_x'
        """
        Free surface elevation, dimensions of domain in x direction
        Line 33, 3rd number of previous nemoh.cal with only 1 body
        """
        self.H5_FREE_SURFACE_DIMENSION_Y = self.post_processing() + 'free_surface_dimension_y'
        """
        Free surface elevation, dimensions of domain in y direction
        Line 33, 4th number of previous nemoh.cal with only 1 body
        """

        self.H5_ENV_VOLUME = self.environment() + 'fluid_mass_volumic'
        """
        Fluid specific volume (KG/M**3)
        Line 2 of previous nemoh.cal
        """
        self.H5_ENV_GRAVITY = self.environment() + 'gravity'
        """
        Gravity  (M/S**2)
        Line 3 of previous nemoh.cal
        """
        self.H5_ENV_DEPTH = self.environment() + 'fluid_depth'
        """
        Water depth (M)
        Line 4 of previous nemoh.cal
        """


        self.H5_ENV_WAVE_POINT = self.environment() + 'wave_measurement_point'
        """
        Wave Point
        Line 5 of previous nemoh.cal
        """

        self.H5_COMPUTE_IRF = self.post_processing() + 'compute_irf'
        """
        Flag controlling the irf computation. (0 for no calculation)
        Line 30 1st number of previous nemoh.cal with only 1 body
        """
        self.H5_IRF_TIME_STEP = self.post_processing() + 'irf_time_step'
        """
        IRF time step. (0 for no calculation)
        Line 30 2nd number of previous nemoh.cal with only 1 body
        """
        self.H5_IRF_DURATION = self.post_processing() + 'irf_duration'
        """
        IRF duration. (0 for no calculation)
        Line 30 3rd number of previous nemoh.cal with only 1 body
        """

        self.H5_BODY_MESH = 'mesh'
        """
        Contains the mesh array of a body.
        Content of the file reference in Line 9 of previous nemoh.cal.
        Example Cylinder.dat
        """
        self.H5_BODY_NUM_POINTS = 'env_num_points'
        """
        The number of points of a body
        Line 10 1st number of previous nemoh.cal
        """
        self.H5_BODY_NUM_PANELS = 'env_num_panels'
        """
        The number of panels of a body
        Line 10 2nd number of previous nemoh.cal
        """


        self.H5_BODIES = self.floating_bodies(False)
        """
        Group to contain all the bodies
        """

        self.body1 = self.H5_BODIES + 'body1' + self.separator()
        self.H5_BODY_BASE = 'body'
        """
        The base group name to use for body inside the self.H5_BODIES group
        """

        self.H5_FREEDOM_DEGREE = 'freedom_degree'
        """
        Freedom degree of a body
        """

        self.H5_GENERALISED_FORCES = 'generalised_forces'
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
        self.H5_OUTPUT = 'output' + self.separator()
        """
        The hdf5 intermediate output group
        """
        self.H5_MESH = self.H5_OUTPUT + 'mesh' + self.separator()
        """
        The group for mesh.
        Contains values previously stored in Mesh
        """

        self.H5_L12 = self.H5_MESH + 'l12' + self.separator()
        """
        The L12 group name
        Contains values previously stored in Mesh/L12.dat
        """
        self.H5_L12_X = self.H5_L12 + 'x'
        """
        Nodes coordinates
        """
        self.H5_L12_P = self.H5_L12 + 'p'
        """
        Connectivities of the mesh
        """
        self.H5_L12_COUNT = self.H5_L12 + 'params'
        """
        The parameters for the mesh like the Symmetry about the xOz plane (1 for yes) (i_sym) variable
        """
        self.H5_L10 = self.H5_MESH + 'l10' + self.separator()
        """
        The L10 group name
        Contains values previously stored in Mesh/L10.dat
        """
        self.H5_L10_CPANEL = self.H5_L10 + 'c_panel'
        """
        To which body belongs the panel
        """
        self.H5_L10_XM = self.H5_L10 + 'xm'
        """
        Centre of panels
        """
        self.H5_L10_N = self.H5_L10 + 'n'
        """
        Normal vectors
        """
        self.H5_L10_A = self.H5_L10 + 'a'
        """
        Area of panel
        """
        self.H5_L10_COUNT = self.H5_L10 + 'param'
        """
        The parameters for the mesh like the number of points, of panels
        """
        self.hdf5_tec = self.H5_MESH + 'mesh_tec'
        self.H5_MESH_INTEGRATION = self.H5_MESH + 'integration'
        """
        The integration results
        """
        self.H5_MESH_FREE_SURFACE = self.H5_MESH + 'free_surface' + self.separator()
        """
        Free surface group name
        """
        self.H5_MESH_FREE_SURFACE_VECTORS = self.H5_MESH_FREE_SURFACE + 'vectors'
        """
        Free surface vectors
        """
        self.H5_MESH_FREE_SURFACE_INDEX = self.H5_MESH_FREE_SURFACE + 'index'
        """
        Free surface indices
        """
        self.H5_MESH_KOCHIN = self.H5_MESH + 'kochin'
        """
        Kochin values
        """


        self.H5_NORMAL_VELOCITY = self.H5_OUTPUT + 'normal_velocity' + self.separator()
        """
        Group name for the normal velocity.
        Contains values previously in NormalVelocities.dat
        """
        self.H5_NORMAL_VELOCITY_W = self.H5_NORMAL_VELOCITY + "w"
        """
        The wave frequency
        """
        self.H5_NORMAL_VELOCITY_BETA = self.H5_NORMAL_VELOCITY + "beta"
        """
        The wave directions
        """

        self.H5_NORMAL_VELOCITY_BETA_RAW = self.H5_NORMAL_VELOCITY + "beta_raw"
        """
        The wave directions
        """

        self.H5_NORMAL_VELOCITY_SWITCH_POTENTIAL = self.H5_NORMAL_VELOCITY + "switch_potential"
        """
        Array of flags controlling whether to show potential for each problem
        """
        self.H5_NORMAL_VELOCITY_SWITCH_FREE_SURFACE = self.H5_NORMAL_VELOCITY + "switch_free_surface"
        """
        Array of flags controlling whether to show free surface for each problem
        """
        self.H5_NORMAL_VELOCITY_SWITCH_KOCHIN = self.H5_NORMAL_VELOCITY + "switch_kochin"
        """
        Array of flags controlling whether to show kochin for each problem
        """
        self.H5_NORMAL_VELOCITY_VELOCITIES = self.H5_NORMAL_VELOCITY + "velocities"
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
        self.H5_RESULTS = 'results' + self.separator()
        """
        The group containing results.
        Previously values inside Results/ directory
        """
        self.H5_RESULTS_FK_FORCES = self.H5_RESULTS + 'fk_forces'
        """
        The froude krylov forces
        Mostly the value in previous results/FKForces.dat
        """
        self.H5_RESULTS_CASE = self.H5_RESULTS + 'case' + self.separator()
        """
        The result case group
        """
        self.H5_RESULTS_CASE_FORCE = self.H5_RESULTS_CASE + 'force'
        """
        The case force
        """
        self.H5_RESULTS_CASE_MOTION = self.H5_RESULTS_CASE + 'motion'
        """
        The case motion
        """

        self.H5_RESULTS_CASE_RADIATION = self.H5_RESULTS_CASE + 'radiation'
        """
        The case radiation
        """

        self.H5_RESULTS_CASE_BETA = self.H5_RESULTS_CASE + 'beta'
        """
        The case wave directions
        """
        self.H5_RESULTS_CASE_W = self.H5_RESULTS_CASE + 'w'
        """
        The case wave frequencies
        """
        self.H5_RESULTS_CASE_THETA = self.H5_RESULTS_CASE + 'theta'
        """
        The case angle
        """
        self.H5_RESULTS_FORCES = self.H5_RESULTS + 'forces'
        """
        The forces.
        Mostly the values in the previous results/Forces.dat
        """
        self.H5_RESULTS_FK_FORCES_RAW = self.H5_RESULTS + 'fk_forces_raw'
        """
        The raw froude krylov forces in complex number
        """
        self.H5_RESULTS_KOCHIN = self.H5_RESULTS + 'kochin'
        """
        The kochin number
        Previously Kochin.*.dat
        """
        self.H5_RESULTS_FREE_SURFACE_PANEL = self.H5_RESULTS + 'free_surface_panel'
        """
        The free surface panel
        Previously  results/free_surface*.dat  1st part
        """
        self.H5_RESULTS_FREE_SURFACE_POINTS = self.H5_RESULTS + 'free_surface_points'
        """
        The free surface points.
        Previously  results/freesurface*.dat 2nd part
        """
        self.H5_RESULTS_POTENTIAL = self.H5_RESULTS + 'potential'
        """
        The potential
        Previous results/potential*.dat
        """

        self.H5_RESULTS_DRIFT_FORCES = self.H5_RESULTS + 'drift_forces'
        """
        The mean drift forces
        """

        self.H5_RESULTS_YAW_MOMENT = self.H5_RESULTS + 'yaw_moment'
        """
        The mean yaw moment
        """

        self.H5_RESULTS_ADDED_MASS = self.H5_RESULTS + 'added_mass'
        """
        The added mass coefficients per wave frequency
        """

        self.H5_RESULTS_ADDED_MASS_INFINITE = self.H5_RESULTS + 'added_mass_infinite'
        """
        The infinite frequency added mass coefficients
        """

        self.H5_RESULTS_ADDED_MASS_ZERO = self.H5_RESULTS + 'added_mass_zero'
        """
        The zero frequency added mass coefficients
        """

        self.H5_RESULTS_RADIATION_DAMPING = self.H5_RESULTS + 'radiation_damping'
        """
        The radiation damping coefficients
        """

        self.H5_RESULTS_EXCITATION_FORCES = self.H5_RESULTS + 'excitation_forces'
        """
        The excitation forces coefficients
        """

        self.H5_RESULTS_VOLUME_DISPLACEMENT = self.H5_RESULTS + 'volume_displacement'
        """
        The volume displacement per body
        """

        self.H5_RESULTS_CENTER_BUOYANCY = self.H5_RESULTS + 'center_buoyancy'
        """
        The center of buoyancy per body
        """

        self.H5_RESULTS_WATER_PLANE_AREA = self.H5_RESULTS + 'water_plane_area'
        """
        The water plane area per body
        """

        self.H5_RESULTS_STIFNESS = self.H5_RESULTS + 'stifness'
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

        self.H5_SOLVER_SWITCH_ODE_INFLUENCE_ATTR = {
            "description": np.string_("Indicate whether or not to use the ode method to compute the influence coefficients")
        }
        """
        Attribute for the self.H5_SOLVER_SWITCH_ODE_INFLUENCE path
        """

        self.H5_SOLVER_USE_DIPOLES_IMPLEMENTATION_ATTR = {
            "description": np.string_("Whether or not to use dipoles Implementation. 1 for using it 0 for no.")
        }
        """
        Whether or not to use dipoles Implementation. 1 for using it 0 for no.
        """

        self.H5_SOLVER_REMOVE_IRREGULAR_FREQUENCIES_ATTR = {
            "description": np.string_("Whether or not to remove irregular frequencies. 1 for using it 0 for no.")
        }
        """
        Whether or not to remove irregular frequencies. 1 for using it 0 for no
        """

        self.H5_SOLVER_IS_INTERIOR_DOMAIN_ATTR = {
            "description": np.string_("Array indicating whether or not a panel is in the body or the interior free surface domain. 1 to indicate the interior, 0 for the body.")
        }
        """
        Array indicating whether or not a panel is in the body or the interior free surface domain. 1 to indicate the interior, 0 for the body.
        """



        self.H5_SOLVER_COMPUTE_YAW_MOMENT_ATTR = {
            "description": np.string_("Whether or not to compute yaw moment. 1 for using it 0 for no.")
        }
        """
        Attribute for the self.H5_SOLVER_COMPUTE_YAW_MOMENT path
        """


        self.H5_SOLVER_COMPUTE_DRIFT_FORCES_ATTR = {
            "description": np.string_("Whether or not to compute drift forces. 1 for using it 0 for no.")
        }
        """
        Attribute for the self.H5_SOLVER_COMPUTE_DRIFT_FORCES path
        """


        self.H5_SOLVER_THIN_PANELS_ATTR = {
            "description": np.string_("Array containing whether a given panel is a conventional or a dipole one's." + 
                "1 for dipoles panel, 0 for conventional")
        }
        """
        Array containing whether a given panel is a conventional or a dipole one's.
        1 for dipoles panel, 0 for conventional
        """

        self.H5_SOLVER_USE_HIGHER_ORDER_ATTR = {
            "description": np.string_("Whether or not to use higher order panel method. 1 for using it 0 for no")
        }
        """
        Attribute for the self.H5_SOLVER_USE_HIGHER_ORDER path
        """

        self.H5_SOLVER_NUM_PANEL_HIGHER_ORDER_ATTR = {
            "description": np.string_("The number of panel per patch in the higher order method")
        }
        """
        The number of panel per patch in the higher order method
        """

        self.H5_SOLVER_B_SPLINE_ORDER_ATTR = {
            "description": np.string_("The order of the B-Spline for the potential in the higher order")
        }
        """
        The order of the B-Spline for the potential in the higher order
        """

        self.H5_SOLVER_GREEN_TABULATION_NUMX_ATTR = {
            "description": np.string_("Number of points in x direction of tabulated data")
        }

        self.H5_SOLVER_GREEN_TABULATION_NUMZ_ATTR = {
            "description": np.string_("Number of points in z direction of tabulated data")
        }

        self.H5_SOLVER_GREEN_TABULATION_SIMPSON_NPOINTS_ATTR = {
            "description": np.string_("Number of sub intervals used to approximate the green function integral using simpson rule")
        }


        self.H5_SOLVER_TYPE_ATTR = {
            "description": np.string_("The solver type. (0) for Direct Gauss (1) for GMRES (2) GMRES with FMM acceleration (2 not implemented yet)")
        }
        """
        The solver type attributes. (0) for Direct Gauss (1) for GMRES (2) GMRES with FMM acceleration (2 not implemented yet)
        Previously line 2 of input.txt
        """

        self.H5_SOLVER_GMRES_RESTART_ATTR = {
            "description": np.string_("The Restart parameter for GMRES.")
        }
        """
        Attributes for The Restart parameter for GMRES.
        Previously line 3 of input.txt
        """

        self.H5_SOLVER_GMRES_STOPPING_ATTR = {
            "description": np.string_("Stopping criterion for GMRES")
        }
        """
        Attributes for the Stopping criterion for GMRES
        Previously line 4 of input.txt
        """
        self.H5_SOLVER_GMRES_MAX_ITERATIONS_ATTR = {
            "description": np.string_("Maximum iterations for GMRES")
        }
        """
        Attributes for the Maximum iterations for GMRES
        Previously line 5 of input.txt
        """

        self.H5_ENVIRONMENT_ATTR = {
            "description": np.string_("Group that contains the environment of the calculations.")
        }
        """
        Attributes for the Group that contains the environment of the calculations.
        Lines 2 to 5 of previous nemoh.cal
        """
        self.H5_NUM_WAVE_FREQUENCIES_ATTR = {
            "description": np.string_("The number of wave frequencies")
        }
        """
        The number of wave frequencies
        Line 27, 1st number of previous nemoh.cal with only 1 body
        """
        self.H5_MIN_WAVE_FREQUENCIES_ATTR = {
            "description": np.string_("The minimum wave frequency, rad/s")
        }
        """
        The minimum wave frequency, rad/s
        Line 27, 2nd number of previous nemoh.cal with only 1 body
        """
        self.H5_MAX_WAVE_FREQUENCIES_ATTR = {
            "description": np.string_("The maximum wave frequency, rad/s")
        }


        """
        The maximum wave frequency, rad/s
        Line 27, 3rd number of previous nemoh.cal with only 1 body
        """

        self.H5_NUM_WAVE_DIRECTIONS_ATTR = {
            "description": np.string_("The number of wave directions")
        }
        """
        The number of wave directions
        Line 28, 1st number of previous nemoh.cal with only 1 body
        """
        self.H5_MIN_WAVE_DIRECTIONS_ATTR = {
            "description": np.string_("The minimum wave direction, degree")
        }
        """
        The minimum wave direction, degree
        Line 28, 2nd number of previous nemoh.cal with only 1 body
        """
        self.H5_MAX_WAVE_DIRECTIONS_ATTR = {
            "description": np.string_("The maximum wave direction, degree")
        }
        """
        The maximum wave direction, degree
        Line 28, 3rd number of previous nemoh.cal with only 1 body
        """

        self.H5_SHOW_PRESSURE_ATTR = {
            "description": np.string_("Flag controlling whether or not to show pressure")
        }
        """
        Flag controlling whether or not to show pressure
        Line 31 of previous nemoh.cal with only 1 body
        """

        self.H5_KOCHIN_NUMBER_ATTR = {
            "description": np.string_("Kochin, Number of directions of calculation (0 for no calculations)")
        }
        """
        Kochin, Number of directions of calculation (0 for no calculations)
        Line 32, 1st number of previous nemoh.cal with only 1 body
        """
        self.H5_KOCHIN_MIN_ATTR = {
            "description": np.string_("Kochin, Minimum directions of calculation")
        }
        """
        Kochin, Minimum directions of calculation
        Line 32, 2nd number of previous nemoh.cal with only 1 body
        """
        self.H5_KOCHIN_MAX_ATTR = {
            "description": np.string_("Kochin, Maximum directions of calculation")
        }
        """
        Kochin, Maximum directions of calculation
        Line 32, 3rd number of previous nemoh.cal with only 1 body
        """

        self.H5_FREE_SURFACE_POINTS_X_ATTR = {
            "description": np.string_("The free surface elevation, Number of points in x direction (0 for no calcutions)")
        }
        """
        Free surface elevation, Number of points in x direction (0 for no calcutions)
        Line 33, 1st number of previous nemoh.cal with only 1 body
        """
        self.H5_FREE_SURFACE_POINTS_Y_ATTR = {
            "description": np.string_("The free surface elevation, Number of points in y direction (0 for no calcutions)")
        }
        """
        Free surface elevation, Number of points in y direction (0 for no calcutions)
        Line 33, 2nd number of previous nemoh.cal with only 1 body
        """
        self.H5_FREE_SURFACE_DIMENSION_X_ATTR = {
            "description": np.string_("The free surface elevation, dimensions of domain in x direction")
        }
        """
        Free surface elevation, dimensions of domain in x direction
        Line 33, 3rd number of previous nemoh.cal with only 1 body
        """
        self.H5_FREE_SURFACE_DIMENSION_Y_ATTR = {
            "description": np.string_("The free surface elevation, dimensions of domain in y direction")
        }
        """
        Free surface elevation, dimensions of domain in y direction
        Line 33, 4th number of previous nemoh.cal with only 1 body
        """

        self.H5_ENV_VOLUME_ATTR = {
            "description": np.string_("Fluid specific volume (KG/M**3)")
        }
        """
        Fluid specific volume (KG/M**3)
        Line 2 of previous nemoh.cal
        """
        self.H5_ENV_GRAVITY_ATTR = {
            "description": np.string_("Gravity  (M/S**2)")
        }
        """
        Gravity  (M/S**2)
        Line 3 of previous nemoh.cal
        """
        self.H5_ENV_DEPTH_ATTR = {
            "description": np.string_("Water depth (M)")
        }
        """
        Water depth (M)
        Line 4 of previous nemoh.cal
        """


        self.H5_ENV_WAVE_POINT_ATTR = {
            "description": np.string_("Wave Point")
        }
        """
        Wave Point
        Line 5 of previous nemoh.cal
        """

        self.H5_COMPUTE_IRF_ATTR = {
            "description": np.string_("Flag controlling the irf computation. (0 for no calculation)")
        }
        """
        Flag controlling the irf computation. (0 for no calculation)
        Line 30 1st number of previous nemoh.cal with only 1 body
        """
        self.H5_IRF_TIME_STEP_ATTR = {
            "description": np.string_("IRF time step. (0 for no calculation)")
        }
        """
        IRF time step. (0 for no calculation)
        Line 30 2nd number of previous nemoh.cal with only 1 body
        """
        self.H5_IRF_DURATION_ATTR = {
            "description": np.string_("IRF duration. (0 for no calculation)")
        }
        """
        IRF duration. (0 for no calculation)
        Line 30 3rd number of previous nemoh.cal with only 1 body
        """

        self.H5_BODY_MESH_ATTR = {
            "description": np.string_("Contains the mesh array of a body.")
        }
        """
        Contains the mesh array of a body.
        Content of the file reference in Line 9 of previous nemoh.cal.
        Example Cylinder.dat
        """
        self.H5_BODY_NUM_POINTS_ATTR = {
            "description": np.string_("The number of points of a body")
        }
        """
        The number of points of a body
        Line 10 1st number of previous nemoh.cal
        """
        self.H5_BODY_NUM_PANELS_ATTR = {
            "description": np.string_("The number of panels of a body")
        }
        """
        The number of panels of a body
        Line 10 2nd number of previous nemoh.cal
        """


        self.H5_BODIES_ATTR = {
            "description": np.string_("Group to contain all the bodies")
        }
        """
        Group to contain all the bodies
        """


        self.H5_BODY_BASE_ATTR = {
            "description": np.string_("The base group name to use for body inside the self.H5_BODIES group")
        }
        """
        The base group name to use for body inside the self.H5_BODIES group
        """

        self.H5_FREEDOM_DEGREE_ATTR = {
            "description": np.string_("Freedom degree of a body")
        }
        """
        Freedom degree of a body
        """

        self.H5_GENERALISED_FORCES_ATTR = {
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
        self.H5_OUTPUT_ATTR = {
            "description": np.string_("The hdf5 intermediate output group")
        }
        """
        The hdf5 intermediate output group
        """
        self.H5_MESH_ATTR = {
            "description": np.string_("The group for mesh.")
        }
        """
        The group for mesh.
        Contains values previously stored in Mesh
        """

        self.H5_L12_ATTR = {
            "description": np.string_("The L12 group name")
        }
        """
        The L12 group name
        Contains values previously stored in Mesh/L12.dat
        """
        self.H5_L12_X_ATTR = {
            "description": np.string_("Nodes coordinates")
        }
        """
        Nodes coordinates
        """
        self.H5_L12_P_ATTR = {
            "description": np.string_("Connectivities of the mesh")
        }
        """
        Connectivities of the mesh
        """
        self.H5_L12_COUNT_ATTR = {
            "description": np.string_("The parameters for the mesh like the Symmetry about the xOz plane (1 for yes) (i_sym) "
                                      "variable")
        }
        """
        The parameters for the mesh like the Symmetry about the xOz plane (1 for yes) (i_sym) variable
        """
        self.H5_L10_ATTR = {
            "description": np.string_("The L10 group name")
        }
        """
        The L10 group name
        Contains values previously stored in Mesh/L10.dat
        """
        self.H5_L10_CPANEL_ATTR = {
            "description": np.string_("To which body belongs the panel")
        }
        """
        To which body belongs the panel
        """
        self.H5_L10_XM_ATTR = {
            "description": np.string_("Centre of panels")
        }
        """
        Centre of panels
        """
        self.H5_L10_N_ATTR = {
            "description": np.string_("Normal vectors")
        }
        """
        Normal vectors
        """
        self.H5_L10_A_ATTR = {
            "description": np.string_("Area of panel")
        }
        """
        Area of panel
        """
        self.H5_L10_COUNT_ATTR = {
            "description": np.string_("The parameters for the mesh like the number of points, of panels")
        }
        """
        The parameters for the mesh like the number of points, of panels
        """
        self.H5_MESH_INTEGRATION_ATTR = {
            "description": np.string_("The integration results")
        }
        """
        The integration results
        """
        self.H5_MESH_FREE_SURFACE_ATTR = {
            "description": np.string_("Free surface group name")
        }
        """
        Free surface group name
        """
        self.H5_MESH_FREE_SURFACE_VECTORS_ATTR = {
            "description": np.string_("Free surface vectors")
        }
        """
        Free surface vectors
        """
        self.H5_MESH_FREE_SURFACE_INDEX_ATTR = {
            "description": np.string_("Free surface indices")
        }
        """
        Free surface indices
        """
        self.H5_MESH_KOCHIN_ATTR = {
            "description": np.string_("Kochin values")
        }
        """
        Kochin values
        """


        self.H5_NORMAL_VELOCITY_ATTR = {
            "description": np.string_("Group name for the normal velocity.")
        }
        """
        Group name for the normal velocity.
        Contains values previously in NormalVelocities.dat
        """
        self.H5_NORMAL_VELOCITY_W_ATTR = {
            "description": np.string_("The wave frequency")
        }
        """
        The wave frequency
        """
        self.H5_NORMAL_VELOCITY_BETA_ATTR = {
            "description": np.string_("The wave directions")
        }
        """
        The wave directions
        """

        self.H5_NORMAL_VELOCITY_BETA_RAW_ATTR = {
            "description": np.string_("The raw wave directions")
        }
        """
        The raw wave directions
        """

        self.H5_NORMAL_VELOCITY_SWITCH_POTENTIAL_ATTR = {
            "description": np.string_("Array of flags controlling whether to show potential for each problem")
        }
        """
        Array of flags controlling whether to show potential for each problem
        """
        self.H5_NORMAL_VELOCITY_SWITCH_FREE_SURFACE_ATTR = {
            "description": np.string_("Array of flags controlling whether to show free surface for each problem")
        }
        """
        Array of flags controlling whether to show free surface for each problem
        """
        self.H5_NORMAL_VELOCITY_SWITCH_KOCHIN_ATTR = {
            "description": np.string_("Array of flags controlling whether to show kochin for each problem")
        }
        """
        Array of flags controlling whether to show kochin for each problem
        """
        self.H5_NORMAL_VELOCITY_VELOCITIES_ATTR = {
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
        self.H5_RESULTS_ATTR = {
            "description": np.string_("The group containing results.")
        }
        """
        The group containing results.
        Previously values inside Results/ directory
        """
        self.H5_RESULTS_FK_FORCES_ATTR = {
            "description": np.string_("The froude krylov forces")
        }
        """
        The froude krylov forces
        Mostly the value in previous results/FKForces.dat
        """
        self.H5_RESULTS_CASE_ATTR = {
            "description": np.string_("The result case group")
        }
        """
        The result case group
        """
        self.H5_RESULTS_CASE_FORCE_ATTR = {
            "description": np.string_("The case force")
        }
        """
        The case force
        """
        self.H5_RESULTS_CASE_MOTION_ATTR = {
            "description": np.string_("The case motion")
        }
        """
        The case motion
        """

        self.H5_RESULTS_CASE_RADIATION_ATTR = {
            "description": np.string_("The case radiation")
        }
        """
        The case motion
        """

        self.H5_RESULTS_CASE_BETA_ATTR = {
            "description": np.string_("The case wave directions")
        }

        """
        The case wave directions
        """
        self.H5_RESULTS_CASE_W_ATTR = {
            "description": np.string_("The case wave frequencies")
        }
        """
        The case wave frequencies
        """
        self.H5_RESULTS_CASE_THETA_ATTR = {
            "description": np.string_("The case angle")
        }
        """
        The case angle
        """
        self.H5_RESULTS_FORCES_ATTR = {
            "description": np.string_("The forces.")
        }
        """
        The forces.
        Mostly the values in the previous results/Forces.dat
        """
        self.H5_RESULTS_FK_FORCES_RAW_ATTR = {
            "description": np.string_("The raw froude krylov forces in complex number")
        }
        """
        The raw froude krylov forces in complex number
        """
        self.H5_RESULTS_KOCHIN_ATTR = {
            "description": np.string_("The kochin number")
        }
        """
        The kochin number
        Previously Kochin.*.dat
        """
        self.H5_RESULTS_FREE_SURFACE_PANEL_ATTR = {
            "description": np.string_("The free surface panel")
        }
        """
        The free surface panel
        Previously  results/free_surface*.dat  1st part
        """
        self.H5_RESULTS_FREE_SURFACE_POINTS_ATTR = {
            "description": np.string_("The free surface points.")
        }
        """
        The free surface points.
        Previously  results/freesurface*.dat 2nd part
        """
        self.H5_RESULTS_POTENTIAL_ATTR = {
            "description": np.string_("The potential")
        }
        """
        The potential
        Previous results/potential*.dat
        """

        self.H5_RESULTS_DRIFT_FORCES_ATTR = {
            "description": np.string_("The mean drift forces")
        }
        """
        The mean drift forces attributes
        """

        self.H5_RESULTS_YAW_MOMENT_ATTR = {
            "description": np.string_("The mean yaw moment")
        }
        """
        The mean yaw moment
        """


        self.H5_RESULTS_ADDED_MASS_ATTR = {
            "description": np.string_("The added mass coefficients per wave frequency")
        }
        """
        The added mass coefficients per wave frequency
        """

        self.H5_RESULTS_ADDED_MASS_INFINITE_ATTR = {
            "description": np.string_("The infinite frequency added mass coefficients")
        }
        """
        The infinite frequency added mass coefficients
        """

        self.H5_RESULTS_ADDED_MASS_ZERO_ATTR = {
            "description": np.string_("The zero frequency added mass coefficients")
        }
        """
        The zero frequency added mass coefficients
        """

        self.H5_RESULTS_RADIATION_DAMPING_ATTR = {
            "description": np.string_("The radiation damping coefficients")
        }
        """
        The radiation damping coefficients
        """

        self.H5_RESULTS_EXCITATION_FORCES_ATTR = {
            "description": np.string_("The excitation forces coefficients")
        }
        """
        The excitation forces coefficients
        """

        self.H5_RESULTS_VOLUME_DISPLACEMENT_ATTR = {
            "description": np.string_("The volume displacement per body")
        }
        """
        The volume displacement per body
        """

        self.H5_RESULTS_CENTER_BUOYANCY_ATTR = {
            "description": np.string_("The center of buoyancy per body")
        }
        """
        The center of buoyancy per body
        """

        self.H5_RESULTS_WATER_PLANE_AREA_ATTR = {
            "description": np.string_("The water plane area per body")
        }
        """
        The water plane area per body
        """

        self.H5_RESULTS_STIFNESS_ATTR = {
            "description": np.string_("The hydrostatic stifness matrix per body")
        }
        """
        The hydrostatic stifness matrix per body
        """


class H5_STRUCTURE(BaseStructure):
    """
    This defines the hdf5 structure
    """
    def top_level_input_group(self):
        """
        Function to get the top level input group
        :return the top level input group
        """
        return 'input' + self.separator()

    def environment(self):
        """
        Function to get environment name
        :return the environment name
        """
        return self.top_level_input_group() + 'calculations' + self.separator() + 'environment' + self.separator()

    def floating_bodies(self, include_separator=True):
        """
        Function to get floating_bodies name
        :return the floating_bodies name
        """
        return self.top_level_input_group() + 'calculations' + self.separator() + 'bodies' + self.separator()

    def load_cases(self):
        """
        Function to get load_cases name
        :return the load_cases name
        """
        return self.environment()

    def post_processing(self):
        """
        Function to get post_processingname
        :return the post_processing name
        """
        return self.environment()

    def calculation(self):
        """
        Function to get calculation name
        :return the calculation name
        """
        return self.top_level_input_group() + 'solver' + self.separator()      



class JSON_STRUCTURE(BaseStructure):
    """
    This define the json structure
    """
    def separator(self):
        """
        Function to get the path separator
        :return the path separator
        """
        return '.'

    def __init__(self):
        super(self.__class__, self).__init__()

        self.BODY_MESHING_FINENESS = self.meshing() + 'fineness'
        self.BODY_MESHING_GRADING = self.meshing() + 'grading'
        self.BODY_MESHING_MAXH = self.meshing() + 'maxh'
        self.BODY_MESHING_MINH = self.meshing() + 'minh'
        self.BODY_MESHING_USE_TOLERANCE = self.meshing() + 'use_tolerance'
        self.BODY_MESHING_TOLERANCE = self.meshing() + 'tolerance'
        self.BODY_SURGE = 'surge'
        self.BODY_SWAY = 'sway'
        self.BODY_HEAVE = 'heave'
        self.BODY_ROLL = 'roll'
        self.BODY_PITCH = 'pitch'
        self.BODY_YAW = 'yaw'

        self.BODY_FORCE_X = 'force_in_x_direction'
        self.BODY_FORCE_Y = 'force_in_y_direction'
        self.BODY_FORCE_Z = 'force_in_z_direction'
        self.BODY_MOMENT_X = 'moment_in_x_direction'
        self.BODY_MOMENT_Y = 'moment_in_y_direction'
        self.BODY_MOMENT_Z = 'moment_in_z_direction'
        self.BODY_MESH_FILE = 'mesh_file'

        self.DEFAULT = 'default'

        self.CONFIGURATION = 'configuration' + self.separator()

        self.LOG_LEVEL = self.CONFIGURATION + 'log_level'
        self.CLEAR_OLD_LOGS = self.CONFIGURATION + 'clear_old_logs'
        self.SIMULATIONS_TO_RUN = 'simulations_to_run'
        self.SIMULATIONS = 'simulations'
        self.PHASE = 'phase'
        self.SIMULATION_DIR = 'simulation_dir'

        # Only name changes of this are supported
        self.VERBOSITY = 'verbosity'