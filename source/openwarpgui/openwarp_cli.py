# -*- coding: utf-8 -*-
"""
Copyright (C) 2014-2016 TopCoder Inc., All Rights Reserved.

This module defines the command line version of OpenWarp.


@author:  TCSASSEMBLER
@version: 1.0
"""

import json
import copy
import collections
import sys
from nemoh.models import MeshFormat
from nemoh import utility

from logutils.queue import QueueListener
import multiprocessing
import logging
from openwarp import services
import os
import jmespath as jp
from nemoh.structure import JSON_STRUCTURE
from nemoh import settings
import subprocess
from openwarp import settings as openwarp_settings


# script file
SCRIPT = os.path.realpath(__file__)

# Default configuration file
DEFAULT_CONFIGURATION_FILE = os.path.join(os.path.dirname(SCRIPT), 'configs', 'default_configuration.json')

# This is the structure of the json configuration file
struct = JSON_STRUCTURE()

def merge_config(d, u):
    """
    This function recursively merge one dict with another
    :param d the base dict
    :param u the dict to update d with
    :return the merged dict
    """
    # if u is a dict and empty we don't update d which is the truth
    for k, v in u.iteritems():
        if isinstance(v, collections.Mapping):
            s = d.get(k, {})
            # this means that if there is mismatch we don't do anything which is correct here
            if isinstance(s, collections.Mapping):
                r = merge_config(d.get(k, {}), v)
                d[k] = r
        else:
            d[k] = u[k]
    return d

def format_list(lst):
    """
    Format a list as a string, ignore if it is string
    :param lst the list to format
    :return the formatted string
    """
    if not isinstance(lst, basestring):
        return " ".join(str(i) for i in lst)
    return lst

def convert_dict_values(d):
    """
    This function makes sure the top level value a dict are str
    :param d the dict to convert
    :return the converted dict
    """
    for key, value in d.iteritems():

        if isinstance(value, bool):
            value = int(value)
        if isinstance(value, (list, tuple)):
            d[key] = [str(i) for i in value]
        elif not isinstance(value, collections.Mapping):
            d[key] = str(value)

    return d


def update_forces(idx, forces):
    """
    This function return the full array for the generalized forces and moment
    :param the index of the force to update
    :param the partial forces value
    :return the full forces
    """

    f = [0 for i in range(7)]
    f[0] = idx/3+1
    f[idx%3+1] = 1
    tmp = forces[-3:]
    f[-len(tmp):] = tmp
    return format_list(f)




def get_floating_body(fb, simulation_dir):
    """
    This function generate a floating body and make convert the mesh file to dat format if not
    already in dat format
    :param simulation_dir the simulation directory (str)
    :param fb : the raw floating body
    :return the updated floating body
    """

    mesh_file = os.path.realpath(jp.search(struct.BODY_MESH_FILE, fb))

    
    if not MeshFormat.is_dat_file(mesh_file):
        meshing_param = {
            "maxh": jp.search(struct.BODY_MESHING_MAXH, fb),
            "minh": jp.search(struct.BODY_MESHING_MINH, fb),
            "fineness": jp.search(struct.BODY_MESHING_FINENESS, fb),
            "grading": jp.search(struct.BODY_MESHING_GRADING, fb),
            "usetolerance": jp.search(struct.BODY_MESHING_USE_TOLERANCE, fb),
            "tolerance": jp.search(struct.BODY_MESHING_TOLERANCE, fb),
            "infile": mesh_file,
            "outfile": os.path.splitext(os.path.basename(mesh_file))[0]
        }
        meshing_dir = os.path.join(simulation_dir, 'meshing')
        utility.mkdir_p(meshing_dir)

        meshing_param = convert_dict_values(meshing_param)
        services.generate_mesh(meshing_dir, services.MeshingParameters(**meshing_param))
        mesh_file = os.path.join(meshing_dir, meshing_param['outfile'] + '-quad.dat')
        

    with open(mesh_file, 'r') as f:
        n_points, n_panels = utility.determine_points_panels(f)

    body_param = {
        'degrees_of_freedom': 6,
        'additional_info_lines': 0,
        'resulting_generalised_forces': 6,
        'surge': update_forces(0, jp.search(struct.BODY_SURGE, fb)),
        'sway': update_forces(1, jp.search(struct.BODY_SWAY, fb)),
        'heave': update_forces(2, jp.search(struct.BODY_HEAVE, fb)),
        'roll_about_cdg': update_forces(3, jp.search(struct.BODY_ROLL, fb)),
        'pitch_about_cdg': update_forces(4, jp.search(struct.BODY_PITCH, fb)),
        'yaw_about_cdg': update_forces(5, jp.search(struct.BODY_YAW, fb)),
        'force_in_x_direction': update_forces(0, jp.search(struct.BODY_FORCE_X, fb)),
        'force_in_y_direction': update_forces(1, jp.search(struct.BODY_FORCE_Y, fb)),
        'force_in_z_direction': update_forces(2, jp.search(struct.BODY_FORCE_Z, fb)),
        'moment_cdg_force_in_x_direction': update_forces(3, jp.search(struct.BODY_MOMENT_X, fb)),
        'moment_cdg_force_in_y_direction': update_forces(4, jp.search(struct.BODY_MOMENT_Y, fb)),
        'moment_cdg_force_in_z_direction': update_forces(5, jp.search(struct.BODY_MOMENT_Z, fb)),
        'mesh_file': mesh_file,
        'points': n_points,
        'panels': n_panels
    }

    body_param = convert_dict_values(body_param)
    return body_param


def generate_floating_bodies(simulation_parameter, simulation_dir):
    """
    This function runs the meshing phase
    :param simulation_parameter  the simulation parameter as dict
    :param simulation_dir the simulation directory (str)
    """
    floating_bodies = jp.search(struct.H5_BODIES, simulation_parameter)
    fbs = []
    for name, body in floating_bodies.iteritems():
        if name == struct.DEFAULT: continue
        fbs.append(get_floating_body(body, simulation_dir))

    return fbs


def run_simulation(simulation_parameter, simulation_dir, queue):
    """
    This function runs the simulation phase
    :param simulation_parameter  the simulation parameter as dict
    :param simulation_dir the simulation directory (str)
    :param queue : the logging listener queue
    """
    wave_point = jp.search(struct.H5_ENV_WAVE_POINT, simulation_parameter)

    simulation_param = {
        'rho': jp.search(struct.H5_ENV_VOLUME, simulation_parameter),
        'g': jp.search(struct.H5_ENV_GRAVITY, simulation_parameter),
        'depth': jp.search(struct.H5_ENV_DEPTH, simulation_parameter),
        'xeff': wave_point[0],
        'yeff': wave_point[1],
        'wave_frequencies': jp.search(struct.H5_NUM_WAVE_FREQUENCIES, simulation_parameter),
        'min_wave_frequencies': jp.search(struct.H5_MIN_WAVE_FREQUENCIES, simulation_parameter),
        'max_wave_frequencies': jp.search(struct.H5_MAX_WAVE_FREQUENCIES, simulation_parameter),
        'wave_directions' : jp.search(struct.H5_NUM_WAVE_DIRECTIONS, simulation_parameter),
        'max_wave_direction': jp.search(struct.H5_MAX_WAVE_DIRECTIONS, simulation_parameter),
        'min_wave_directions': jp.search(struct.H5_MIN_WAVE_DIRECTIONS, simulation_parameter),
        'indiq_solver': jp.search(struct.H5_SOLVER_TYPE, simulation_parameter),
        'ires': jp.search(struct.H5_SOLVER_GMRES_RESTART, simulation_parameter),
        'tol_gmres': jp.search(struct.H5_SOLVER_GMRES_STOPPING, simulation_parameter),
        'max_iterations': jp.search(struct.H5_SOLVER_GMRES_MAX_ITERATIONS, simulation_parameter),
        'save_potential': 1,
        'green_tabulation_numx': jp.search(struct.H5_SOLVER_GREEN_TABULATION_NUMX, simulation_parameter),
        'green_tabulation_numz' : jp.search(struct.H5_SOLVER_GREEN_TABULATION_NUMZ, simulation_parameter),
        'green_tabulation_simpson_npoints': jp.search(struct.H5_SOLVER_GREEN_TABULATION_SIMPSON_NPOINTS, simulation_parameter),
        'use_ode_influence_coefficients': jp.search(struct.H5_SOLVER_SWITCH_ODE_INFLUENCE, simulation_parameter),
        'use_higher_order': jp.search(struct.H5_SOLVER_USE_HIGHER_ORDER, simulation_parameter),
        'num_panel_higher_order': jp.search(struct.H5_SOLVER_NUM_PANEL_HIGHER_ORDER, simulation_parameter),
        'b_spline_order': jp.search(struct.H5_SOLVER_B_SPLINE_ORDER, simulation_parameter),
        'use_dipoles_implementation': jp.search(struct.H5_SOLVER_USE_DIPOLES_IMPLEMENTATION, simulation_parameter),
        'thin_panels': format_list(jp.search(struct.H5_SOLVER_THIN_PANELS , simulation_parameter)),
        'compute_drift_forces' : jp.search(struct.H5_SOLVER_COMPUTE_DRIFT_FORCES, simulation_parameter),
        'compute_yaw_moment': jp.search(struct.H5_SOLVER_COMPUTE_YAW_MOMENT, simulation_parameter),
        'remove_irregular_frequencies' : jp.search(struct.H5_SOLVER_REMOVE_IRREGULAR_FREQUENCIES, simulation_parameter)
    }

    simulation_param = convert_dict_values(simulation_param)
    simulation_param['floating_bodies'] = generate_floating_bodies(simulation_parameter, simulation_dir)

    return services.simulate(simulation_dir, services.construct_simulation_parameters(simulation_param), queue)



def post_process(simulation_parameter, simulation_dir, queue):
    """
    This function runs the post processing phase
    :param simulation_parameter  the simulation parameter as dict
    :param simulation_dir the simulation directory (str)
    :param queue : the logging listener queue
    """
    post_process_param = {
        'irf': [int(jp.search(struct.H5_COMPUTE_IRF, simulation_parameter)), 
                jp.search(struct.H5_IRF_TIME_STEP, simulation_parameter),
                jp.search(struct.H5_IRF_DURATION, simulation_parameter)],
        'show_pressure': jp.search(struct.H5_SHOW_PRESSURE, simulation_parameter),
        'kochin_function': [jp.search(struct.H5_KOCHIN_NUMBER, simulation_parameter),
                            jp.search(struct.H5_KOCHIN_MIN, simulation_parameter),
                            jp.search(struct.H5_KOCHIN_MAX, simulation_parameter)],
        'free_surface_elevation': [jp.search(struct.H5_FREE_SURFACE_POINTS_X, simulation_parameter),
                                    jp.search(struct.H5_FREE_SURFACE_POINTS_Y, simulation_parameter),
                                    jp.search(struct.H5_FREE_SURFACE_DIMENSION_X, simulation_parameter),
                                    jp.search(struct.H5_FREE_SURFACE_DIMENSION_Y, simulation_parameter)]
    }

    post_process_param = convert_dict_values(post_process_param)

    return services.postprocess(simulation_dir, services.construct_postprocess_parameters(post_process_param), queue)



def apply_configuration(config):
    """
    This function applies the simulations configurations
    :param config the configuration file
    """
    log_level = str(jp.search(struct.LOG_LEVEL, config)).upper()
    log_level = 10 if log_level == 'DEBUG' else 20
    configuration_parameter = {

        'logging_level': log_level,
        'clear_log_flag': jp.search(struct.CLEAR_OLD_LOGS, config)
    }
    configuration_parameter = convert_dict_values(configuration_parameter)

    services.apply_configuration(services.ConfigurationParameters(**configuration_parameter))



def run_phases(simulation_parameter, simulation_dir, phase, queue):
    """
    This function run the simulation, meshing and or postprocessing phases
    :param simulation_parameter  the simulation parameter as dict
    :param simulation_dir the simulation directory (str)
    :param phase the phase to run (str). If None all phases are run
    :param queue : the logging listener queue
    """


    if phase == 'MESHING':
        generate_floating_bodies(simulation_parameter, simulation_dir)
    elif phase == 'POSTPROCESSING':
        post_process(simulation_parameter, simulation_dir, queue)
    elif phase == 'SIMULATION':
        run_simulation(simulation_parameter, simulation_dir, queue)
    elif not phase:
        run_simulation(simulation_parameter, simulation_dir, queue)
        post_process(simulation_parameter, simulation_dir, queue)




def run(user_config, queue):
    """
    This method takes a user configuration parse it and fill in all default
    :param user_config the user configuration dictionary
    :param queue : the logging listener queue
    :return the filled in configuration dictionary
    """

    # Get application based default configuration
    with open(DEFAULT_CONFIGURATION_FILE, 'rt') as f:
        default_config = json.load(f)
    # Merge default configuration with user configuration
    config = merge_config(copy.deepcopy(default_config), copy.deepcopy(user_config))
    
    # Get all simulations
    all_simulations = jp.search(struct.SIMULATIONS, config)

    # Get the simulations to run
    simulations = all_simulations.keys()
    if jp.search(struct.SIMULATIONS_TO_RUN, config): simulations = jp.search(struct.SIMULATIONS_TO_RUN, config)
    if struct.DEFAULT in simulations: simulations.remove(struct.DEFAULT)

    
    # Get the default simulation parameter
    default_simulation_parameters = all_simulations[struct.DEFAULT]

    # Get the phase of the simulations

    phase = jp.search(struct.PHASE, config)

    if phase == 'APPLY_CONFIGURATION' or not phase:
        apply_configuration(config)

    logger = logging.getLogger(__name__)

    logger.info('merged configuration ' + json.dumps(config) + '\n')

    verbosity = jp.search(struct.VERBOSITY, config)

    if verbosity == 0:
        return False

    
    # Update configuration of each floating body for each simulation
    for simulation in simulations:

        # Merge default simulation parameter with this simulation
        simulation_parameter = merge_config(copy.deepcopy(default_simulation_parameters),
        copy.deepcopy(all_simulations[simulation]))

        # Get the floating bodies
        floating_bodies = jp.search(struct.H5_BODIES, simulation_parameter)
        # Merge the default floating boadies parameter with this body parameter
        if struct.DEFAULT in floating_bodies:
            default_floating_bodies = floating_bodies[struct.DEFAULT]
            for name, body in floating_bodies.iteritems():
                if name == struct.DEFAULT: continue
                floating_bodies[name] = merge_config(copy.deepcopy(default_floating_bodies), copy.deepcopy(body))


        logger.info('Running simulation ' + str(simulation) + ' with parameter ' + json.dumps(simulation_parameter) + '\n')

        # Now we have the correct parameter to run in simulation_parameter
        simulation_dir = os.path.realpath(jp.search(struct.SIMULATION_DIR, simulation_parameter))
        # Create the simulation directory if not exists
        utility.mkdir_p(simulation_dir)
        run_phases(simulation_parameter, simulation_dir, phase, queue)

    return True



if __name__ == '__main__':
    """
    Main module, only called when the script is started from the command line
    """
    log_file_found = utility.setup_logging(default_conf_path=settings.LOGGING_CONFIGURATION_FILE,
     logging_path=openwarp_settings.LOG_FILE)
    # Compile python module if it was not compiled.
    # This should always output to the terminal no matter the verbosity level
    subprocess.call(['python', 'setup.py', 'build_ext', '--inplace'], cwd='nemoh')

    logger = logging.getLogger(__name__)

    if len(sys.argv) <= 1:
        utility.log_and_print(logger, 'Error: No configurations file given to the CLI. Usage of script is: openwarp_cli configuration1 .. ')

    queue = multiprocessing.Queue(-1)
    ql = QueueListener(queue, *logging.getLogger().handlers)
    ql.start()

    for i in range(1, len(sys.argv)):
        path = sys.argv[i]
        user_config = None
        # Checking if it is a valid path
        if os.path.exists(path):
            utility.log_and_print(logger, 'Processing configuration file at ' + path + '\n')
            with open(path, 'rt') as f:
                user_config = json.load(f)
        else: # Check if it is json string
            try:
                user_config = json.loads(path)
            except Exception as e:
                user_config = None
        
        if user_config:
            logger.info('Found configuration ' + json.dumps(user_config) + '\n')
            # If there is no log file found, then we must show every thing on terminal
            if not log_file_found : user_config[struct.VERBOSITY] = 2

            # Disable output on terminal if verbosity is disabled
            if not run(user_config, queue):
                DEVNULL = open(os.devnull, 'w')
                user_config[struct.VERBOSITY] = 2
                subprocess.call((sys.executable, SCRIPT, json.dumps(user_config)), stdout=DEVNULL, stderr=DEVNULL)
        else:
            utility.log_and_print(logger, 'Error: Configuration file/json ' + str(path) + ' not found/not valid')


    ql.stop()