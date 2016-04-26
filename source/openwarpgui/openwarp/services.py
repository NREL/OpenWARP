# -*- coding: utf-8 -*-
"""
This Python module provides various service functions.

Updated since version 1.1:
    1. Added support for postprocess and visualization.
    2. Added file path validation for parameters of all related methods.

Updated since version 1.2: Merge Code and Update GUI
    1. Integrate New Nemoh using hdf5 and python.

Updated since version 1.3: OpenWarp - Add Logging Functionality
    1. Added support for logging
"""
__author__ = "caoweiquan322, yedtoss"
__copyright__ = "Copyright (C) 2014-2016 TopCoder Inc. All rights reserved."
__version__ = "1.3"

import collections
import uuid
from settings import *
import os
import time
import subprocess
from multiprocessing import Process, Manager
import multiprocessing as mp
import logging
from openwarp import helper
from nemoh import utility
from nemoh import preprocessor
from nemoh import postprocessor
from nemoh import solver
import warnings
from nemoh.utility import Silence

import fnmatch
import h5py

# This class represents parameters used in the meshing process.
# This class is a subclass of "tuple", and is created using collections.namedtuple factory function.
MeshingParameters = collections.namedtuple('MeshingParameters',
                                           'infile outfile maxh minh fineness grading usetolerance tolerance')

# This class represents parameters used in configuring the application.
# This class is a subclass of "tuple", and is created using collections.namedtuple factory function.
ConfigurationParameters = collections.namedtuple('ConfigurationParameters',
                                                 'logging_level clear_log_flag')

# This class represents parameters used in the simulation process.
# This class is a subclass of "tuple", and is created using collections.namedtuple factory function.
SimulationParameters = collections.namedtuple('SimulationParameters', 'rho g depth xeff yeff wave_frequencies ' +
                                              'min_wave_frequencies max_wave_frequencies wave_directions ' +
                                              'max_wave_direction min_wave_directions floating_bodies ' +
                                              'indiq_solver ires tol_gmres max_iterations save_potential ' +
                                              'green_tabulation_numx green_tabulation_numz ' +
                                              'green_tabulation_simpson_npoints use_ode_influence_coefficients ' +
                                              'use_higher_order num_panel_higher_order b_spline_order ' +
                                              'use_dipoles_implementation thin_panels compute_drift_forces ' +
                                              'compute_yaw_moment remove_irregular_frequencies')

# This class represents a floating body used in the SimulationParameters.
# This class is a subclass of "tuple", and is created using collections.namedtuple factory function.
FloatingBody = collections.namedtuple('FloatingBody', 'mesh_file points panels degrees_of_freedom surge sway ' +
                                      'heave roll_about_cdg pitch_about_cdg yaw_about_cdg ' +
                                      'resulting_generalised_forces force_in_x_direction force_in_y_direction ' +
                                      'force_in_z_direction moment_cdg_force_in_x_direction ' +
                                      'moment_cdg_force_in_y_direction moment_cdg_force_in_z_direction ' +
                                      'additional_info_lines')

# This class represents parameters used in the post-proessing.
# This class is a subclass of "tuple", and is created using collections.namedtuple factory function.
PostprocessingParameters = collections.namedtuple('PostprocessingParameters', 'irf show_pressure ' +
                                                  'kochin_function free_surface_elevation')

# The pre-defined config file name used by MESH_GENERATOR_BIN.
_CONFIG_FILE_NAME = 'config.txt'

# The pre-defined stdout log file name.
_LOG_FILE_NAME = 'log.txt'

class ServiceError(Exception):
    '''
    This exception indicates a service error.
    It will be raised by methods of this module.
    '''
    pass


def _set_log_level(log_level):
        """
        helper method to change level of logs depending of log_level

        :param log_level: integer or string (10 for DEBUG, 20 for INFO)

        :return: A message indicating success or not
        """
        # setting logging level:
        if log_level and str(log_level) in ["10", "20"]:
            # Setting the root logger to that level
            level = int(log_level)
            logging.getLogger(__name__).info("Setting logging level to " +
                                             logging.getLevelName(level))
            logging.getLogger().setLevel(level)
            out = "Logging successfully set to level " + logging.getLevelName(level)
            logging.getLogger(__name__).info(out)
            return out

        else:
            # nofifying user and setting up loggers to debug by default
            logger = logging.getLogger(__name__)
            out = ("Logging level unknown! Should be 10 (DEBUG) or 20 (INFO). Keeping current level of " +
                   logging.getLevelName(logging.getLogger().getEffectiveLevel()))
            logger.warning(out)

            return out


def _clear_log():
    """
    clearing log on server
    :param None
    :return A message indicating success or not
    """
    num = 1
    MAX_NUM = 10000
    while True:
        path = os.path.join(LOG_FILE, str(num))
        if num > MAX_NUM:
            return num-1
        if os.path.exists(path):
            os.remove(path)
            num += 1
        else:
            return num-1

    return num-1


def apply_configuration(params):
    # The logger object for logging.
    logger = logging.getLogger(__name__)
    signature = __name__ + '.apply_configuration(params)'
    helper.log_entrance(logger, signature,
                        {'params': params})



    # Checking parameters
    helper.check_type_value(params, 'params', ConfigurationParameters, False)
    helper.check_not_none_nor_empty(params.logging_level, 'params.logging_level')
    helper.check_not_none_nor_empty(params.clear_log_flag, 'params.clear_log_flag')

    try:
        output = "\n" + _set_log_level(params.logging_level)

        if params.clear_log_flag == "true":
            num = _clear_log()
            out = str(num) + " logs files were successfully cleared"
            logging.getLogger(__name__).info(out)
            output += "\n" + out
        else:
            out = "clearing logs were not enabled"
            output += "\n" + out
            logging.getLogger(__name__).info(out)

        output += ("\n\n" + "The simulation and running time logs are saved in " +
                   LOG_FILE)

        helper.log_exit(logger, signature, [output])

        return output
    except Exception as e:
        # No logging when raising an exception. Another function will log it
        raise ServiceError('Error occurs when generating mesh. Caused by:\n' + unicode(str(e)))





def generate_mesh(meshing_dir, params):
    '''
    Launch Mesh Generator to generate mesh.

    @param meshing_dir: the meshing directory
    @param params: the meshing parameters
    @return: the mesh generation log content
    @raise TypeError: if any input parameter is not of required type
    @raise ValueError: if any input parameter is None/empty, or any field of MeshingParameters is not
                        of valid value
    @raise ServiceError: if error occurs during generating mesh
    '''
    # The logger object for logging.
    logger = logging.getLogger(__name__)
    signature = __name__ + '.generate_mesh()'
    helper.log_entrance(logger, signature,
                        {'meshing_dir': meshing_dir,
                         'params': params})
    # Checking parameters
    helper.check_not_none_nor_empty(meshing_dir, 'meshing_dir')
    helper.check_is_directory(meshing_dir, 'meshing_dir')
    helper.check_type_value(params, 'params', MeshingParameters, False)
    helper.check_not_none_nor_empty(params.infile, 'params.infile')
    helper.check_is_file(params.infile, 'params.infile')
    helper.check_not_none_nor_empty(params.outfile, 'params.outfile')
    helper.check_not_none_nor_empty(params.maxh, 'params.maxh')
    helper.check_not_none_nor_empty(params.minh, 'params.minh')
    helper.check_not_none_nor_empty(params.fineness, 'params.fineness')
    helper.check_not_none_nor_empty(params.grading, 'params.grading')
    helper.check_not_none_nor_empty(params.usetolerance, 'params.usetolerance')
    if params.usetolerance == '1':
        helper.check_not_none_nor_empty(params.tolerance, 'params.tolerance')

    try:
        config_file_path = os.path.join(meshing_dir, _CONFIG_FILE_NAME)
        log_file_path = os.path.join(meshing_dir, _LOG_FILE_NAME)
        # Generate config.txt according to given parameters
        with open(config_file_path, 'w') as f:
            f.write('\n'.join("%s: %s" % item for item in vars(params).items() if item[1] is not None))

        # Launch mesh generator
        with open(log_file_path, 'w') as log_file:
            logger.debug('Start mesh generator in subprocess.')
            subprocess.call(MESH_GENERATOR_BIN, cwd=meshing_dir, stdout=log_file)
            logger.debug('End mesh generator in subprocess.')
        
        # Read and return the log file content
        with open(log_file_path, 'r') as log_file:
            ret = log_file.read().splitlines()
            helper.log_exit(logger, signature, [ret])
            return ret
    except Exception as e:
        helper.log_exception(logger, signature, e)
        raise ServiceError('Error occurs when generating mesh. Caused by:\n' + unicode(str(e)))

def simulate(simulation_dir, params):
    '''
    Run simulation.

    @param simulation_dir: the simulation directory
    @param params: the simulation parameters
    @return: the simulation log content
    @raise TypeError: if any input parameter is not of required type
    @raise ValueError: if any input parameter is None/empty, or any field of SimulationParameters is not
                        of valid value
    @raise ServiceError: if any other error occurred when launching the simulation
    '''
    # The logger object for logging.
    logger = logging.getLogger(__name__)
    signature = __name__ + '.simulate()'
    helper.log_entrance(logger, signature,
                        {'simulation_dir': simulation_dir,
                         'params': params})
    # Checking parameters
    helper.check_not_none_nor_empty(simulation_dir, 'simulation_dir')
    helper.check_is_directory(simulation_dir, 'simulation_dir')
    helper.check_type_value(params, 'params', SimulationParameters, False)
    helper.check_not_none_nor_empty(params.rho, 'params.rho')
    helper.check_not_none_nor_empty(params.g, 'params.g')
    helper.check_not_none_nor_empty(params.depth, 'params.depth')
    helper.check_not_none_nor_empty(params.xeff, 'params.xeff')
    helper.check_not_none_nor_empty(params.yeff, 'params.yeff')
    helper.check_not_none_nor_empty(params.wave_frequencies, 'params.wave_frequencies')
    helper.check_not_none_nor_empty(params.min_wave_frequencies, 'params.min_wave_frequencies')
    helper.check_not_none_nor_empty(params.max_wave_frequencies, 'params.max_wave_frequencies')
    helper.check_not_none_nor_empty(params.wave_directions, 'params.wave_directions')
    helper.check_not_none_nor_empty(params.min_wave_directions, 'params.min_wave_directions')
    helper.check_not_none_nor_empty(params.max_wave_direction, 'params.max_wave_direction')
    helper.check_not_none_nor_empty(params.indiq_solver, 'params.indiq_solver')
    helper.check_not_none_nor_empty(params.ires, 'params.ires')
    helper.check_not_none_nor_empty(params.tol_gmres, 'params.tol_gmres')
    helper.check_not_none_nor_empty(params.max_iterations, 'params.max_iterations')
    helper.check_not_none_nor_empty(params.save_potential, 'params.save_potential')
    helper.check_not_none_nor_empty(params.green_tabulation_numx, 'params.green_tabulation_numx')
    helper.check_not_none_nor_empty(params.green_tabulation_numz, 'params.green_tabulation_numz')
    helper.check_not_none_nor_empty(params.green_tabulation_simpson_npoints, 'params.green_tabulation_simpson_npoints')
    helper.check_not_none_nor_empty(params.use_ode_influence_coefficients, 'params.use_ode_influence_coefficients')
    helper.check_not_none_nor_empty(params.use_higher_order, 'params.use_higher_order')
    helper.check_not_none_nor_empty(params.num_panel_higher_order, 'params.num_panel_higher_order')
    helper.check_not_none_nor_empty(params.b_spline_order, 'params.b_spline_order')
    helper.check_not_none_nor_empty(params.use_dipoles_implementation, 'params.use_dipoles_implementation')
    helper.check_not_none_nor_empty(params.thin_panels, 'params.thin_panels')
    helper.check_not_none_nor_empty(params.compute_drift_forces, 'params.compute_drift_forces')
    helper.check_not_none_nor_empty(params.remove_irregular_frequencies, 'params.remove_irregular_frequencies')
    helper.check_not_none_nor_empty(params.compute_yaw_moment, 'params.compute_yaw_moment')
    
    helper.check_type_value(params.floating_bodies, 'params.floating_bodies', list, True)
    if params.floating_bodies is not None:
        for body in params.floating_bodies:
            helper.check_type_value(body, 'params.floating_bodies item', FloatingBody, False)
            helper.check_not_none_nor_empty(body.mesh_file, 'body.mesh_file')
            helper.check_not_none_nor_empty(body.points, 'body.points')
            helper.check_not_none_nor_empty(body.panels, 'body.panels')
            helper.check_not_none_nor_empty(body.degrees_of_freedom, 'body.degrees_of_freedom')
            helper.check_not_none_nor_empty(body.resulting_generalised_forces, 'body.resulting_generalised_forces')
            helper.check_not_none_nor_empty(body.additional_info_lines, 'body.additional_info_lines')

    try:
        # Write the hdf5 inputs according to given parameters
        # Bug solving in old h5py version: Creating the file first
        hdf5_path = os.path.join(simulation_dir, 'db.hdf5')
        utility.touch(hdf5_path)
        with h5py.File(hdf5_path, "a") as hdf5_data:
            utility.write_calculations(params, hdf5_data)
        
        # Launch preProcessor and Solver
        # A prepared 'results' folder is necessary for the Nemoh software suite
        os.mkdir(os.path.join(simulation_dir, 'results'))
        simulation_log_path = os.path.join(simulation_dir, 'simulation_log.txt')
        custom_config = {
            'HDF5_FILE': hdf5_path,
            'NEMOH_CALCULATIONS_FILE': None,
            'NEMOH_INPUT_FILE': None,
            'MESH_TEC_FILE': os.path.join(simulation_dir, 'mesh', 'mesh.tec'),
            'FK_FORCE_TEC_FILE': os.path.join(simulation_dir, 'results', 'fkforce.tec'),
            'RADIATION_COEFFICIENTS_TEC_FILE': os.path.join(simulation_dir, 'results', 'radiationcoefficients.tec'),
            'DIFFRACTION_FORCE_TEC_FILE': os.path.join(simulation_dir, 'results', 'diffractionforce.tec'),
            'EXCITATION_FORCE_TEC_FILE': os.path.join(simulation_dir, 'results', 'excitationforce.tec'),
            'IRF_TEC_FILE': os.path.join(simulation_dir, 'results', 'irf.tec'),
            'WAVE_FIELD_TEC_FILE': os.path.join(simulation_dir, 'results', 'WaveField.tec'),
            'GREEN_TABULATION_NUMX' : int(params.green_tabulation_numx),
            'GREEN_TABULATION_NUMZ' : int(params.green_tabulation_numz),
            'GREEN_TABULATION_SIMPSON_NPOINTS' : int(params.green_tabulation_simpson_npoints),
            'USE_ODE_INFLUENCE_COEFFICIENTS': bool(int(params.use_ode_influence_coefficients)),
            'USE_HIGHER_ORDER' : bool(int(params.use_higher_order)),
            'NUM_PANEL_HIGHER_ORDER' : int(params.num_panel_higher_order),
            'B_SPLINE_ORDER': int(params.b_spline_order),
            'USE_DIPOLES_IMPLEMENTATION': bool(int(params.use_dipoles_implementation)),
            'THIN_PANELS': [int(i) for i in params.thin_panels.split()],
            'COMPUTE_DRIFT_FORCES' : bool(int(params.compute_drift_forces)),
            'COMPUTE_YAW_MOMENT': bool(int(params.compute_yaw_moment)),
            'REMOVE_IRREGULAR_FREQUENCIES' : bool(int(params.remove_irregular_frequencies))
        }

        logger.debug('Start preProcessor function.')
        run_thread(preprocessor.preprocess, (custom_config,), simulation_log_path)


        logger.debug('End preProcessor function.')
        logger.debug('Start solver function.')
        output = run_thread(solver.solve, (custom_config,), None)



        with open(simulation_log_path, 'a') as log_file:
            log_file.write(output)
        logger.debug('End solver function.')
        with open(simulation_log_path, 'r') as log_file:
            ret = log_file.read().splitlines()
            helper.log_exit(logger, signature, [ret])
            return ret
    except Exception as e:
        helper.log_exception(logger, signature, e)
        raise ServiceError('Error occurs when doing simulation. Caused by:\n' + unicode(str(e)))

def postprocess(simulation_dir, params):
    '''
    Run post-processing.

    @param simulation_dir: the simulation directory
    @param params: the post-processing parameters
    @return: the post-processing log content
    @raise TypeError: if any input parameter is not of required type
    @raise ValueError: if any input parameter is None/empty, or any field of PostprocessingParameters is not
                        of valid value
    @raise ServiceError: if error occurs during launching the post-processing
    '''
    # The logger object for logging.
    logger = logging.getLogger(__name__)
    signature = __name__ + '.postprocess()'
    helper.log_entrance(logger, signature,
                        {'simulation_dir': simulation_dir,
                         'params': params})
    # Checking parameters
    helper.check_not_none_nor_empty(simulation_dir, 'simulation_dir')
    helper.check_is_directory(simulation_dir, 'simulation_dir')
    helper.check_type_value(params, 'params', PostprocessingParameters, False)
    helper.check_type_value(params.irf, 'params.irf', list, False)
    for irf_item in params.irf:
        helper.check_not_none_nor_empty(irf_item, 'irf_item')
    helper.check_not_none_nor_empty(params.show_pressure, 'params.show_pressure')
    helper.check_type_value(params.kochin_function, 'params.kochin_function', list, False)
    for kochin_function_item in params.kochin_function:
        helper.check_not_none_nor_empty(kochin_function_item, 'kochin_function_item')
    helper.check_type_value(params.free_surface_elevation, 'params.free_surface_elevation', list, False)
    for elevation_item in params.free_surface_elevation:
        helper.check_not_none_nor_empty(elevation_item, 'elevation_item')

    try:
        with h5py.File(os.path.join(simulation_dir, 'db.hdf5'), "a") as hdf5_data:
            utility.write_postprocessing_section(params, hdf5_data)

        # Launch postProcessor
        postprocessing_log_path = os.path.join(simulation_dir, 'postprocessing_log.txt')
        custom_config = {
            'HDF5_FILE': os.path.join(simulation_dir, 'db.hdf5'),
            'NEMOH_CALCULATIONS_FILE': None,
            'NEMOH_INPUT_FILE': None,
            'MESH_TEC_FILE': os.path.join(simulation_dir, 'mesh', 'mesh.tec'),
            'FK_FORCE_TEC_FILE': os.path.join(simulation_dir, 'results', 'fkforce.tec'),
            'RADIATION_COEFFICIENTS_TEC_FILE': os.path.join(simulation_dir, 'results', 'radiationcoefficients.tec'),
            'DIFFRACTION_FORCE_TEC_FILE': os.path.join(simulation_dir, 'results', 'diffractionforce.tec'),
            'EXCITATION_FORCE_TEC_FILE': os.path.join(simulation_dir, 'results', 'excitationforce.tec'),
            'IRF_TEC_FILE': os.path.join(simulation_dir, 'results', 'irf.tec'),
            'WAVE_FIELD_TEC_FILE': os.path.join(simulation_dir, 'results', 'WaveField.tec'),
            'GREEN_TABULATION_NUMX' : 328,
            'GREEN_TABULATION_NUMZ' : 46,
            'GREEN_TABULATION_SIMPSON_NPOINTS' : 251,
            'USE_ODE_INFLUENCE_COEFFICIENTS': False,
            'USE_HIGHER_ORDER' : False,
            'NUM_PANEL_HIGHER_ORDER' : 1,
            'B_SPLINE_ORDER': 1,
            'USE_DIPOLES_IMPLEMENTATION': False,
            'THIN_PANELS': [-1],
            'COMPUTE_DRIFT_FORCES' : False,
            'COMPUTE_YAW_MOMENT': False,
            'REMOVE_IRREGULAR_FREQUENCIES' : False
        }
        logger.debug('Start postProcessor function.')
        run_thread(postprocessor.postprocess, (custom_config,), postprocessing_log_path)
        logger.debug('End postProcessor in subprocess.')

        with open(postprocessing_log_path, 'r') as log_file:
            ret = log_file.read().splitlines()
            helper.log_exit(logger, signature, [ret])
            return ret
    except Exception as e:
        helper.log_exception(logger, signature, e)
        raise ServiceError('Error occurs when running postprocess. Caused by:\n' + unicode(str(e)))


def visualize(simulation_dir):
    '''
    Launch ParaView to visualize simulation results.

    @param simulation_dir: the simulation directory
    @raise TypeError: if any input parameter is not of required type
    @raise ValueError: if any input parameter is None/empty
    @raise ServiceError: if error occurs during launching the ParaView
    '''
    # The logger object for logging.
    logger = logging.getLogger(__name__)
    signature = __name__ + '.visualize()'
    helper.log_entrance(logger, signature, {'simulation_dir': simulation_dir})
    # Checking parameters
    helper.check_not_none_nor_empty(simulation_dir, 'simulation_dir')
    helper.check_is_directory(simulation_dir, 'simulation_dir')

    try:
        # Filter files to be opened in ParaView
        files = []
        for f in os.listdir(os.path.join(simulation_dir, 'results')):
            for ext in VISUALIZATION_FILE_EXTENSIONS:
                if fnmatch.fnmatch(f, '*.' + ext):
                    files.append(os.path.join(simulation_dir, 'results', f))

        # Check if there's tec/vtk/stl file to visualize
        if len(files) == 0:
            raise ServiceError('There is no accepted file to visualize.')
        logger.debug('List of files to load:')
        logger.debug(str(files))

        # Prepare script to run by ParaView
        paraview_script = os.path.join(os.path.join(simulation_dir, 'results'), 'load_data.py')
        prepare_paraview_script(paraview_script, files)

        # Launch ParaView without waiting for the ParaView to exit
        logger.debug('Start launching ParaView in subprocess.')
        subprocess.Popen([PARAVIEW_BIN, '--script=' + paraview_script + ''])
        logger.debug('End launching ParaView in subprocess.')
        helper.log_exit(logger, signature, None)
    except Exception as e:
        helper.log_exception(logger, signature, e)
        raise ServiceError('Error occurs when launching the ParaView. Caused by:\n' + unicode(str(e)))

def prepare_paraview_script(script_path, files):
    '''
    Prepare a script to be run by ParaView from a template.

    @param script_path: path of the new script to create
    @param files: a list of data files path
    @raise Exception: to its caller if any error occurs
    '''
    # Since this is a inner function, no entrance/exit information would be logged.
    with open(PARAVIEW_SCRIPT_TEMPLATE, 'r') as fin:
        with open(script_path, 'w') as fout:
            for line in fin.readlines():
                fout.write(line.rstrip().replace('<parameter_files>', str(files)) + '\n')


def wrapper_io(func, fd, args, return_dict):
    """
    Run a function while redirecting its output to a file descriptor
    Args:
        func: A python function to run
        fd: a file descriptor
        args: A tuple containing argument for the function
        return_dict: Dictionary where to put the result of the function

    """
    return_dict["output"] = ''
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        if fd:
            with Silence(stdout=fd, stderr=os.devnull, mode='a'):
                return_dict["output"] = func(*args)
        else:
            return_dict["output"] = func(*args)


def run_thread(func, args, fd):
    """
    Run a python function in a thread and wait for it to complete.
    Redirect its output to fd
    Args:
        func: A python function to run
        args: A tuple containing argument for the function
        fd: a file descriptor
    """

    manager = Manager()
    return_dict = manager.dict()

    p = Process(target=wrapper_io, args=(func, fd, args, return_dict))
    p.start()
    p.join()

    return return_dict["output"]



def writeline_if_not_none(fout, data):
    '''
    Write one line to the specified file if data is not None.

    @param fout: the file object to write line in
    @param data: the data to write as line
    '''
    # Since this is a inner function, no entrance/exit information would be logged.
    if data is not None:
        fout.write(str(data) + '\n')

def prepare_dir(prefix):
    '''
    Prepare a directory, the directory will be a sub-directory of USER_DATA_DIRECTORY with current timestamp
    prefixed given prefix as the directory name.

    @param prefix: the directory prefix
    @return: the meshing/simulation directory full path
    @raise TypeError: if any input parameter is not of required type
    @raise ValueError: if any input parameter is None/empty
    @raise ServiceError: if any error occurred when preparing the directory
    '''
    # The logger object for logging.
    logger = logging.getLogger(__name__)
    signature = __name__ + '.prepare_dir()'
    helper.log_entrance(logger, signature, {'prefix': prefix})
    # Checking parameters
    helper.check_not_none_nor_empty(prefix, 'prefix')
    
    try:
        # Create a directory for this run (sub-directory name in format simulation_YYYYMMDDhhmmss)
        # We should consider adding some more uuid suffix to allow more concurrent requests within 1 SINGLE second.
        run_dir = os.path.join(USER_DATA_DIRECTORY, prefix + time.strftime('%Y%m%d%H%M%S') + '_' + uuid.uuid1().hex)
        os.makedirs(run_dir)
        helper.log_exit(logger, signature, [run_dir])
        return run_dir
    except Exception as e:
        helper.log_exception(logger, signature, e)
        raise ServiceError('Error occurs when preparing the directory. Caused by:\n' + unicode(str(e)))
