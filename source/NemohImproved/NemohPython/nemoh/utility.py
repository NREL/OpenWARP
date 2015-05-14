#!/usr/bin/env python
"""
This module defines utility methods used by the application

Changes in version 1.1:
    Added possibility to run the code with custom settings

Changes in version 1.2 (Irregular Frequencies Assembly)
        Added require_dataset function to allow dataset to be resized
        automatically.
        Applied some bug fixes to allow the shape of hdf5 file dataset 
        to be automatically resized.
"""
import numpy as np
import sys
import structure
import settings
import os
import errno
from models import TEnvironment


__author__ = "yedtoss, TCSASSEMBLER"
__copyright__ = "Copyright (C) 2014-2015 TopCoder Inc. All rights reserved."
__version__ = "1.2"


EPS = 1e-7
"""
The epsilon value 1e-7
"""

EPS_6 = 1e-6
"""
The epsilon value 1e-6
"""

II = complex(0, 1)
"""
Complex number identity
"""


def cih(k, z, h):
    """
    Computes COSH(k*(z+h))/COSH(kh)
    Args:
        k: float, number
        z: float, number
        h: float, number
    """

    if 0 < k*h <= 20:
        result = np.cosh(k*(z+h))/np.cosh(k*h)
    else:
        result = np.exp(k*z)

    return result


def sih(k, z, h):
    """
    Computes SINH(k*(z+h))/COSH(kh)
    Args:
        k: float, number
        z: float, number
        h: float, number
    """

    if 0 < k*h <= 20:
        result = np.sinh(k*(z+h))/np.cosh(k*h)
    else:
        result = np.exp(k*z)

    return result


def require_dataset(hdf5_data, path, shape, dtype, maxshape=(None)):
    """
    Create or update a dataset, making sure that its shape is resized
    if needed
    Args:
        hdf5_data: object, an already opened hdf5 file
        path: string, the path to the dataset
        shape: tuple of integers, the shape of the dataset
        dtype: string or int, the type of the dataset
        maxshape: tuple of integers, the maximum shape to which the dataset can be 
                  resized to. (Unused currently)
    Returns:
        The dataset newly created or updated.
    """
    dset = hdf5_data.get(path, default = None)
    # Dataset not existing
    if dset is None:
        maxshape = [None for i in xrange(len(shape))]
        dset = hdf5_data.create_dataset(path, shape, dtype, maxshape=tuple(maxshape))
    else:
        # Dataset is already existing
        dset.resize(shape)

    return dset




def get_setting(default_setting, custom_config, config_key):
    """
    Args:
        default_setting:  The default setting to use
        custom_config: dict, the custom configuration dictionary
        config_key: string, the key of the setting in the custom_config variable

    Returns:
        default_settings if config_key is not in custom_config,
        otherwise return custom_config[config_key]
    """

    if custom_config and config_key in custom_config:
        return custom_config[config_key]
    else:
        return default_setting


def mkdir_p(path):
    """
    Create a directory if not exists
    Args:
        path: The path of the directory to create
    Raises:
        Exception if the directory can not be created when it does not exist
    """
    try:
        os.makedirs(path)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise


def validate_string(s, name=''):
    """
    Validate a required string
    Args:
        s: string, the string to validate
        name: it's name
    """
    assert (s is not None), name + ' settings should not be None'
    assert (isinstance(s, str)), name + ' settings should be a string'
    assert (s != ''), name + ' settings should be not be empty'


def validate_file(inp, name=''):
    """
    Validate a required string
    Args:
        s: string, the string to validate
        name: it's name
    """
    validate_string(inp, name)
    assert (os.path.exists(inp)), name + ' settings with value ' + inp + ' should exist.'


def read_environment(hdf5_data):
    """
    Read the environment from the hdf5 file
    Args:
        hdf5_data: object, the hdf5 opened storage

    Returns:
        The environment
    """
    environment = TEnvironment()
    environment.rho = hdf5_data.get(structure.H5_ENV_VOLUME)[0]
    environment.g = hdf5_data.get(structure.H5_ENV_GRAVITY)[0]
    environment.depth = hdf5_data.get(structure.H5_ENV_DEPTH)[0]
    wave_point = hdf5_data.get(structure.H5_ENV_WAVE_POINT)
    environment.x_eff = wave_point[0]
    environment.y_eff = wave_point[1]
    return environment


def compute_wave_number(w, environment):
    """
    Computes the wave number
    Args:
        w: float, the wave frequency
        environment: object, the environment
    Returns:
        The wave number
    """

    wave_number = w*w/environment.g
    x0= wave_number*environment.depth
    n_item_x=10000

    if 0 < x0 <= 20:
        xg=0.
        xd=x0
        n_ite=0

        while n_ite < n_item_x and (x0 - xg*np.tanh(xg))*(x0 - xd*np.tanh(xd)) > 0:
            xg = xd
            xd = 2*xd
            n_ite += 1


        n_ite = 0
        # Weird this will never happen. is the n_ite = 0 statement correct?
        if n_ite >= n_item_x:
            print('Error: unable to find the wavenumber')
            sys.exit()

        xc = 0.5*(xd+xg)

        while n_ite < n_item_x and np.abs(xd - xg)/ np.abs(xc) >= 1e-6:
            xc=0.5*(xd+xg)
            if (x0 - xg*np.tanh(xg)) * (x0 - xc*np.tanh(xc)) > 0:
                xg = xc
            else:
                xd = xc
            n_ite += 1

        if n_ite >= n_item_x:
            print('Error: unable to find the wavenumber')
            sys.exit()

        wave_number = xc/environment.depth

    return wave_number


def set_hdf5_attributes(dset, attributes):
    """
    Set the attributes for a hdf5 group or dataset
    Args:
        dset: instance of hdf5.File, hdf5.Group of hdf5.DataSet, the hdf5 group or dataset to set
        attributes: dict, dictionary of attributes to set
    Returns:
        The group or attributes with attributes set
    """
    for key in attributes.iterkeys():
        dset.attrs[key] = attributes[key]

    return dset


def convert_input(filename, hdf5_data):
    """
    Convert nemoh old input file to hdf5_data
    Args:
        filename: string, path to the old nemoh input file
        hdf5_data: object, an already opened hdf5 file
    """
    x1 = []
    with open(filename, 'r') as inp:
        for line in inp:
            x1.append(line)
    idx = 1
    dset = require_dataset(hdf5_data, structure.H5_SOLVER_TYPE, (1,), dtype=settings.NEMOH_INT)
    dset[0] = int(float(x1[idx].split()[0]))
    set_hdf5_attributes(dset, structure.H5_SOLVER_TYPE_ATTR)
    idx += 1

    dset = require_dataset(hdf5_data, structure.H5_SOLVER_GMRES_RESTART, (1,), dtype=settings.NEMOH_INT)
    dset[0] = int(float(x1[idx].split()[0]))
    set_hdf5_attributes(dset, structure.H5_SOLVER_GMRES_RESTART_ATTR)
    idx += 1

    dset = require_dataset(hdf5_data, structure.H5_SOLVER_GMRES_STOPPING, (1,), dtype=settings.NEMOH_FLOAT)
    dset[0] = float(x1[idx].split()[0])
    set_hdf5_attributes(dset, structure.H5_SOLVER_GMRES_STOPPING_ATTR)
    idx += 1

    dset = require_dataset(hdf5_data, structure.H5_SOLVER_GMRES_MAX_ITERATIONS, (1,), dtype=settings.NEMOH_INT)
    dset[0] = int(float(x1[idx].split()[0]))
    set_hdf5_attributes(dset, structure.H5_SOLVER_GMRES_MAX_ITERATIONS_ATTR)


def write_calculations(params, hdf5_data):
    """
    Write calculations parameter from params to hdf5_data
    Args:
        params: SimulationParameters object, params to write
        hdf5_data: object, an already opened hdf5 file
    """

    if params.rho is not None:
        dset = require_dataset(hdf5_data, structure.H5_ENV_VOLUME, (1,), dtype=settings.NEMOH_FLOAT)
        dset[0] = float(params.rho)
        set_hdf5_attributes(dset, structure.H5_ENV_VOLUME_ATTR)

    if params.g is not None:
        dset = require_dataset(hdf5_data, structure.H5_ENV_GRAVITY, (1,), dtype=settings.NEMOH_FLOAT)
        dset[0] = float(params.g)
        set_hdf5_attributes(dset, structure.H5_ENV_GRAVITY_ATTR)

    if params.depth is not None:
        dset = require_dataset(hdf5_data, structure.H5_ENV_DEPTH, (1,), dtype=settings.NEMOH_FLOAT)
        dset[0] = float(params.depth)
        set_hdf5_attributes(dset, structure.H5_ENV_DEPTH_ATTR)

    if (params.xeff is not None) and (params.yeff is not None):
        dset = require_dataset(hdf5_data, structure.H5_ENV_WAVE_POINT, (2,), dtype=settings.NEMOH_FLOAT)
        dset[0] = float(params.xeff)
        dset[1] = float(params.yeff)
        set_hdf5_attributes(dset, structure.H5_ENV_WAVE_POINT_ATTR)

    if params.floating_bodies is not None:
        num_bodies = len(params.floating_bodies)
        i = 0
        for fb in params.floating_bodies:
            i += 1
            body = structure.H5_BODIES + structure.H5_BODY_BASE + str(i) + '/'
            mesh_x = []
            with open(fb.mesh_file, 'r') as mesh_file:
                for line in mesh_file:
                    mesh_x.append(line)

            num_points = int(float(fb.points))
            num_panels = int(float(fb.panels))
            dset = require_dataset(hdf5_data, body + structure.H5_BODY_NUM_POINTS, (1, ), dtype=settings.NEMOH_INT)
            dset[0] = num_points
            set_hdf5_attributes(dset, structure.H5_BODY_NUM_POINTS_ATTR)

            dset = require_dataset(hdf5_data, body + structure.H5_BODY_NUM_PANELS, (1, ), dtype=settings.NEMOH_INT)
            dset[0] = num_panels
            set_hdf5_attributes(dset, structure.H5_BODY_NUM_PANELS_ATTR)

            mesh_idx = 0
            dset = require_dataset(hdf5_data, body + structure.H5_BODY_MESH, (num_points+num_panels+1, 4),
                                             dtype=settings.NEMOH_FLOAT)
            mesh_x2 = mesh_x[mesh_idx].split()
            set_hdf5_attributes(dset, structure.H5_BODY_MESH_ATTR)

            dset[0, 0] = int(float(mesh_x2[0]))
            dset[0, 1] = int(float(mesh_x2[1]))

            for j in range(1, num_points+num_panels+1):
                mesh_idx += 1
                mesh_x2 = mesh_x[mesh_idx].split()
                dset[j, :] = [float(x) for x in mesh_x2[:4]]

                if j == num_points:
                    mesh_idx += 1

            num = int(float(fb.degrees_of_freedom))
            dset = require_dataset(hdf5_data, body + structure.H5_FREEDOM_DEGREE, (num, 7), dtype=settings.NEMOH_FLOAT)
            set_hdf5_attributes(dset, structure.H5_FREEDOM_DEGREE_ATTR)

            x1 = [fb.surge, fb.sway, fb.heave, fb.roll_about_cdg, fb.pitch_about_cdg, fb.yaw_about_cdg]
            for j in range(len(x1)):
                if x1[j]:
                    x2 = x1[j].split()
                    dset[j, :] = np.array([float(x) for x in x2[:7]])

            num = int(float(fb.resulting_generalised_forces))
            dset = require_dataset(hdf5_data, body + structure.H5_GENERALISED_FORCES, (num, 7),
                                             dtype=settings.NEMOH_FLOAT)
            set_hdf5_attributes(dset, structure.H5_GENERALISED_FORCES_ATTR)
            x1 = [fb.force_in_x_direction, fb.force_in_y_direction, fb.force_in_z_direction,
                  fb.moment_cdg_force_in_x_direction, fb.moment_cdg_force_in_y_direction,
                  fb.moment_cdg_force_in_z_direction]
            for j in range(len(x1)):
                if x1[j]:
                    x2 = x1[j].split()
                    dset[j, :] = [float(x) for x in x2[:7]]

    if params.wave_frequencies is not None:
        dset = require_dataset(hdf5_data, structure.H5_NUM_WAVE_FREQUENCIES, (1,), dtype=settings.NEMOH_INT)
        set_hdf5_attributes(dset, structure.H5_NUM_WAVE_FREQUENCIES_ATTR)
        dset[0] = int(float(params.wave_frequencies))

    if params.min_wave_frequencies is not None:
        dset = require_dataset(hdf5_data, structure.H5_MIN_WAVE_FREQUENCIES, (1,), dtype=settings.NEMOH_FLOAT)
        set_hdf5_attributes(dset, structure.H5_MIN_WAVE_FREQUENCIES_ATTR)
        dset[0] = float(params.min_wave_frequencies)

    if params.max_wave_frequencies is not None:
        dset = require_dataset(hdf5_data, structure.H5_MAX_WAVE_FREQUENCIES, (1,), dtype=settings.NEMOH_FLOAT)
        set_hdf5_attributes(dset, structure.H5_MAX_WAVE_FREQUENCIES_ATTR)
        dset[0] = float(params.max_wave_frequencies)

    if params.wave_directions is not None:
        dset = require_dataset(hdf5_data, structure.H5_NUM_WAVE_DIRECTIONS, (1,), dtype=settings.NEMOH_INT)
        set_hdf5_attributes(dset, structure.H5_NUM_WAVE_DIRECTIONS_ATTR)
        dset[0] = int(params.wave_directions)

    if params.min_wave_directions is not None:
        dset = require_dataset(hdf5_data, structure.H5_MIN_WAVE_DIRECTIONS, (1,), dtype=settings.NEMOH_FLOAT)
        set_hdf5_attributes(dset, structure.H5_MIN_WAVE_DIRECTIONS_ATTR)
        dset[0] = float(params.min_wave_directions)

    if params.max_wave_direction is not None:
        dset = require_dataset(hdf5_data, structure.H5_MAX_WAVE_DIRECTIONS, (1,), dtype=settings.NEMOH_FLOAT)
        set_hdf5_attributes(dset, structure.H5_MAX_WAVE_DIRECTIONS_ATTR)
        dset[0] = float(params.max_wave_direction)

    x1 = ['1 0.1 10.', '0', '181. 0. 180.', '1 2 1000. 2.']
    idx = 0
    x2 = x1[idx].split()

    dset = require_dataset(hdf5_data, structure.H5_COMPUTE_IRF, (1,), dtype=settings.NEMOH_INT)
    set_hdf5_attributes(dset, structure.H5_COMPUTE_IRF_ATTR)
    dset[0] = int(x2[0])

    dset = require_dataset(hdf5_data, structure.H5_IRF_TIME_STEP, (1,), dtype=settings.NEMOH_FLOAT)
    set_hdf5_attributes(dset, structure.H5_IRF_TIME_STEP_ATTR)
    dset[0] = float(x2[1])
    dset = require_dataset(hdf5_data, structure.H5_IRF_DURATION, (1,), dtype=settings.NEMOH_FLOAT)
    set_hdf5_attributes(dset, structure.H5_IRF_DURATION_ATTR)
    dset[0] = float(x2[2])

    idx += 1
    x2 = x1[idx].split()
    dset = require_dataset(hdf5_data, structure.H5_SHOW_PRESSURE, (1,), dtype=settings.NEMOH_INT)
    set_hdf5_attributes(dset, structure.H5_SHOW_PRESSURE_ATTR)
    dset[0] = int(x2[0])

    idx += 1
    x2 = x1[idx].split()
    dset = require_dataset(hdf5_data, structure.H5_KOCHIN_NUMBER, (1,), dtype=settings.NEMOH_FLOAT)
    set_hdf5_attributes(dset, structure.H5_KOCHIN_NUMBER_ATTR)
    dset[0] = float(x2[0])
    dset = require_dataset(hdf5_data, structure.H5_KOCHIN_MIN, (1,), dtype=settings.NEMOH_FLOAT)
    set_hdf5_attributes(dset, structure.H5_KOCHIN_MIN_ATTR)
    dset[0] = float(x2[1])
    dset = require_dataset(hdf5_data, structure.H5_KOCHIN_MAX, (1,), dtype=settings.NEMOH_FLOAT)
    set_hdf5_attributes(dset, structure.H5_KOCHIN_MAX_ATTR)
    dset[0] = float(x2[2])

    idx += 1
    x2 = x1[idx].split()
    dset = require_dataset(hdf5_data, structure.H5_FREE_SURFACE_POINTS_X, (1,), dtype=settings.NEMOH_INT)
    set_hdf5_attributes(dset, structure.H5_FREE_SURFACE_POINTS_X_ATTR)
    dset[0] = int(x2[0])
    dset = require_dataset(hdf5_data, structure.H5_FREE_SURFACE_POINTS_Y, (1,), dtype=settings.NEMOH_INT)
    set_hdf5_attributes(dset, structure.H5_FREE_SURFACE_POINTS_Y_ATTR)
    dset[0] = int(x2[1])
    dset = require_dataset(hdf5_data, structure.H5_FREE_SURFACE_DIMENSION_X, (1,), dtype=settings.NEMOH_FLOAT)
    set_hdf5_attributes(dset, structure.H5_FREE_SURFACE_DIMENSION_X_ATTR)
    dset[0] = float(x2[2])
    dset = require_dataset(hdf5_data, structure.H5_FREE_SURFACE_DIMENSION_Y, (1,), dtype=settings.NEMOH_FLOAT)
    set_hdf5_attributes(dset, structure.H5_FREE_SURFACE_DIMENSION_Y_ATTR)
    dset[0] = float(x2[3])

    if params.indiq_solver is not None:
        dset = require_dataset(hdf5_data, structure.H5_SOLVER_TYPE, (1,), dtype=settings.NEMOH_INT)
        dset[0] = int(float(params.indiq_solver))
        set_hdf5_attributes(dset, structure.H5_SOLVER_TYPE_ATTR)

    if params.ires is not None:
        dset = require_dataset(hdf5_data, structure.H5_SOLVER_GMRES_RESTART, (1,), dtype=settings.NEMOH_INT)
        dset[0] = int(float(params.ires))
        set_hdf5_attributes(dset, structure.H5_SOLVER_GMRES_RESTART_ATTR)

    if params.tol_gmres is not None:
        dset = require_dataset(hdf5_data, structure.H5_SOLVER_GMRES_STOPPING, (1,), dtype=settings.NEMOH_FLOAT)
        dset[0] = float(params.tol_gmres)
        set_hdf5_attributes(dset, structure.H5_SOLVER_GMRES_STOPPING_ATTR)

    if params.max_iterations is not None:
        dset = require_dataset(hdf5_data, structure.H5_SOLVER_GMRES_MAX_ITERATIONS, (1,), dtype=settings.NEMOH_INT)

        dset[0] = int(float(params.max_iterations))
        set_hdf5_attributes(dset, structure.H5_SOLVER_GMRES_MAX_ITERATIONS_ATTR)


def write_postprocessing_section(params, hdf5_data):
    """
    Write calculations parameter from params to hdf5_data
    Args:
        params: SimulationParameters object, params to write
        hdf5_data: object, an already opened hdf5 file
    """

    if params.irf is not None:
        x2 = (' '.join(params.irf)).split()
        dset = require_dataset(hdf5_data, structure.H5_COMPUTE_IRF, (1,), dtype=settings.NEMOH_INT)
        set_hdf5_attributes(dset, structure.H5_COMPUTE_IRF_ATTR)
        dset[0] = int(float(x2[0]))

        dset = require_dataset(hdf5_data, structure.H5_IRF_TIME_STEP, (1,), dtype=settings.NEMOH_FLOAT)
        set_hdf5_attributes(dset, structure.H5_IRF_TIME_STEP_ATTR)
        dset[0] = float(x2[1])
        dset = require_dataset(hdf5_data, structure.H5_IRF_DURATION, (1,), dtype=settings.NEMOH_FLOAT)
        set_hdf5_attributes(dset, structure.H5_IRF_DURATION_ATTR)
        dset[0] = float(x2[2])

    if params.show_pressure is not None:
        dset = require_dataset(hdf5_data, structure.H5_SHOW_PRESSURE, (1,), dtype=settings.NEMOH_INT)
        set_hdf5_attributes(dset, structure.H5_SHOW_PRESSURE_ATTR)
        dset[0] = int(float(x2[0]))

    if params.kochin_function is not None:
        x2 = (' '.join(params.kochin_function)).split()
        dset = require_dataset(hdf5_data, structure.H5_KOCHIN_NUMBER, (1,), dtype=settings.NEMOH_FLOAT)
        set_hdf5_attributes(dset, structure.H5_KOCHIN_NUMBER_ATTR)
        dset[0] = float(x2[0])
        dset = require_dataset(hdf5_data, structure.H5_KOCHIN_MIN, (1,), dtype=settings.NEMOH_FLOAT)
        set_hdf5_attributes(dset, structure.H5_KOCHIN_MIN_ATTR)
        dset[0] = float(x2[1])
        dset = require_dataset(hdf5_data, structure.H5_KOCHIN_MAX, (1,), dtype=settings.NEMOH_FLOAT)
        set_hdf5_attributes(dset, structure.H5_KOCHIN_MAX_ATTR)
        dset[0] = float(x2[2])

    if params.free_surface_elevation:
        x2 = (' '.join(params.free_surface_elevation)).split()
        dset = require_dataset(hdf5_data, structure.H5_FREE_SURFACE_POINTS_X, (1,), dtype=settings.NEMOH_INT)
        set_hdf5_attributes(dset, structure.H5_FREE_SURFACE_POINTS_X_ATTR)
        dset[0] = int(x2[0])
        dset = require_dataset(hdf5_data, structure.H5_FREE_SURFACE_POINTS_Y, (1,), dtype=settings.NEMOH_INT)
        set_hdf5_attributes(dset, structure.H5_FREE_SURFACE_POINTS_Y_ATTR)
        dset[0] = int(x2[1])
        dset = require_dataset(hdf5_data, structure.H5_FREE_SURFACE_DIMENSION_X, (1,), dtype=settings.NEMOH_FLOAT)
        set_hdf5_attributes(dset, structure.H5_FREE_SURFACE_DIMENSION_X_ATTR)
        dset[0] = float(x2[2])
        dset = require_dataset(hdf5_data, structure.H5_FREE_SURFACE_DIMENSION_Y, (1,), dtype=settings.NEMOH_FLOAT)
        set_hdf5_attributes(dset, structure.H5_FREE_SURFACE_DIMENSION_Y_ATTR)
        dset[0] = float(x2[3])



def convert_calculations(filename, hdf5_data):
    """
    Convert nemoh old calculations file to hdf5_data
    Args:
        filename: string, path to the old nemoh calculations file
        hdf5_data: object, an already opened hdf5 file
    """
    x1 = []

    with open(filename, 'r') as inp:
        for line in inp:
            x1.append(line)

    idx = 1
    dset = require_dataset(hdf5_data, structure.H5_ENV_VOLUME, (1,), dtype=settings.NEMOH_FLOAT)
    dset[0] = float(x1[idx].split()[0])
    set_hdf5_attributes(dset, structure.H5_ENV_VOLUME_ATTR)
    idx += 1

    dset = require_dataset(hdf5_data, structure.H5_ENV_GRAVITY, (1,), dtype=settings.NEMOH_FLOAT)
    dset[0] = float(x1[idx].split()[0])
    set_hdf5_attributes(dset, structure.H5_ENV_GRAVITY_ATTR)
    idx += 1


    dset = require_dataset(hdf5_data, structure.H5_ENV_DEPTH, (1,), dtype=settings.NEMOH_FLOAT)
    dset[0] = float(x1[idx].split()[0])
    set_hdf5_attributes(dset, structure.H5_ENV_DEPTH_ATTR)
    idx += 1

    dset = require_dataset(hdf5_data, structure.H5_ENV_WAVE_POINT, (2,), dtype=settings.NEMOH_FLOAT)
    x2 = x1[idx].split()
    dset[0] = float(x2[0])
    dset[1] = float(x2[1])
    set_hdf5_attributes(dset, structure.H5_ENV_WAVE_POINT_ATTR)

    idx = 6

    num_bodies = int(x1[idx].split()[0])

    for i in range(num_bodies):

        body = structure.H5_BODIES + structure.H5_BODY_BASE + str(i+1) + '/'
        idx += 2

        mesh_x = []

        mesh_path = os.path.join(os.path.abspath(os.path.dirname(filename)), str(x1[idx].split()[0]).strip(' \t\n\r'))

        with open(mesh_path, 'r') as mesh_file:
            for line in mesh_file:
                mesh_x.append(line)

        idx += 1
        x2 = x1[idx].split()

        num_points = int(x2[0])
        num_panels = int(x2[1])
        dset = require_dataset(hdf5_data, body + structure.H5_BODY_NUM_POINTS, (1, ), dtype=settings.NEMOH_INT)
        dset[0] = num_points
        set_hdf5_attributes(dset, structure.H5_BODY_NUM_POINTS_ATTR)

        dset = require_dataset(hdf5_data, body + structure.H5_BODY_NUM_PANELS, (1, ), dtype=settings.NEMOH_INT)
        dset[0] = num_panels
        set_hdf5_attributes(dset, structure.H5_BODY_NUM_PANELS_ATTR)

        mesh_idx = 0
        dset = require_dataset(hdf5_data, body + structure.H5_BODY_MESH, (num_points+num_panels+1, 4), dtype=settings.NEMOH_FLOAT)
        mesh_x2 = mesh_x[mesh_idx].split()
        set_hdf5_attributes(dset, structure.H5_BODY_MESH_ATTR)

        dset[0, 0] = int(mesh_x2[0])
        dset[0, 1] = int(mesh_x2[1])

        for j in range(1, num_points+num_panels+1):
            mesh_idx += 1
            mesh_x2 = mesh_x[mesh_idx].split()
            dset[j, :] = [float(x) for x in mesh_x2[:4]]

            if j == num_points:
                mesh_idx += 1

        idx += 1
        num = int(x1[idx].split()[0])
        dset = require_dataset(hdf5_data, body + structure.H5_FREEDOM_DEGREE, (num, 7), dtype=settings.NEMOH_FLOAT)
        set_hdf5_attributes(dset, structure.H5_FREEDOM_DEGREE_ATTR)
        for j in range(num):
            idx += 1
            x2 = x1[idx].split()
            dset[j, :] = np.array([float(x) for x in x2[:7]])

        idx += 1
        num = int(x1[idx].split()[0])
        dset = require_dataset(hdf5_data, body + structure.H5_GENERALISED_FORCES, (num, 7), dtype=settings.NEMOH_FLOAT)
        set_hdf5_attributes(dset, structure.H5_GENERALISED_FORCES_ATTR)
        for j in range(num):
            idx += 1
            x2 = x1[idx].split()
            dset[j, :] = [float(x) for x in x2[:7]]

        idx += 1
        num = int(x1[idx].split()[0])
        for j in range(num):
            idx += 1

    idx += 2
    x2 = x1[idx].split()

    dset = require_dataset(hdf5_data, structure.H5_NUM_WAVE_FREQUENCIES, (1,), dtype=settings.NEMOH_INT)
    set_hdf5_attributes(dset, structure.H5_NUM_WAVE_FREQUENCIES_ATTR)
    dset[0] = int(x2[0])
    dset = require_dataset(hdf5_data, structure.H5_MIN_WAVE_FREQUENCIES, (1,), dtype=settings.NEMOH_FLOAT)
    set_hdf5_attributes(dset, structure.H5_MIN_WAVE_FREQUENCIES_ATTR)
    dset[0] = float(x2[1])
    dset = require_dataset(hdf5_data, structure.H5_MAX_WAVE_FREQUENCIES, (1,), dtype=settings.NEMOH_FLOAT)
    set_hdf5_attributes(dset, structure.H5_MAX_WAVE_FREQUENCIES_ATTR)
    dset[0] = float(x2[2])
    idx += 1
    x2 = x1[idx].split()
    dset = require_dataset(hdf5_data, structure.H5_NUM_WAVE_DIRECTIONS, (1,), dtype=settings.NEMOH_INT)
    set_hdf5_attributes(dset, structure.H5_NUM_WAVE_DIRECTIONS_ATTR)
    dset[0] = int(x2[0])

    dset = require_dataset(hdf5_data, structure.H5_MIN_WAVE_DIRECTIONS, (1,), dtype=settings.NEMOH_FLOAT)
    set_hdf5_attributes(dset, structure.H5_MIN_WAVE_DIRECTIONS_ATTR)
    dset[0] = float(x2[1])

    dset = require_dataset(hdf5_data, structure.H5_MAX_WAVE_DIRECTIONS, (1,), dtype=settings.NEMOH_FLOAT)
    set_hdf5_attributes(dset, structure.H5_MAX_WAVE_DIRECTIONS_ATTR)
    dset[0] = float(x2[2])

    idx += 2
    x2 = x1[idx].split()

    dset = require_dataset(hdf5_data, structure.H5_COMPUTE_IRF, (1,), dtype=settings.NEMOH_INT)
    set_hdf5_attributes(dset, structure.H5_COMPUTE_IRF_ATTR)
    dset[0] = int(x2[0])

    dset = require_dataset(hdf5_data, structure.H5_IRF_TIME_STEP, (1,), dtype=settings.NEMOH_FLOAT)
    set_hdf5_attributes(dset, structure.H5_IRF_TIME_STEP_ATTR)
    dset[0] = float(x2[1])
    dset = require_dataset(hdf5_data, structure.H5_IRF_DURATION, (1,), dtype=settings.NEMOH_FLOAT)
    set_hdf5_attributes(dset, structure.H5_IRF_DURATION_ATTR)
    dset[0] = float(x2[2])

    idx += 1
    x2 = x1[idx].split()
    dset = require_dataset(hdf5_data, structure.H5_SHOW_PRESSURE, (1,), dtype=settings.NEMOH_INT)
    set_hdf5_attributes(dset, structure.H5_SHOW_PRESSURE_ATTR)
    dset[0] = int(x2[0])

    idx += 1
    x2 = x1[idx].split()
    dset = require_dataset(hdf5_data, structure.H5_KOCHIN_NUMBER, (1,), dtype=settings.NEMOH_FLOAT)
    set_hdf5_attributes(dset, structure.H5_KOCHIN_NUMBER_ATTR)
    dset[0] = float(x2[0])
    dset = require_dataset(hdf5_data, structure.H5_KOCHIN_MIN, (1,), dtype=settings.NEMOH_FLOAT)
    set_hdf5_attributes(dset, structure.H5_KOCHIN_MIN_ATTR)
    dset[0] = float(x2[1])
    dset = require_dataset(hdf5_data, structure.H5_KOCHIN_MAX, (1,), dtype=settings.NEMOH_FLOAT)
    set_hdf5_attributes(dset, structure.H5_KOCHIN_MAX_ATTR)
    dset[0] = float(x2[2])


    idx += 1
    x2 = x1[idx].split()
    dset = require_dataset(hdf5_data, structure.H5_FREE_SURFACE_POINTS_X, (1,), dtype=settings.NEMOH_INT)
    set_hdf5_attributes(dset, structure.H5_FREE_SURFACE_POINTS_X_ATTR)
    dset[0] = int(x2[0])
    dset = require_dataset(hdf5_data, structure.H5_FREE_SURFACE_POINTS_Y, (1,), dtype=settings.NEMOH_INT)
    set_hdf5_attributes(dset, structure.H5_FREE_SURFACE_POINTS_Y_ATTR)
    dset[0] = int(x2[1])
    dset = require_dataset(hdf5_data, structure.H5_FREE_SURFACE_DIMENSION_X, (1,), dtype=settings.NEMOH_FLOAT)
    set_hdf5_attributes(dset, structure.H5_FREE_SURFACE_DIMENSION_X_ATTR)
    dset[0] = float(x2[2])
    dset = require_dataset(hdf5_data, structure.H5_FREE_SURFACE_DIMENSION_Y, (1,), dtype=settings.NEMOH_FLOAT)
    set_hdf5_attributes(dset, structure.H5_FREE_SURFACE_DIMENSION_Y_ATTR)
    dset[0] = float(x2[3])

