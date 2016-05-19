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

Changes in version 1.3 (OpenWarp - Add Logging Functionality)
       Added support for logging.

Changes in version 1.4 (OPENWARP - FIX WAVE FREQUENCY AND DIRECTION CRASH BUG):
    1. Removed the class handling the logic to capture the output from a child
       process 
    (taken from http://code.activestate.com/recipes/577564-context-manager-for-low-level-redirection-of-stdou/).
     Using now the capturer module
"""
import numpy as np
import sys
import structure
import settings
import os
import errno
from models import TEnvironment
import inspect
from collections import namedtuple
import json
import logging.config
import h5py
from logutils.queue import QueueHandler
import logging


__author__ = "yedtoss, TCSASSEMBLER"
__copyright__ = "Copyright (C) 2014-2016 TopCoder Inc. All rights reserved."
__version__ = "1.4"


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
    dset = None
    path_exists = path in hdf5_data
    if path_exists:
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
            raise ValueError('Unable to find the wavenumber after ' + str(n_ite) + ' iterations')

        xc = 0.5*(xd+xg)

        while n_ite < n_item_x and np.abs(xd - xg)/ np.abs(xc) >= 1e-6:
            xc=0.5*(xd+xg)
            if (x0 - xg*np.tanh(xg)) * (x0 - xc*np.tanh(xc)) > 0:
                xg = xc
            else:
                xd = xc
            n_ite += 1

        if n_ite >= n_item_x:
            raise ValueError('Unable to find the wavenumber after ' + str(n_ite) + ' iterations')

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


# http://stackoverflow.com/questions/2677185/how-can-i-read-a-functions-signature-including-default-argument-values
DefaultArgSpec = namedtuple('DefaultArgSpec', 'has_default default_value')
def _get_default_arg(args, defaults, arg_index):
    """ Method that determines if an argument has default value or not,
    and if yes what is the default value for the argument

    :param args: array of arguments, eg: ['first_arg', 'second_arg', 'third_arg']
    :param defaults: array of default values, eg: (42, 'something')
    :param arg_index: index of the argument in the argument array for which,
    this function checks if a default value exists or not. And if default value
    exists it would return the default value. Example argument: 1
    :return: Tuple of whether there is a default or not, and if yes the default
    value, eg: for index 2 i.e. for "second_arg" this function returns (True, 42)
    """
    if not defaults:
        return DefaultArgSpec(False, None)

    args_with_no_defaults = len(args) - len(defaults)

    if arg_index < args_with_no_defaults:
        return DefaultArgSpec(False, None)
    else:
        value = defaults[arg_index - args_with_no_defaults]
        if (type(value) is str):
            value = '"%s"' % value
        return DefaultArgSpec(True, value)

def get_method_sig(method):
    """ Given a function, it returns a string that pretty much looks how the
    function signature would be written in python.

    :param method: a python method
    :return: A string similar describing the pythong method signature.
    eg: "my_method(first_argArg, second_arg=42, third_arg='something')"
    """

    # The return value of ArgSpec is a bit weird, as the list of arguments and
    # list of defaults are returned in separate array.
    # eg: ArgSpec(args=['first_arg', 'second_arg', 'third_arg'],
    # varargs=None, keywords=None, defaults=(42, 'something'))
    argspec = inspect.getargspec(method)
    arg_index=0
    args = []

    # Use the args and defaults array returned by argspec and find out
    # which arguments has default
    for arg in argspec.args:
        default_arg = _get_default_arg(argspec.args, argspec.defaults, arg_index)
        if default_arg.has_default:
            args.append("%s=%s" % (arg, default_arg.default_value))
        else:
            args.append(arg)
        arg_index += 1
    return "%s(%s)" % (method.__name__, ", ".join(args))


def log_entrance(logger, signature, parasMap):
    """
    Logs for entrance into public methods at DEBUG level.

    :param logger: the logger object
    :param signature: the method signature
    :param parasMap: the passed parameters
    """
    logger.debug('[Entering method ' + signature + ']')
    if parasMap is not None and len(parasMap.items()) > 0:
        paraStr = '[Input parameters['
        for (k,v) in parasMap.items():
            paraStr += (str(k) + ':' + str(v) + ', ')
        paraStr += ']]'
        logger.debug(paraStr)


def log_exit(logger, signature, parasList):
    """
    Logs for exit from public methods at DEBUG level.

    :param logger: the logger object
    :param signature: the method signature
    :param parasList: the objects to return
    """
    logger.debug('[Exiting method ' + signature + ']')
    if parasList is not None and len(parasList) > 0:
        logger.debug('[Output parameter ' + str(parasList) + ']')


def log_exception(logger, signature, e):
    """
    Logging exception at ERROR level.

    :param logger: the logger object
    :param signature: the method signature
    :param e: the error
    """
    # This will log the traceback.
    logger.error('[Error in method ' + signature + ': Details ' + str(e) + ']', exc_info=True)
    return e


def validate_str(val, allow_none=False, allow_empty=False):
    """
    Check if the given value is a string 
    :param val: the given value to check
    :param allow_none: whether the string is allowed to be None
    :param allow_empty: whether the string is allowed to be empty
    :return False: if val is not of type string and allow_none is False
    :return False: if val is None or empty string
    :return True: otherwise
    """

    if val is None:
        if not allow_none:
            return False
    else:

        if not isinstance(val, str) and not isinstance(val, unicode):
            return False

        elif len(val.strip()) == 0 and not allow_empty:
            return False

    return True


def check_str(val, name, allow_none=False, allow_empty=False):
    """
    Check if the given value is a string 
    :param val: the given value to check
    :param name: name of val
    :param allow_none: whether the string is allowed to be None
    :param allow_empty: whether the string is allowed to be empty
    :raise TypeError: if val is not of type string and allow_none is False
    :raise ValueError: if val is None or empty string
    """

    if val is None:
        if not allow_none:
            raise ValueError(name + ' of value ' + str(val) + ' should not be None.')
    else:

        if not isinstance(val, str) and not isinstance(val, unicode):
            raise TypeError(name + ' of value ' + str(val) + ' should be a string.' + ' but is of type ' + type(val).__name__)

        elif len(val.strip()) == 0 and not allow_empty:
            raise ValueError(name + ' of value ' + str(val) + ' should not empty string.')


def check_type_value(val, name, expected_type, allow_none=False, print_value=True, none_msg=''):
    """
    Check if the given value is of expected type. And also check if the val is None.

    :param val: the given value to check
    :param name: name of val
    :param expected_type: the expected type
    :param allow_none: whether the val is allowed to be None
    :param print_value: whether or not to print the value name in case of error
    :param location: The location of the potential hdf5 value to check
    :raise TypeError: if val is not of expected type
    :raise ValueError: if val is None while not allow None
    """
    message = name

    if print_value:
        message += ' of value ' + str(val)

    if val is None and not allow_none:
        raise ValueError(message + ' should not be None.' + none_msg)
    if not isinstance(val, expected_type):
        raise TypeError(message  + ' should be of type ' + str(expected_type) + '.' + ' but is of type ' + type(val).__name__)

    return val


def check_group_type(val, name='The hdf5 group', allow_none=False, print_value=True, location=''):
    """
    Check if the given value is an hdf5 group. And also check if the val is None.

    :param val: the given value to check
    :param name: name of val
    :param print_value: whether or not to print the value name in case of error
    :param location: The location of the potential hdf5 value to check
    :param allow_none: whether the val is allowed to be None
    :raise TypeError: if val is not of expected type
    :raise ValueError: if val is None while not allow None
    """
    none_msg = name + ' was not found in the hdf5 file at its location ' + location
    return check_type_value(val, name, h5py._hl.group.Group,
                     allow_none=allow_none, print_value=print_value, none_msg=none_msg)


def check_dataset_type(val, name='The hdf5 dataset', allow_none=False, print_value=True, location=''):
    """
    Check if the given value is an hdf5 dataset. And also check if the val is None.
    :param val: the given value to check
    :param name: name of val
    :param print_value: whether or not to print the value name in case of error
    :param location: The location of the potential hdf5 value to check
    :param allow_none: whether the val is allowed to be None
    :raise TypeError: if val is not of expected type
    :raise ValueError: if val is None while not allow None
    """
    none_msg = name + ' was not found in the hdf5 file at its location ' + location
    return check_type_value(val, name, h5py._hl.dataset.Dataset,
                     allow_none=allow_none, print_value=print_value, none_msg=none_msg)


def get_dataset(hdf5_data, path_attribute):
    """
    Check if the array exists in the hdf5 structure and return it
    :param hdf5_data: the hdf5 data
    :param path_attribute: the attribute explaining the role of the array to check
    :return: the array
    :raise ValueError if the array is not found or is not a valid dataset
    """
    path = getattr(structure, path_attribute)
    dset = hdf5_data.get(path)
    default_name = {
        "description": path
    }
    name = str(getattr(structure, path + "_ATTR", default_name)["description"])
    check_dataset_type(dset, name=name, location=path)
    return dset


def get_1d_array(logger, hdf5_data, path_attribute,  expected_dim=-1):
    """
    Get a 1D hdf5 array after checking it exists in the hdf5 structure.
    If the path to the array is not found of is not of expected shape, an error is raised
    Also, if the first dimension of the array does not have the expected number of elements,
    an error is raised
    :param logger: the logger
    :param hdf5_data: the hdf5 data
    :param path_attribute: the attribute containing the description of the role of the array to check
    :param expected_dim: the expected number of elements of the first dimension
    :return: the array if there is no error
    :raise ValueError if the array is not a 1D array or the first dimension has lower number of elements than expected
    """
    path = getattr(structure, path_attribute)
    dset = hdf5_data.get(path)
    default_name = {
        "description": path
    }
    name = str(getattr(structure, path + "_ATTR", default_name)["description"])
    check_dataset_type(dset, name=name, location=path)

    check_array_ndim(dset, name=name, expected_ndim=1)

    if expected_dim > -1:
        check_array_dim(logger, dset, name=name, expected_dim=expected_dim, dim_idx=0)

    return dset

def get_2d_array(logger, hdf5_data, path_attribute,  expected_dim=-1):
    """
    Get a 2D hdf5 array after checking it exists in the hdf5 structure.
    If the path to the array is not found of is not of expected shape, an error is raised
    Also, if the second dimension of the array does not have the expected number of elements,
    an error is raised
    :param logger: the logger
    :param hdf5_data: the hdf5 data
    :param path_attribute: the attribute containing the description of the role of the array to check
    :param expected_dim: the expected number of elements of the second dimension
    :return: the array if there is no error
    :raise ValueError if the array is not a 2D array or the second dimension has lower number of elements than expected
    """
    path = getattr(structure, path_attribute)
    dset = hdf5_data.get(path)
    default_name = {
        "description": path
    }
    name = str(getattr(structure, path + "_ATTR", default_name)["description"])
    check_dataset_type(dset, name=name, location=path)

    check_array_ndim(dset, name=name, expected_ndim=2)

    if expected_dim > -1:
        check_array_dim(logger, dset, name=name, expected_dim=expected_dim, dim_idx=1)

    return dset


def check_value(is_valid, error_msg):
    """
    This function raises an error is is_valid is False
    :param is_valid: whether or not to raise and error
    :param error_msg: the message of the error
    :return: None
    :raise ValueError if is_valid is False
    """
    if not is_valid:
        raise ValueError(error_msg)


def check_array_ndim(arr, name, expected_ndim=2):
    """
    Check that the number of dimensions of a hdf5 array is of expected size
    :param arr: the array to check
    :param name: the name of the array in the hdf5 structure
    :param expected_ndim: the expected number of dimensions
    :return the array
    :raise ValueError if the array is not of the expected number of dimensions
    """
    # We just need to check that the number of dimensions if equal to expected_ndim
    ndim = len(arr.shape)
    check_value(is_valid=(ndim == expected_ndim), error_msg=
                'The number of dimension of ' + name
                + ' is ' + str(ndim) + ' but is expected to be ' + str(expected_ndim))

    return arr


def check_array_dim(logger, arr, name, expected_dim, dim_idx):
    """
    Check that the number of elements of a given dimension of an hdf5 array is of expected size
    Log a warning in case the number of elements in greater than expected.
    If the number of elements is lower than expected an error is raised
    :param logger: the logger
    :param arr: the array to check
    :param name: the name of the array in the hdf5 structure
    :param expected_dim: the expected number of elements
    :param dim_idx: the dimension to check
    :param logger: the logger to use
    :return the array arr if there is no error
    :raise ValueError if the number of elements is lower than expected_dim
    """
    dim = arr.shape[dim_idx]
    # The second  needs to be equal to expected_dim.
    # We raise an error if it less than expected_dim
    dim_msg = str(dim_idx)
    if dim_idx == 0:
        dim_msg = 'The first dimension of '
    elif dim_idx == 1:
        dim_msg = 'The second dimension of '

    check_value(is_valid=(dim >= expected_dim), error_msg=
                dim_msg + name + ' is ' + str(dim) + ' but is expected to be ' + str(expected_dim))
    # If it is greater than expected_dim we warn the user
    if dim > expected_dim:
        logger.warn(dim_msg + name + 'is ' + str(dim) + ' but is expected to be ' + str(expected_dim))

    return arr


def check_array_shape(logger, arr, name, expected_shape):
    """
    This function checks whether or not an hdf5 array is of expected shape
    :param logger: the logger to use
    :param arr: the array to check
    :param name: the name of the array in the hdf5 structure
    :param expected_shape:
    :return: the array arr if there is no error
    """
    shape = arr.shape
    check_array_ndim(arr, name, len(expected_shape))
    for i in range(len(shape)):
        check_array_dim(logger, arr, name, expected_shape[i], i)

    return arr


def check_path_exists(val, name):
    """
    Check if the given value is a legal file or directory path.

    @param val: the given value to check
    @param name: name of val
    @raise ValueError: if val is not a legal file or directory path
    """
    check_str(val, name)
    if not os.path.exists(val):
        raise ValueError(name + ' of value ' + val + '" does not exist.')


def check_is_directory(val, name):
    """
    Check if the given value is a legal directory path.

    @param val: the given value to check
    @param name: name of val
    @raise ValueError: if val is not a legal directory path
    """
    check_path_exists(val, name)
    if not os.path.isdir(val):
        raise ValueError(name + ' of value ' + val + '" is not a legal directory.')


def check_is_file(val, name):
    """
    Check if the given value is a legal file path.

    @param val: the given value to check
    @param name: name of val
    @raise ValueError: if val is not a legal file path
    """
    
    check_path_exists(val, name)
    if not os.path.isfile(val):
        raise ValueError(name + ' of value: ' + val + '" is not a legal file.')


def setup_logging(
    default_conf_path='logging.json', 
    default_level=logging.INFO,
    env_key='LOG_CFG',
    logging_path=None
):
    """
    Setup logging configuration.
    This function reads the logging configuration from a file.
    It then setup the default logger according to the configuration file
    If the file does not exist, we will try to get it from the environment variable denoted by env_key
    :param default_conf_path: the location of the logging configuration file
    :param default_level: the default logging level to use if no configuration file is found
    :param env_key:  the environment key where to retrieve the logging configuration
    :param logging_path: the location where to store the logs
    :return:
    """
    path = default_conf_path
    value = os.getenv(env_key, None)
    if value:
        path = value
    if os.path.exists(path):
        print('Found logging configuration file at ' + default_conf_path + '\n')
        with open(path, 'rt') as f:
            config = json.load(f)

            if logging_path and 'handlers' in config:
                logging_path = os.path.abspath(logging_path)
                print('Writing log at ' + logging_path + '\n')
                mkdir_p(os.path.abspath(os.path.dirname(logging_path)))
                for key, value in config['handlers'].iteritems():
                    if 'filename' in value:
                        value['filename'] = logging_path

        logging.config.dictConfig(config)
    else:
        print('Could not find logging configuration at '+ default_conf_path + '\n')
        print('Using default logging option on console' + '\n')
        logging.basicConfig(level=default_level)

    logging.captureWarnings(capture=True)


def log_and_print(logger, message):
    """
    This function log and prints to the console
    :param logger:  the logger
    :param message:  the message to log
    :return: None
    """
    print(message + '\n')
    logger.info(message)


def get_abs(s):
    """
    This function return the absolute path of a given relative path
    :param s: the path
    :return: the absolute path
    """
    return os.path.abspath(s)


def read_environment(hdf5_data):
    """
    Read the environment from the hdf5 file
    Args:
        hdf5_data: object, the hdf5 opened storage

    Returns:
        The environment
    """

    signature = __name__ + '.read_environment(hdf5_data)'
    logger = logging.getLogger(__name__)
    log_entrance(logger, signature,
                 {})

    environment = TEnvironment()

    dset = hdf5_data.get(structure.H5_ENV_VOLUME)
    check_dataset_type(dset, name=str(structure.H5_ENV_VOLUME_ATTR['description']), location=structure.H5_ENV_VOLUME)
    environment.rho = dset[0]

    dset = hdf5_data.get(structure.H5_ENV_GRAVITY)
    check_dataset_type(dset, name=str(structure.H5_ENV_GRAVITY_ATTR['description']), location=structure.H5_ENV_GRAVITY)
    environment.g = dset[0]

    dset = hdf5_data.get(structure.H5_ENV_DEPTH)
    check_dataset_type(dset, name=str(structure.H5_ENV_DEPTH_ATTR['description']), location=structure.H5_ENV_DEPTH)
    environment.depth = dset[0]

    wave_point = hdf5_data.get(structure.H5_ENV_WAVE_POINT)
    name = str(structure.H5_ENV_WAVE_POINT_ATTR['description'])
    check_dataset_type(wave_point, name=name,
                       location=structure.H5_ENV_WAVE_POINT)
    check_array_ndim(wave_point, name=name,
                     expected_ndim=1)
    check_array_dim(logger, wave_point, name=name, expected_dim=2, dim_idx=0)

    environment.x_eff = wave_point[0]
    environment.y_eff = wave_point[1]
    log_exit(logger, signature, [str(environment)])
    return environment


def touch(fname, times=None):
    """
    This function changes the access time of a file and create it if it does not exists
    :param fname: the path to the file
    :param times: the times to use as access time
    :return: None
    """
    with open(fname, 'a'):
        os.utime(fname, times)

def setup_subprocess_logging(queue, logger):
    """
    This function setup a sub process to log into a queue listened
    by the calling process
    :param queue: the queue to log into
    :param logger: the logger to set up
    """
    # Let's setup a queue handler for the log
    h = QueueHandler(queue) # Just the one handler needed
    logger.handlers = []
    logger.addHandler(h)
    logger.setLevel(logging.DEBUG) # Accepting all logs here, parent process will filter them out
    logging.captureWarnings(capture=True)