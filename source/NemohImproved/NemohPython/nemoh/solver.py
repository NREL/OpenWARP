#!/usr/bin/env python
"""
This is the main program for nemoh solver

Changes in version 1.1:
    Added possibility to run the code with custom settings

Changes in version 1.2 (Drift forces and QTF Implementation of Nemoh)
                  Added parameter to store drift forces and yaw moment computed from
                  Nemoh Fortran library

Changes in version 1.3 (Code Acceleration of the Calculation of Influence Coefficients of Nemoh):
    Added switch influence to ode hdf5 settings

Changes in version 1.4 (Implementation of Higher Order Panel Methods):
    Added logic to send handle additionnal Higher order panel variables to be sent to
    Nemoh Fortran
    Disable computation of kochin function, yaw moments and drift forces in the higher Order
    panel method

Changes in version 1.5 (Dipoles Implementation in NEMOH):
    Added logic to send handle additionnal dipoles implementation variables to be sent to
    Nemoh Fortran
    Disable computation of kochin function, yaw moments and drift forces in the dipoles
    implementation

Changes in version 1.6 (Hydrodynamic Data Exporter Assembly v1.0)
       Stored mesh properties like center of buoyancy, volume displacement into hdf5 file.

Changes in version 1.7 (Irregular Frequencies Assembly)
       Added logic to send handle additionnal variables to be sent to Nemoh Fortran.
       Those additional variables determine whether or not Irregular frequencies should be
       removed and if so, the panels which are in the interior of the free surface.
       Applied some bug fixes to allow the shape of hdf5 file dataset 
       to be automatically resized.

Changes in version 1.8 (OpenWarp - Add Logging Functionality)
       Added support for logging.

Changes in version 1.9 (OPENWARP - FIX WAVE FREQUENCY AND DIRECTION CRASH BUG):
    1. Changed the way we do logging from this module when it is run
    as a child process. Using the capturer module now.
"""

import solver_fortran
import structure
import settings
import utility

import numpy as np
import sys
import h5py
import logging
import cStringIO
from contextlib import contextmanager
import contextlib
import os
import StringIO
import tempfile
import json

__author__ = "yedtoss"
__copyright__ = "Copyright (C) 2014-2016 TopCoder Inc. All rights reserved."
__version__ = "1.9"


def init_data():
    """
    Initialize the data to send to the nemoh fortran wrapper
    """

    # No need to log

    d = 2

    n_points = d
    n_panels = d
    n_bodies = d
    n_problems = d
    nbc_panels = d
    n_integration = d
    n_theta = d
    nfs_points = d
    nfs_panels = d
    i_sym = d
    n_potentials = 5*n_points*(1 + (i_sym == 1)) +9*n_panels*(1 + (i_sym == 1))

    data = {
        "rho": 0,
        "g": 0,
        "depth": 0,
        "xeff": 0,
        "yeff": 0,
        "zeff": 0,
        "indiq_solver": 10,
        "max_iterations": 10,
        "restart_param": 10,
        "tol_gmres": 10,
        "i_sym": i_sym,
        "n_panels": n_panels,
        "n_points": n_points,
        "n_bodies": n_bodies,
        "n_problems": n_problems,
        "nbc_panels": nbc_panels,
        "n_integration": n_integration,
        "n_theta": n_theta,
        "nfs_points": nfs_points,
        "nfs_panels": nfs_panels,
        "mesh_p": np.ones((4, n_panels), np.intc, order="F"),
        "mesh_x": np.ones((3, n_points), order="F"),
        "mesh_cpanel": np.ones((n_panels, ), np.intc, order="F"),
        "mesh_xm": np.ones((3, n_panels), order="F"),
        "mesh_n": np.ones((3, n_panels), order="F"),
        "mesh_a": np.ones((n_panels, ), order="F"),
        "bc_normal_velocity": np.ones((nbc_panels, n_problems), np.complex, order="F"),
        "bc_omega": np.ones((n_problems, )),
        "bc_switch_potential": np.ones((n_problems, ), np.intc, order="F"),
        "bc_switch_freesurface": np.ones((n_problems, ), np.intc, order="F"),
        "bc_switch_kochin": np.ones((n_problems, ), np.intc, order="F"),
        "bc_switch_type": np.ones((n_problems, ), np.intc, order="F"),
        "nds": np.ones((n_integration, nbc_panels), order="F"),
        "theta": np.ones((n_theta, ), order="F"),
        "meshfs_p": np.ones((4, nfs_panels), np.intc, order="F"),
        "meshfs_x": np.ones((3, nfs_points), order="F"),
        "out_phi": np.ones((n_problems, 1+ nfs_points), np.complex, order="F"),
        "out_pressure": np.ones((n_problems, nbc_panels), np.complex, order="F"),
        "out_hkochin": np.ones((n_problems, n_theta), np.complex, order="F"),
        "line": np.ones((n_integration, n_problems*2), order="F"),
        "out_potential": np.zeros((n_problems, n_potentials), dtype='f', order="F"),
        "n_potentials": n_potentials,
        "drift_forces": np.ones((n_problems, n_theta, 2), order="F"),
        "yaw_moment": np.ones((n_problems, n_theta), order="F"),
        "fast_influence_switch" : np.zeros((n_problems, ), np.intc, order="F"),
        "use_higher_order": 0

    }

    return data


def write_result(hdf5_data, data):
    """
    Write the result from nemoh fortran to the hdf5
    Args:
        hdf5_data: object the hdf5 opened data
        data: the data sent from nemoh fortran
    """
    signature = __name__ + '.write_result(hdf5_data, data)'
    logger = logging.getLogger(__name__)
    # data is too huge for logging
    utility.log_entrance(logger, signature,
                        {'hdf5_data': hdf5_data})

    dset = utility.require_dataset(hdf5_data, structure.H5_RESULTS_FORCES, data["line"].shape, dtype='f')
    utility.set_hdf5_attributes(dset, structure.H5_RESULTS_FORCES_ATTR)
    dset[:, :] = data["line"].astype(copy=False, dtype='f')

    temp = np.array(data["out_potential"], dtype='f')
    count_skip = 0
    for i in range(data["n_problems"]):
        if data["bc_switch_potential"][i] != 1:
            temp[i, :] = 0
            count_skip += 1
    if count_skip == data["n_problems"]:
        temp = np.zeros((0, 0), dtype='f')

    dset = utility.require_dataset(hdf5_data, structure.H5_RESULTS_POTENTIAL, temp.shape, dtype='f')
    utility.set_hdf5_attributes(dset, structure.H5_RESULTS_POTENTIAL_ATTR)
    dset[:, :] = temp

    kochin = np.zeros((data["n_theta"], 3, data["n_problems"]), dtype='f')
    count_skip = 0

    for i in range(data["n_problems"]):
        if data["bc_switch_kochin"][i] == 1:
            for j in range(data["n_theta"]):

                kochin[j, 0, i] = data["theta"][j]
                kochin[j, 1, i] = np.abs(data["out_hkochin"][i, j])
                kochin[j, 2, i] = np.arctan2(np.imag(data["out_hkochin"][i, j]), np.real(data["out_hkochin"][i, j]))
        else:
            count_skip += 1

    if count_skip ==  data["n_problems"]:
        kochin = np.zeros((0, 0, 0), dtype='f')


    dset = utility.require_dataset(hdf5_data, structure.H5_RESULTS_KOCHIN, kochin.shape, dtype='f')
    utility.set_hdf5_attributes(dset, structure.H5_RESULTS_KOCHIN_ATTR)
    dset[:, :, :] = kochin


    temp = np.zeros((data["nfs_points"], 6, data["n_problems"]), dtype='f')
    count_skip = 0

    for i in range(data["n_problems"]):
        if data["bc_switch_freesurface"][i] == 1:
            for j in range(data["nfs_points"]):

                temp[j, 0, i] = data["meshfs_x"][0, j]
                temp[j, 1, i] = data["meshfs_x"][1, j]
                temp[j, 2, i] = np.abs(data["out_phi"][i, j])
                temp[j, 3, i] = np.arctan2(np.imag(data["out_phi"][i, j]), np.real(data["out_phi"][i, j]))
                temp[j, 4, i] = -np.imag(data["out_phi"][i, j])
                temp[j, 5, i] = -np.real(data["out_phi"][i, j])
        else:
            count_skip += 1

    if count_skip ==  data["n_problems"]:
        temp = np.zeros((0, 0, 0), dtype='f')

    dset = utility.require_dataset(hdf5_data, structure.H5_RESULTS_FREE_SURFACE_POINTS, temp.shape, dtype='f')
    utility.set_hdf5_attributes(dset, structure.H5_RESULTS_FREE_SURFACE_POINTS_ATTR)
    dset[:, :, :] = temp

    temp = np.zeros((data["nfs_panels"], 4, data["n_problems"]), dtype='f')
    count_skip = 0

    for i in range(data["n_problems"]):
        if data["bc_switch_freesurface"][i] == 1:
            for j in range(data["nfs_panels"]):
                temp[j,:, i ] = data["meshfs_p"]
        else:
            count_skip += 1
    if count_skip ==  data["n_problems"]:
        temp = np.zeros((0, 0, 0), dtype='f')

    dset = utility.require_dataset(hdf5_data, structure.H5_RESULTS_FREE_SURFACE_PANEL, temp.shape, dtype='f')
    utility.set_hdf5_attributes(dset, structure.H5_RESULTS_FREE_SURFACE_PANEL_ATTR)
    dset[:, :, :] = temp


    dset = utility.require_dataset(hdf5_data, structure.H5_RESULTS_DRIFT_FORCES, data["drift_forces"].shape, dtype='f')
    utility.set_hdf5_attributes(dset, structure.H5_RESULTS_DRIFT_FORCES_ATTR)
    dset[:, :, :] = data["drift_forces"]

    dset = utility.require_dataset(hdf5_data, structure.H5_RESULTS_YAW_MOMENT, data["yaw_moment"].shape, dtype='f')
    utility.set_hdf5_attributes(dset, structure.H5_RESULTS_YAW_MOMENT_ATTR)
    dset[:, :] = data["yaw_moment"]

    dset = utility.require_dataset(hdf5_data, structure.H5_RESULTS_CENTER_BUOYANCY, data["center_buoyancy"].shape, dtype='f')
    utility.set_hdf5_attributes(dset, structure.H5_RESULTS_CENTER_BUOYANCY_ATTR)
    dset[:, :] = data["center_buoyancy"]

    dset = utility.require_dataset(hdf5_data, structure.H5_RESULTS_VOLUME_DISPLACEMENT, data["displacement"].shape, dtype='f')
    utility.set_hdf5_attributes(dset, structure.H5_RESULTS_VOLUME_DISPLACEMENT_ATTR)
    dset[:] = data["displacement"]

    dset = utility.require_dataset(hdf5_data, structure.H5_RESULTS_WATER_PLANE_AREA, data["waterplane_area"].shape, dtype='f')
    utility.set_hdf5_attributes(dset, structure.H5_RESULTS_WATER_PLANE_AREA_ATTR)
    dset[:] = data["waterplane_area"]

    dset = utility.require_dataset(hdf5_data, structure.H5_RESULTS_STIFNESS, data["stifness"].shape, dtype='f')
    utility.set_hdf5_attributes(dset, structure.H5_RESULTS_STIFNESS_ATTR)
    dset[:, :, :] = data["stifness"]

    utility.log_exit(logger, signature, [None])


def run(hdf5_data):
    """
    Run the solver
    Args:
        hdf5_data: object, the hdf5 opened storage
    Returns:
        the output of the fortran function as string if successful
    """

    signature = __name__ + '.run(hdf5_data)'
    logger = logging.getLogger(__name__)
    utility.log_entrance(logger, signature,
                        {'hdf5_data': hdf5_data})

    data = init_data()
    data["log_level"] = logging.getLogger().getEffectiveLevel()

    dset = utility.get_1d_array(logger, hdf5_data, "H5_L10_COUNT",  expected_dim=4)

    offset = 1

    l10_i_sym = int(dset[0])

    n_points = int(dset[1])
    n_panels = int(dset[2])
    n_bodies = int(dset[3])

    mesh_cpanel = utility.get_dataset(hdf5_data, 'H5_L10_CPANEL')
    mesh_xm = utility.get_dataset(hdf5_data, 'H5_L10_XM')
    mesh_n = utility.get_dataset(hdf5_data, 'H5_L10_N')
    mesh_a = utility.get_dataset(hdf5_data, 'H5_L10_A')


    dset = utility.get_1d_array(logger, hdf5_data, "H5_L12_COUNT",  expected_dim=2)

    i_sym = int(dset[1])

    if l10_i_sym != i_sym or int(dset[0]) != 2:
        raise ValueError('Stopping because the mesh file format is not correct.'
                         'The symmetry about xoz axis is inconsistent')


    data["i_sym"] = i_sym

    data["mesh_p"] = np.asarray(utility.get_dataset(hdf5_data, 'H5_L12_P'), order='F', dtype='i')
    data["mesh_x"] = np.asarray(utility.get_dataset(hdf5_data, 'H5_L12_X'), order='F', dtype='f')

    data["n_points"] = n_points
    data["n_panels"] = n_panels
    data["n_bodies"] = n_bodies

    data["mesh_cpanel"] = np.asarray(mesh_cpanel, order='F', dtype='i')
    data["mesh_xm"] = np.asarray(mesh_xm, order='F', dtype='f')
    data["mesh_n"] = np.asarray(mesh_n, order='F', dtype='f')
    data["mesh_a"] = np.asarray(mesh_a, order='F', dtype='f')

    dset = utility.get_dataset(hdf5_data, 'H5_NORMAL_VELOCITY_W')
    bc_omega = np.asarray(dset, order='F', dtype='f')
    n_problems = bc_omega.shape[0]

    data["bc_omega"] = bc_omega
    data["n_problems"] = n_problems
    dset = utility.get_dataset(hdf5_data, 'H5_NORMAL_VELOCITY_BETA')
    data["bc_switch_type"] = np.asarray(dset, order='F', dtype='i')

    dset = utility.get_dataset(hdf5_data, 'H5_NORMAL_VELOCITY_SWITCH_POTENTIAL')
    data["bc_switch_potential"] = np.asarray(dset, order='F', dtype='i')

    dset = utility.get_dataset(hdf5_data, 'H5_NORMAL_VELOCITY_SWITCH_FREE_SURFACE')
    data["bc_switch_freesurface"] = np.asarray(dset, order='F', dtype='i')

    dset = utility.get_dataset(hdf5_data, 'H5_NORMAL_VELOCITY_SWITCH_KOCHIN')
    data["bc_switch_kochin"] = np.asarray(dset, order='F', dtype='i')

    dset = utility.get_dataset(hdf5_data, 'H5_NORMAL_VELOCITY_VELOCITIES')
    data["bc_normal_velocity"] = np.asarray(dset, order='F', dtype='F')
    data["nbc_panels"] = data["bc_normal_velocity"].shape[0]


    data["rho"] = utility.get_1d_array(logger, hdf5_data, "H5_ENV_VOLUME",  expected_dim=1)[0]
    data["g"] = utility.get_1d_array(logger, hdf5_data, "H5_ENV_GRAVITY",  expected_dim=1)[0]
    data["depth"] = utility.get_1d_array(logger, hdf5_data, "H5_ENV_DEPTH",  expected_dim=1)[0]
    dset = utility.get_1d_array(logger, hdf5_data, "H5_ENV_WAVE_POINT",  expected_dim=2)
    data["xeff"] = dset[0]
    data["y_eff"] = dset[1]

    data["indiq_solver"] = utility.get_1d_array(logger, hdf5_data, "H5_SOLVER_TYPE",  expected_dim=1)[0]
    data["max_iterations"] = utility.get_1d_array(logger, hdf5_data, "H5_SOLVER_GMRES_MAX_ITERATIONS", 1)[0]
    data["restart_param"] = utility.get_1d_array(logger, hdf5_data, "H5_SOLVER_GMRES_RESTART",  expected_dim=1)[0]
    data["tol_gmres"] = utility.get_1d_array(logger, hdf5_data, "H5_SOLVER_GMRES_STOPPING",  expected_dim=1)[0]



    data["nds"] = np.asarray(utility.get_dataset(hdf5_data, 'H5_MESH_INTEGRATION'), order='F', dtype='f')
    data["n_integration"] = data["nds"].shape[0]

    data["use_higher_order"] = utility.get_1d_array(logger, hdf5_data, "H5_SOLVER_USE_HIGHER_ORDER", 1)[0]

    data["num_panel_higher_order"] = utility.get_1d_array(logger, hdf5_data, "H5_SOLVER_NUM_PANEL_HIGHER_ORDER", 1)[0]

    data["b_spline_order"] = utility.get_1d_array(logger, hdf5_data, "H5_SOLVER_B_SPLINE_ORDER",  expected_dim=1)[0]

    data["theta"] = np.asarray(utility.get_dataset(hdf5_data, 'H5_MESH_KOCHIN'), order='F', dtype='f')
    data["n_theta"] = data["theta"].shape[0]


    data["meshfs_x"] = np.asarray(utility.get_2d_array(logger, hdf5_data, "H5_MESH_FREE_SURFACE_VECTORS"), dtype='f')
    data["nfs_points"] = data["meshfs_x"].shape[1]

    data["meshfs_p"] = np.asarray(utility.get_2d_array(logger, hdf5_data, "H5_MESH_FREE_SURFACE_INDEX"),
                                  order='F', dtype='i') + offset
    data["nfs_panels"] = data["meshfs_p"].shape[1]

    data["out_phi"] = np.zeros((n_problems, 1+data["nfs_points"]), dtype='F', order="F")
    data["out_pressure"] = np.zeros((n_problems, data["nbc_panels"]), dtype='F', order="F")
    data["out_hkochin"] = np.zeros((n_problems, data["n_theta"]), dtype='F', order="F")
    data["line"] = np.zeros((data["n_integration"], n_problems*2), order="F", dtype='f')
    data["drift_forces"] = np.zeros((n_problems, data["n_theta"], 2), order="F", dtype='f')
    data["yaw_moment"] = np.zeros((n_problems, data["n_theta"]), order="F", dtype='f')
    data["center_buoyancy"] = np.zeros((n_bodies, 3), order="F", dtype='f')
    data["displacement"] = np.zeros((n_bodies), order="F", dtype='f')
    data["waterplane_area"] = np.zeros((n_bodies), order="F", dtype='f')
    data["stifness"] = np.zeros((n_bodies, 6, 6), order="F", dtype='f')

    n_potentials = 5*n_points*(1 + (i_sym == 1)) +9*n_panels*(1 + (i_sym == 1))
    data["out_potential"] = np.zeros((n_problems, n_potentials), dtype='f', order="F")
    data["n_potentials"] = n_potentials

    data["n_tabulatedx"] = int(utility.get_1d_array(logger, hdf5_data, "H5_SOLVER_GREEN_TABULATION_NUMX", 1)[0])
    data["n_tabulatedz"] = int(utility.get_1d_array(logger, hdf5_data, "H5_SOLVER_GREEN_TABULATION_NUMZ", 1)[0])
    data["n_points_simpson"] = int(utility.get_1d_array(logger, hdf5_data,
                                                        "H5_SOLVER_GREEN_TABULATION_SIMPSON_NPOINTS", 1)[0])

    dset = utility.get_dataset(hdf5_data, 'H5_SOLVER_SWITCH_ODE_INFLUENCE')
    data["fast_influence_switch"] = np.asarray(dset, order='F', dtype='i')

    data["is_interior_domain"] = np.zeros((n_panels), dtype='i', order="F")

    remove_irregular_frequencies = utility.get_1d_array(logger, hdf5_data,
                                                        "H5_SOLVER_REMOVE_IRREGULAR_FREQUENCIES", 1)[0]

    if remove_irregular_frequencies:
        # Bug??? Previous code used dset = hdf5_data.get(structure.H5_SOLVER_REMOVE_IRREGULAR_FREQUENCIES)
        dset = utility.get_dataset(hdf5_data, 'H5_SOLVER_IS_INTERIOR_DOMAIN')
        data["is_interior_domain"] = np.asarray(dset, order='F', dtype='i')


    dset = utility.get_dataset(hdf5_data, 'H5_RESULTS_CASE_BETA')
    data["beta"] = np.asarray(dset, order='F', dtype='f')
    data["n_beta"] = data["beta"].shape[0]

    dset = utility.get_dataset(hdf5_data, 'H5_RESULTS_CASE_RADIATION')
    data["rad_case"] = np.asarray(dset, order='F', dtype='f')
    data["n_radiation"] = data["rad_case"].shape[0]

    dset = utility.get_dataset(hdf5_data, 'H5_SOLVER_THIN_PANELS')
    data["is_thin_body"] = np.asarray(dset, order='F', dtype='i')

    dset = utility.get_dataset(hdf5_data, 'H5_SOLVER_USE_DIPOLES_IMPLEMENTATION')
    data["use_dipoles_implementation"] = dset[0]

    data["remove_irregular_frequencies"] = utility.get_dataset(hdf5_data, 'H5_SOLVER_REMOVE_IRREGULAR_FREQUENCIES')[0]

    dset = utility.get_1d_array(logger, hdf5_data, "H5_SOLVER_COMPUTE_YAW_MOMENT", 1)
    data["compute_yaw_moment"] = dset[0]

    dset = utility.get_1d_array(logger, hdf5_data, "H5_SOLVER_COMPUTE_DRIFT_FORCES", 1)
    data["compute_drift_forces"] = dset[0]

    # Disable kochin, yaw moments and drift forces
    if data["use_higher_order"] == 1 or data["use_dipoles_implementation"] == 1:
        data["n_theta"] = 0
        logger.info('Disabling koching, yaw monment and drift forces computation as '
                    'not supported when higher order panel or dipoles implementation is '
                    'enabled')


    #with CaptureOutput() as capturer:
    solver_fortran.run_solver(data)
    write_result(hdf5_data, data)
    


def solve(custom_config):
    """
    Configure and then run the solver

    Args:
        custom_config, dict The custom configuration dictionary
    Returns:
        the output of the fortran function as string if successful
    """

    signature = __name__ + '.solve(custom_config)'
    logger = logging.getLogger(__name__)
    utility.log_entrance(logger, signature,
                        {'custom_config': custom_config})

    if not custom_config:
        custom_config = {}

    hdf5_file = utility.get_setting(settings.HDF5_FILE, custom_config, 'HDF5_FILE')
    utility.check_is_file(hdf5_file, 'The path to the hdf5 file configured by HDF5_FILE')

    n_tabulatedx = utility.get_setting(settings.GREEN_TABULATION_NUMX, custom_config,
                                       'GREEN_TABULATION_NUMX')

    n_tabulatedz = utility.get_setting(settings.GREEN_TABULATION_NUMZ, custom_config,
                                       'GREEN_TABULATION_NUMZ')

    

    n_points_simpson = utility.get_setting(settings.GREEN_TABULATION_SIMPSON_NPOINTS, custom_config,
                                           'GREEN_TABULATION_SIMPSON_NPOINTS')

    with h5py.File(hdf5_file, "a") as hdf5_db:
        if n_tabulatedx and n_tabulatedx > 0:
            dset = utility.require_dataset(hdf5_db, structure.H5_SOLVER_GREEN_TABULATION_NUMX, (1, ), dtype='i')
            dset[:] = n_tabulatedx

        if n_tabulatedz and n_tabulatedz > 0:
            dset = utility.require_dataset(hdf5_db, structure.H5_SOLVER_GREEN_TABULATION_NUMZ, (1, ), dtype='i')
            dset[:] = n_tabulatedz

        if n_points_simpson and n_points_simpson > 0:
            dset = utility.require_dataset(hdf5_db, structure.H5_SOLVER_GREEN_TABULATION_SIMPSON_NPOINTS, (1, ), dtype='i')
            dset[:] = n_points_simpson

        return run(hdf5_db)


def run_as_process(custom_config, queue):
    utility.setup_subprocess_logging(queue, logging.getLogger())
    return solve(custom_config)

if __name__ == '__main__':
    custom_config = {}
    # Allow custom configuration to be passed from command line
    if len(sys.argv) > 1:
        custom_config = json.loads(sys.argv[1])

    # Allow logging setup to be disabled from command line
    if len(sys.argv) < 3:
        utility.setup_logging(default_conf_path=settings.LOGGING_CONFIGURATION_FILE, logging_path=settings.LOG_FILE)
    try:
        solve({})
    except Exception as e:
        # exc_info=True means the stack trace will be printed automatically
        logging.getLogger(__name__).error('Program halted due to a fatal error whose detail is as follow: ',
                                          exc_info=True)