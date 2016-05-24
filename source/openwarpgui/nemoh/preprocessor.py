#!/usr/bin/env python
"""
This is the main program for the pre processor. It reads and prepares the Mesh and
calculation cases.  (radiation and diffraction; set of body conditions).

Changes in version 1.1:
    Added possibility to run the code with custom settings

Changes in version 1.2 (Code Acceleration of the Calculation of Influence Coefficients of Nemoh):
    Added switch influence to ode hdf5 settings

Changes in version 1.3 (Implementation of Higher Order Panel Methods):
    Added logic to store USE_HIGHER_ORDER, NUM_PANEL_HIGHER_ORDER, B_SPLINE_ORDER settings
    in hdf5 file.

Changes in version 1.4 (Dipoles Implementation in NEMOH):
    Added logic to store USE_DIPOLES_IMPLEMENTATION and THIN_PANELS settings
    in hdf5 file.

Changes in version 1.5 (Hydrodynamic Data Exporter Assembly v1.0)
       Added parameters controlling wether or not to compute the drift forces or yaw moment

Changes in version 1.6 (Irregular Frequencies Assembly)
       Added logic to discretize the interior of the free surface when the newly added settings
       to remove irregular frequencies is on.
       Applied some bug fixes to allow the shape of hdf5 file dataset 
       to be automatically resized.

Changes in version 1.7 (OpenWarp - Add Logging Functionality)
       Added support for logging.

Changes in version 1.8 (OPENWARP - FIX WAVE FREQUENCY AND DIRECTION CRASH BUG):
    1. Corrected the fact that we were computing the normal velocities beta
    using the wrong shape.

    2. Changed the way we do logging from this module when it is run
    as a child process.
"""

import utility
import numpy as np
import math
import sys
import h5py
import structure

from models import TMesh
from models import TCase
from utility import cih
from utility import sih
import settings
import os
from scipy.spatial import Delaunay
import logging


__author__ = "yedtoss"
__copyright__ = "Copyright (C) 2014-2016 TopCoder Inc. All rights reserved."
__version__ = "1.8"


def read_mesh(hdf5_data, custom_config):
    """
    Read the mesh data from the hdf5 file
    Args:
        hdf5_data: object, the hdf5 opened file

    Return:
        the mesh data
    """
    # Getting the logger here and not globally following recommendation from http://victorlin.me/posts/2012/08/26/good-logging-practice-in-python
    logger = logging.getLogger(__name__)
    signature = __name__ + '.read_mesh(hdf5_data, custom_config)'
    # No need to log the parameter of the method here as it will only be duplicate.
    # This function is never called directly by the user and always call from the preprocess function
    # which already logs the configuration.
    utility.log_entrance(logger, signature,
                        {})
    n_points=0
    n_panels=0
    bodies = hdf5_data.get(structure.H5_BODIES).values()
    n_bodies = len(bodies)

    interior_mesh_points = np.empty((3, 0))
    interior_mesh_panels = np.empty((4, 0))
    interior_c_panels = np.empty((0))
    interior_n_points = 0
    interior_n_panels = 0
    remove_irregular_frequencies = utility.get_setting(settings.REMOVE_IRREGULAR_FREQUENCIES, custom_config,
                                       'REMOVE_IRREGULAR_FREQUENCIES')
    for c in range(n_bodies):
        body = bodies[c]
        dset = body.get(structure.H5_BODY_NUM_POINTS)
        utility.check_dataset_type(dset, name='The number of points for body ' + str(c), location=structure.H5_BODY_NUM_POINTS)
        n_points += dset[0]

        dset = body.get(structure.H5_BODY_NUM_PANELS)
        utility.check_dataset_type(dset, name='The number of panels for body ' + str(c), location=structure.H5_BODY_NUM_PANELS)
        n_panels += dset[0]

    mesh = TMesh(n_points=n_points, n_panels=n_panels, n_bodies=n_bodies)
    logger.info('Found ' + str(n_points) + ' points and '
                + str(n_panels) + ' panels with ' + str(n_bodies) + ' bodies for the mesh')

    n_points = 0
    n_panels = 0

    logger.info('Retrieving the panels coordinates as well as the indices between the body and the problems')
    for c in range(n_bodies):
        body = bodies[c]

        # Both are already checked and logged        
        m = body.get(structure.H5_BODY_NUM_POINTS)[0]
        n = body.get(structure.H5_BODY_NUM_PANELS)[0]

        mesh_arr = body.get(structure.H5_BODY_MESH)
        utility.check_dataset_type(mesh_arr, name='The mesh body ' + str(c), location=structure.H5_BODY_MESH)
        utility.check_array_shape(logger, mesh_arr, name='the mesh array for body number ' + str(c), expected_shape = (m+n+1, 4))

        ns = mesh_arr[0, 1]

        if c > 0 and (ns != mesh.i_sym):
            raise ValueError('There is an inconsistency in the mesh files regarding the xOz symmetries:'
                'The symmetry detected in the body of the first mesh is different from the one on body ' + str(c))
        else:
            mesh.i_sym = int(ns)

        

        for i in range(m):
            mesh.x[:, n_points + i] = np.array(mesh_arr[i + 1, 1:4])

        if remove_irregular_frequencies:
            # If we have to remove frequencies, then we need to discretize the free surface
            int_mesh = generate_mesh(np.asarray(mesh_arr[1:m, 1:4]))
            interior_mesh_points = np.concatenate((interior_mesh_points, int_mesh["x"]), axis=1)
            interior_mesh_panels = np.concatenate((interior_mesh_panels, int_mesh["p"]+mesh.n_points+interior_n_points), axis=1)
            interior_c_panels = np.concatenate((interior_c_panels, c*np.ones(int_mesh["n_panels"])), axis=0)
            interior_n_points += int_mesh["n_points"]
            interior_n_panels += int_mesh["n_panels"]

        for i in range(m, m+n):
            mesh.p[:, n_panels+i-m] = np.array(mesh_arr[i + 1, 0:4]) - 1
            for j in range(4):
                mesh.p[j, n_panels + i-m] += n_points
            mesh.c_panel[n_panels+i-m] = c

        n_points += m
        n_panels += n
        mesh.last_panel[c] = n_panels

    if remove_irregular_frequencies:
        # If we have to remove frequencies, then we need to extend the mesh so
        # that it contains the panels of the free surface too
        mesh_interior = TMesh(n_points=n_points +interior_n_points , n_panels=n_panels + interior_n_panels, n_bodies=n_bodies)
        mesh_interior.x[:, 0:n_points] = mesh.x
        mesh_interior.x[:, n_points:] = interior_mesh_points
        mesh_interior.p[:, 0:n_panels] = mesh.p
        mesh_interior.p[:, n_panels:] = interior_mesh_panels
        mesh_interior.last_panel = mesh.last_panel
        mesh_interior.c_panel[0:n_panels] = mesh.c_panel
        mesh_interior.c_panel[n_panels: ] = interior_c_panels
        mesh_interior.i_sym = mesh.i_sym
        mesh = mesh_interior


        is_interior_domain = np.zeros((n_panels + interior_n_panels))
        is_interior_domain[n_panels:] = 1

        dset = utility.require_dataset(hdf5_data, structure.H5_SOLVER_IS_INTERIOR_DOMAIN, is_interior_domain.shape, dtype='i')
        utility.set_hdf5_attributes(dset, structure.H5_SOLVER_IS_INTERIOR_DOMAIN_ATTR)
        dset[:] = is_interior_domain

        n_panels += interior_n_panels
        n_points += interior_n_points




    logger.info('Computing surface, center of gravity and normal vector of the mesh panels')
    for i in range(mesh.n_panels):
        u = mesh.x[:, mesh.p[1, i]] - mesh.x[:, mesh.p[0, i]]
        v = mesh.x[:, mesh.p[3, i]] - mesh.x[:, mesh.p[1, i]]
        w1 = np.cross(u, v)
        a1 = 0.5*np.linalg.norm(w1)

        u = mesh.x[:, mesh.p[3, i]] - mesh.x[:, mesh.p[2, i]]
        v = mesh.x[:, mesh.p[1, i]] - mesh.x[:, mesh.p[2, i]]
        w2 = np.cross(u, v)
        a2 = 0.5*np.linalg.norm(w2)

        mesh.a[i]= a1+a2

        if mesh.a[i] < utility.EPS:
            raise ValueError('Error: surface of panel ' + str(i) + ' is too small (' + str(mesh.a[i]) + ')')

        mesh.xm[:, i] = (1./3)*(mesh.x[:, mesh.p[0, i]] + mesh.x[:, mesh.p[1, i]] + mesh.x[:, mesh.p[3, i]])*a1/mesh.a[i]

        mesh.xm[:, i] += (1./3)*(mesh.x[:, mesh.p[1, i]] + mesh.x[:, mesh.p[2, i]] + mesh.x[:, mesh.p[3, i]])*a2/mesh.a[i]

        u = w1 + w2

        mesh.n[:, i] = u/np.linalg.norm(u)


    utility.log_exit(logger, signature, [str(mesh)])

    return mesh


def generate_mesh(raw_points):
    """
    Given a list of points corresponding to the discretization of the body domain,
    determine the points belonging to the plan where the body touches the water (or 
        free surface); then use the free surface to generate triangle meshing.

    No debug/info logging here to avoid log pollution as the function is expected to be called many times 
    internally.

    Args:
        raw_points: The 2D array containing the list of points of shape
                    (n_points, 3)

    Return:
        A dictionary containing the points and triangles panel.
    """
    # Get points in the waterline plane. The water plane is z=0 and we allow a tolerance of 1e-3
    points = raw_points[np.abs(raw_points[:, 2]) < 1e-3]

    # Generate a triangle mesh from the waterline segments such that each triangle angle is not
    # too small
    tri_mesh = Delaunay(points[:, 0:2])

    n_panels = tri_mesh.simplices.shape[0]

    # Get the points of the interior of the free surface
    x = points[:, :]

    x[:, 2] = 0

    # Get the meshing connectivity
    p = np.zeros((n_panels, 4))

    p[:, 0:3] = tri_mesh.simplices

    p[:, 3] = tri_mesh.simplices[:, 0]

    return {"n_points" : x.shape[0],
            "n_panels":  n_panels,
            "x": x.transpose(),
            "p": p.transpose()
    }


    


def write_mesh_l12(mesh, hdf5_data):
    """
    Write the l12 data to hdf5 from the mesh
    Args:
        mesh: object, the mesh
        hdf5_data: object, the hdf5 opened file
    """

    # Getting the logger here and not globally following recommendation from http://victorlin.me/posts/2012/08/26/good-logging-practice-in-python
    logger = logging.getLogger(__name__)
    signature = __name__ + '.write_mesh_l12(mesh, hdf5_data)'
    # No need to log the hdf5_data parameter of the method here as it will only be duplicate.
    # This function is never called directly by the user and always call from the preprocess function
    # which already logs the configuration.
    utility.log_entrance(logger, signature,
                        {"mesh" : str(mesh)})

    dset = utility.require_dataset(hdf5_data, structure.H5_L12_COUNT, (2, ), dtype='i')
    utility.set_hdf5_attributes(dset, structure.H5_L12_COUNT_ATTR)
    dset[0] = 2
    dset[1] = int(mesh.i_sym)

    dset = utility.require_dataset(hdf5_data, structure.H5_L12_X, mesh.x.shape, dtype='f')
    utility.set_hdf5_attributes(dset, structure.H5_L12_X_ATTR)
    dset[:, :] = mesh.x



    dset = utility.require_dataset(hdf5_data, structure.H5_L12_P, mesh.p.shape, dtype='i')
    utility.set_hdf5_attributes(dset, structure.H5_L12_P_ATTR)
    dset[:, :] = mesh.p + 1


    utility.log_exit(logger, signature, [None])


def write_mesh_l10(mesh, hdf5_data):
    """
    Write the l10 data to hdf5 from the mesh
    Args:
        mesh: object, the mesh
        hdf5_data: object, the hdf5 opened file
    """

    # Getting the logger here and not globally following recommendation from http://victorlin.me/posts/2012/08/26/good-logging-practice-in-python
    logger = logging.getLogger(__name__)
    signature = __name__ + '.write_mesh_l10(mesh, hdf5_data)'
    # No need to log the hdf5_data parameter of the method here as it will only be duplicate.
    # This function is never called directly by the user and always call from the preprocess function
    # which already logs the configuration.
    utility.log_entrance(logger, signature,
                        {"mesh" : str(mesh)})

    dset = utility.require_dataset(hdf5_data, structure.H5_L10_COUNT, (4, ), dtype='i')
    utility.set_hdf5_attributes(dset, structure.H5_L10_COUNT_ATTR)

    dset[0] = mesh.i_sym
    dset[1] = mesh.n_points
    dset[2] = mesh.n_panels
    dset[3] = mesh.n_bodies

    dset = utility.require_dataset(hdf5_data, structure.H5_L10_CPANEL, mesh.c_panel.shape, dtype='i')
    utility.set_hdf5_attributes(dset, structure.H5_L10_CPANEL_ATTR)
    dset[:] = mesh.c_panel + 1

    dset = utility.require_dataset(hdf5_data, structure.H5_L10_XM, mesh.xm.shape, dtype='f')
    utility.set_hdf5_attributes(dset, structure.H5_L10_XM_ATTR)
    dset[:, :] = mesh.xm

    dset = utility.require_dataset(hdf5_data, structure.H5_L10_N, mesh.n.shape, dtype='f')
    utility.set_hdf5_attributes(dset, structure.H5_L10_N_ATTR)
    dset[:, :] = mesh.n

    dset = utility.require_dataset(hdf5_data, structure.H5_L10_A, mesh.a.shape, dtype='f')
    utility.set_hdf5_attributes(dset, structure.H5_L10_A_ATTR)
    dset[:] = mesh.a
    utility.log_exit(logger, signature, [None])


def write_mesh_tec(mesh, mesh_tec_file):
    """
    Export the mesh to tec file
    Args:
        mesh: object, the mesh
        mesh_tec_file: string, the path to the mesh tec file to save
    """
    # Getting the logger here and not globally following recommendation from http://victorlin.me/posts/2012/08/26/good-logging-practice-in-python
    logger = logging.getLogger(__name__)
    signature = __name__ + '.write_mesh_tec(mesh, mesh_tec_file)'
    utility.log_entrance(logger, signature,
                        {"mesh" : str(mesh),
                        "mesh_tec_file": mesh_tec_file})

    utility.mkdir_p(os.path.abspath(os.path.dirname(mesh_tec_file)))
    logger.info('Converting the mesh file in a tecplot format at ' + str(mesh_tec_file))
    with open(mesh_tec_file, 'w') as inp:
        inp.write('VARIABLES="X" "Y" "Z" "NX" "NY" "NZ" "A"\n')
        inp.write('ZONE N=\t' + str(mesh.n_points) + '\t, E=\t' + str(mesh.n_panels) + '\t, F=FEPOINT,ET=QUADRILATERAL\n')
        for i in range(mesh.n_points):
            s = str(mesh.x[0, i]) + '\t' + str(mesh.x[1, i]) + '\t' + str(mesh.x[2, i]) + '\t0.\t0.\t0.\t0.\n'
            inp.write(s)

        for i in range(mesh.n_panels):
            s = str(mesh.p[0, i] + 1) + '\t' + str(mesh.p[1, i] + 1) + '\t' + str(mesh.p[2, i] + 1) + '\t' + str(mesh.p[3, i] + 1) + '\n'
            inp.write(s)

        inp.write('ZONE t="normales", F=POINT, I=\t' + str(mesh.n_panels) + '\n')

        for i in range(mesh.n_panels):
            s = str(mesh.xm[0, i]) + '\t' + str(mesh.xm[1, i]) + '\t' + str(mesh.xm[2, i]) + '\t'
            s += str(mesh.n[0, i]) + '\t' + str(mesh.n[1, i]) + '\t' + str(mesh.n[2, i]) + '\t'
            s += str(mesh.a[i]) + '\n'
            inp.write(s)

    utility.log_and_print(logger, 'The mesh is converted to tec format in '
                          + utility.get_abs(mesh_tec_file))
    utility.log_exit(logger, signature, [None])


def write_fk_force_tec(int_case, fk_force, w, beta, filename):
    """
    Writes the froude krylov forces to .tec format
    Args:
        int_case: 1D array, the integration cases
        fk_forces: 3D array, the froudkrylov forces
        w: 1D array, represents the wave frequencies omega
        beta: 1D array, represents the wave directions beta
        filename: string, the path to the file where to save the forces
    """
    # Getting the logger here and not globally following recommendation from http://victorlin.me/posts/2012/08/26/good-logging-practice-in-python
    logger = logging.getLogger(__name__)
    signature = __name__ + '.write_fk_force_tec(int_case, fk_force, w, beta, filename)'
    utility.log_entrance(logger, signature,
                        {"int_case" : str(int_case),
                        "filename": filename})

    utility.mkdir_p(os.path.abspath(os.path.dirname(filename)))
    n_integration = len(int_case)
    n_beta = len(beta)
    n_w = len(w)
    logger.info('Converting the froude krylov forces in a tecplot format at ' + str(filename))
    with open(filename, 'w') as inp:
        inp.write('VARIABLES="w (rad/s)"\n')
        for k in range(n_integration):
            s = '"abs(F\t' + str(int_case[k].body + 1) + '\t' + str(k+1) + ')" "angle(F\t'
            s += str(int_case[k].body + 1) + '\t' + str(k+1) + ')"\n'
            inp.write(s)

        for c in range(n_beta):
            inp.write('Zone t="FKforce - beta =\t' + str(beta[c]*180./np.pi) + '",I=\t' + str(n_w) + ',F=POINT\n')
            for i in range(n_w):
                s = str(w[i]) + '\t'
                for k in range(n_integration):
                    val = str(np.arctan2(np.imag(fk_force[i, c, k]), np.real(fk_force[i, c, k])))
                    s += str(np.abs(fk_force[i, c, k])) + '\t' + val + '\t'
                inp.write(s)
                inp.write('\n')

    utility.log_and_print(logger, 'The fk forces were converted to tec format in '
                          + utility.get_abs(filename))

    utility.log_exit(logger, signature, [None])


def compute_nds(mesh, c, i_case, direction, axis):
    """
    Compute the integration nds
    Args:
        mesh: object The mesh
        c: int, the panel index
        i_case: int, the integration case
        direction 1D array of length 3: The direction (x, y or z)
        axis 1D array of length 3:  The axis coordinate

    No debug/info logging here to avoid log pollution as the function is expected to be called many times 
        internally.

    Returns:
        the integration array nds
    """
    nds = np.zeros(mesh.n_panels*2**mesh.i_sym, settings.NEMOH_FLOAT)
    vel = np.copy(direction[0:3])
    if i_case == 1:
        for i in range(mesh.n_panels):
            if mesh.c_panel[i] == c:
                #vel = np.copy(direction[0:3])
                nds[i] = - mesh.a[i] * (mesh.n[0, i] *vel[0] + mesh.n[1, i] *vel[1] + mesh.n[2, i] *vel[2])
            else:
                nds[i]=0.

            if mesh.i_sym == 1:
                if mesh.c_panel[i] == c:
                    #vel = np.copy(direction[0:3])
                    nds[i+ mesh.n_panels] = -mesh.a[i]*(mesh.n[0, i]*vel[0]-mesh.n[1, i]*vel[1] + mesh.n[2, i]*vel[2])
                else:
                    nds[i+ mesh.n_panels] = 0.

    elif i_case == 2:
        for i in range(mesh.n_panels):
            if mesh.c_panel[i] == c:
                vel[0] = direction[1]*(mesh.xm[2, i] - axis[2]) - direction[2]*(mesh.xm[1, i] - axis[1])
                vel[1] = direction[2]*(mesh.xm[0, i] - axis[0]) - direction[0]*(mesh.xm[2, i] - axis[2])
                vel[2] = direction[0]*(mesh.xm[1, i] - axis[1]) - direction[1]*(mesh.xm[0, i] - axis[0])
                nds[i] = - mesh.a[i] * (mesh.n[0, i] *vel[0] + mesh.n[1, i] *vel[1] + mesh.n[2, i] *vel[2])
            else:
                nds[i]=0.

            if mesh.i_sym == 1:
                if mesh.c_panel[i] == c:
                    vel[0] = direction[1]*(mesh.xm[2, i] - axis[2]) - direction[2]*(-mesh.xm[1, i] - axis[1])
                    vel[1] = direction[2]*(mesh.xm[0, i] - axis[0]) - direction[0]*(mesh.xm[2, i] - axis[2])
                    vel[2] = direction[0]*(-mesh.xm[1, i] - axis[1]) - direction[1]*(mesh.xm[0, i] - axis[0])
                    nds[i+ mesh.n_panels] = -mesh.a[i]*(mesh.n[0, i]*vel[0] - mesh.n[1, i]*vel[1] + mesh.n[2, i]*vel[2])
                else:
                    nds[i+ mesh.n_panels] = 0.

    elif i_case == 3:
        raise NotImplementedError('Force case 3 is not implemented yet')
    else:
        raise RuntimeError('The radiation case of index ' + str(i_case) + ' is unknown')

    return nds


def compute_radiation_condition(mesh, c, i_case,  direction, axis):
    """
    Compute the radiation condition
    Args:
        mesh: object The mesh
        c: int, the panel index
        i_case: int, the integration case
        direction 1D array of length 3: The direction (x, y or z)
        axis 1D array of length 3:  The axis coordinate

    No debug/info logging here to avoid log pollution as the function is expected to be called many times 
        internally.

    Returns:
        the radiation condition array n_vel
    """
    n_vel = np.zeros(mesh.n_panels*2**mesh.i_sym, settings.NEMOH_COMPLEX)
    vel = np.copy(direction[0:3])
    if i_case == 1:
        for i in range(mesh.n_panels):
            if mesh.c_panel[i] == c:
                vel = np.copy(direction[0:3])
                n_vel[i] = complex(np.sum(mesh.n[:, i].flatten()*vel.flatten()), 0)
            else:
                n_vel[i] = complex(0, 0)

            if mesh.i_sym == 1:
                if mesh.c_panel[i] == c:
                    vel = np.copy(direction[0:3])
                    #nn = mesh.n[:, i]
                    #nn[1] *= -1
                    #n_vel[i + mesh.n_panels] = complex(np.sum(nn.flatten()*vel.flatten()), 0)
                    n_vel[i + mesh.n_panels] = complex(mesh.n[0,i]*vel[0]-mesh.n[1,i]*vel[1]+ mesh.n[2,i]*vel[2], 0)
                else:
                    n_vel[i+ mesh.n_panels] = complex(0, 0)

    elif i_case == 2:
        for i in range(mesh.n_panels):
            if mesh.c_panel[i] == c:
                vel[0] = direction[1]*(mesh.xm[2, i] - axis[2]) - direction[2]*(mesh.xm[1, i] - axis[1])
                vel[1] = direction[2]*(mesh.xm[0, i] - axis[0]) - direction[0]*(mesh.xm[2, i] - axis[2])
                vel[2] = direction[0]*(mesh.xm[1, i] - axis[1]) - direction[1]*(mesh.xm[0, i] - axis[0])
                n_vel[i] = complex(np.sum(mesh.n[:, i].flatten()*vel.flatten()), 0)
            else:
                n_vel[i] = complex(0, 0)

            if mesh.i_sym == 1:
                if mesh.c_panel[i] == c:
                    vel[0] = direction[1]*(mesh.xm[2, i] - axis[2]) - direction[2]*(-mesh.xm[1, i] - axis[1])
                    vel[1] = direction[2]*(mesh.xm[0, i] - axis[0]) - direction[0]*(mesh.xm[2, i] - axis[2])
                    vel[2] = direction[0]*(-mesh.xm[1, i] - axis[1]) - direction[1]*(mesh.xm[0, i] - axis[0])
                    #nn = mesh.n[:, i]
                    #nn[1] *= -1
                    #n_vel[i+ mesh.n_panels] = complex(np.sum(nn.flatten()*vel.flatten()), 0)
                    n_vel[i + mesh.n_panels] = complex(mesh.n[0,i]*vel[0]-mesh.n[1,i]*vel[1]+ mesh.n[2,i]*vel[2], 0)
                else:
                    n_vel[i + mesh.n_panels] = complex(0, 0)

    elif i_case == 3:
        raise NotImplementedError('Force case 3 is not implemented yet')
    else:
        raise RuntimeError('The radiation case of index ' + str(i_case) + ' is unknown')

    return n_vel


def compute_one_wave(k, w, beta, wt, environment):
    """
    Calculate the complex potential, pressure and fluid velocities for a regular wave eta=sin(k*wbar-wt)
    Args:
        k: float, the wave number
        w: float, the wave frequency
        beta: float, the wave direction
        wt: 1D array of length 3, the wave position
        environment: object, the environment

    No debug/info logging here to avoid log pollution as the function is expected to be called many times 
        internally.

    Returns
        A dictionary containing the potential, pressure and fluid velocities
    """

    x = wt[0]
    y = wt[1]
    z = wt[2]

    w_bar = (x-environment.x_eff)*np.cos(beta)+(y-environment.y_eff)*np.sin(beta)
    phi = -environment.g/w*cih(k, z, environment.depth)*np.exp(utility.II*k*w_bar)
    p = -environment.rho*environment.g*utility.II*cih(k, z, environment.depth)*np.exp(utility.II*k*w_bar)
    vx = -environment.g/w*utility.II*k*np.cos(beta)*cih(k, z, environment.depth)*np.exp(utility.II*k*w_bar)
    vy = -environment.g/w*utility.II*k*np.sin(beta)*cih(k, z, environment.depth)*np.exp(utility.II*k*w_bar)
    vz = -environment.g/w*k*sih(k, z, environment.depth)*np.exp(utility.II*k*w_bar)

    return {"phi": phi, "p": p, "vx": vx, "vy": vy, "vz": vz, "v": np.array([vx, vy, vz])}


def compute_wave(mesh, w, beta, environment):
    """
    Calculate the array of complex potential, pressure and fluid velocities for a wave
    Args:
        mesh: object, the mesh
        w: float, the wave frequency
        beta: float, the wave direction
        environment: object, the environment

    No debug/info logging here to avoid log pollution as the function is expected to be called many times 
        internally.
    Returns
        A dictionary containing the pressure and fluid velocities
    """

    n_vel = np.zeros(mesh.n_panels*2**mesh.i_sym, settings.NEMOH_COMPLEX)
    pressure = np.zeros(mesh.n_panels*2**mesh.i_sym, settings.NEMOH_COMPLEX)

    k_wave = utility.compute_wave_number(w, environment)

    for i in range(2**mesh.i_sym*mesh.n_panels):
        if i < mesh.n_panels:
            wbar = (mesh.xm[0, i] - environment.x_eff)*np.cos(beta) + (mesh.xm[1, i] - environment.y_eff)*np.sin(beta)

            pressure[i] = -environment.g/w*np.exp(utility.II*k_wave*wbar)

            n_vel[i] = pressure[i]*(utility.II*k_wave*(np.cos(beta)*mesh.n[0,i]+ \
                    np.sin(beta)*mesh.n[1,i])*cih(k_wave,mesh.xm[2, i], environment.depth)+ \
                    k_wave*mesh.n[2,i]*sih(k_wave,mesh.xm[2,i],environment.depth))
            pressure[i] *= cih(k_wave,mesh.xm[2,i], environment.depth)
            one_wave = compute_one_wave(k_wave, w, beta, mesh.xm[0:3, i], environment)
            # This makes previous pressure[i] statement useless
            pressure[i] = one_wave["p"]
            n_vel[i] = np.sum(one_wave["v"].flatten()*mesh.n[:, i].flatten())

        else:
            wbar=(mesh.xm[0, i-mesh.n_panels]-environment.x_eff)*np.cos(beta)+ \
                    (-mesh.xm[1, i-mesh.n_panels]-environment.y_eff)*np.sin(beta)
            pressure[i] = -environment.g/w*np.exp(utility.II*k_wave*wbar)
            n_vel[i] = pressure[i]*(utility.II*k_wave*(np.cos(beta)*mesh.n[0,i-mesh.n_panels]+\
                    np.sin(beta)*(-1.*mesh.n[1,i-mesh.n_panels]))*cih(k_wave, mesh.xm[2, i-mesh.n_panels],\
                    environment.depth)+k_wave*mesh.n[2,i-mesh.n_panels]*sih(k_wave,mesh.xm[2,i-mesh.n_panels],\
                    environment.depth))
            pressure[i] *= cih(k_wave, mesh.xm[2, i-mesh.n_panels], environment.depth)
            #one_wave = compute_one_wave(k_wave, w, beta, mesh.xm[0:3, i-mesh.n_panels], environment)
            #xm = mesh.xm[0:3, i-mesh.n_panels]
            #xm[1] = - xm[1]

            xm = np.array([mesh.xm[0, i-mesh.n_panels], -mesh.xm[1, i-mesh.n_panels], mesh.xm[2, i-mesh.n_panels]])
            one_wave = compute_one_wave(k_wave, w, beta, xm, environment)
            pressure[i] = one_wave["p"]
            #one_wave["v"][1] *= -1
            #n_vel[i] = np.sum(one_wave["v"]*mesh.n[:, i-mesh.n_panels])
            vx = one_wave["v"][0]
            vy = one_wave["v"][1]
            vz = one_wave["v"][2]
            n_vel[i] = vx*mesh.n[0, i-mesh.n_panels] - vy*mesh.n[1, i-mesh.n_panels] + vz*mesh.n[2, i-mesh.n_panels]

    return {"n_vel": n_vel, "pressure": pressure}


def run(hdf5_data, custom_config):
    """
    This function run the preprocessor
    Args:
        hdf5_data: object, the hdf5 opened file
        custom_config, dict The custom configuration dictionary
    """

    # Getting the logger here and not globally following recommendation from http://victorlin.me/posts/2012/08/26/good-logging-practice-in-python
    logger = logging.getLogger(__name__)
    signature = __name__ + '.run(hdf5_data, custom_config)'
    # No need to log the parameter of the method here as it will only be duplicate.
    # This function is never called directly by the user and always call from the preprocess function
    # which already logs the configuration.
    utility.log_entrance(logger, signature,
                        {})

    n_radiation = 0
    n_integration = 0

    bodies = hdf5_data.get(structure.H5_BODIES)
    utility.check_group_type(bodies, name='The bodies group', location=structure.H5_BODIES)

    logger.info("Processing the bodies group: " + str(bodies))
    bodies = bodies.values()
    for body in bodies:
        utility.check_group_type(body, name='The sub-body group', location=structure.H5_BODIES)
        dset = body.get(structure.H5_FREEDOM_DEGREE)
        utility.check_dataset_type(dset, name='The freedom degree', location=structure.H5_FREEDOM_DEGREE)
        n_radiation += dset.shape[0]

        dset = body.get(structure.H5_GENERALISED_FORCES)
        utility.check_dataset_type(dset, name='The generalised forces', location=structure.H5_GENERALISED_FORCES)
        n_integration += dset.shape[0]

    logger.info("Solving " +  str(n_radiation) + " problems  with " + str(n_integration) + " forces")

    logger.info("Processing wave frequencies information")

    dset = hdf5_data.get(structure.H5_NUM_WAVE_FREQUENCIES)
    utility.check_dataset_type(dset, name='The number of wave frequencies', location=structure.H5_NUM_WAVE_FREQUENCIES)
    n_w = dset[0]

    dset = hdf5_data.get(structure.H5_MIN_WAVE_FREQUENCIES)
    utility.check_dataset_type(dset, name='The minimum wave frequency', location=structure.H5_MIN_WAVE_FREQUENCIES)
    w_min = dset[0]

    dset = hdf5_data.get(structure.H5_MAX_WAVE_FREQUENCIES)
    utility.check_dataset_type(dset, name='The maximum wave frequency', location=structure.H5_MAX_WAVE_FREQUENCIES)
    w_max = dset[0]

    w = np.zeros(n_w, settings.NEMOH_FLOAT)
    if n_w > 1:
        for j in range(n_w):
            w[j] = w_min+(w_max-w_min)*j/(n_w-1)
    else:
        w[0] = w_min
    logger.info(' Using ' + str(n_w) + ' equally spaced wave frequencies from ' + str(w[0]) + ' to ' + str(w[n_w-1]))


    logger.info("Processing wave directions information")

    dset = hdf5_data.get(structure.H5_NUM_WAVE_DIRECTIONS)
    utility.check_dataset_type(dset, name='The number of wave directions', location=structure.H5_NUM_WAVE_DIRECTIONS)
    n_beta = dset[0]

    dset = hdf5_data.get(structure.H5_MIN_WAVE_DIRECTIONS)
    utility.check_dataset_type(dset, name='The minimum wave direction', location=structure.H5_MIN_WAVE_DIRECTIONS)
    beta_min = dset[0]

    dset = hdf5_data.get(structure.H5_MAX_WAVE_DIRECTIONS)
    utility.check_dataset_type(dset, name='The maximum wave direction', location=structure.H5_MAX_WAVE_DIRECTIONS)
    beta_max = dset[0]

    beta = np.zeros(n_beta, settings.NEMOH_FLOAT)
    if n_beta > 1:
        for j in range(n_beta):
            beta[j] = (beta_min+(beta_max-beta_min)*j/(n_beta-1))*math.pi/180.
    else:
        beta[0] = beta_min * math.pi/180.
    logger.info(' Using ' + str(n_beta) + str(' equally spaced wave directions from  ') + str(beta[0]) + ' to ' + str(beta[n_beta-1]))


    dset = hdf5_data.get(structure.H5_SHOW_PRESSURE)
    utility.check_dataset_type(dset, name='The switch for showing pressure', location=structure.H5_SHOW_PRESSURE)
    switch_potential = dset[0] >= 1

    dset = hdf5_data.get(structure.H5_KOCHIN_NUMBER)
    utility.check_dataset_type(dset, name='The number of direction for the computation of far field coefficients (Kochin function)',
     location=structure.H5_KOCHIN_NUMBER)
    n_theta = dset[0]

    dset = hdf5_data.get(structure.H5_KOCHIN_MIN)
    utility.check_dataset_type(dset, name='The minimum number of direction for the computation of far field coefficients (Kochin function)',
     location=structure.H5_KOCHIN_MIN)
    theta_min = dset[0]

    dset = hdf5_data.get(structure.H5_KOCHIN_MAX)
    utility.check_dataset_type(dset, name='The maximum  number of direction for the computation of far field coefficients (Kochin function)',
     location=structure.H5_KOCHIN_MAX)
    theta_max = dset[0]

    switch_kochin = n_theta > 0

    if switch_kochin:
        logger.info('The computation of far field coefficients (Kochin function) is enabled with '
         + str(n_theta) + ' directions from ' + str(theta_min) + ' to ' + str(theta_max))
    else:
        logger.info('The computation of far field coefficients (Kochin function) is disabled')

    dset = hdf5_data.get(structure.H5_FREE_SURFACE_POINTS_X)
    utility.check_dataset_type(dset, name=str(structure.H5_FREE_SURFACE_POINTS_X_ATTR['description']),
                               location=structure.H5_FREE_SURFACE_POINTS_X)
    n_x = dset[0]

    dset = hdf5_data.get(structure.H5_FREE_SURFACE_POINTS_Y)
    utility.check_dataset_type(dset, name=str(structure.H5_FREE_SURFACE_POINTS_Y_ATTR['description']),
                               location=structure.H5_FREE_SURFACE_POINTS_Y)
    n_y = dset[0]

    dset = hdf5_data.get(structure.H5_FREE_SURFACE_DIMENSION_X)
    utility.check_dataset_type(dset, name=str(structure.H5_FREE_SURFACE_DIMENSION_X_ATTR['description']),
                               location=structure.H5_FREE_SURFACE_DIMENSION_X)
    l_x = dset[0]

    dset = hdf5_data.get(structure.H5_FREE_SURFACE_DIMENSION_Y)
    utility.check_dataset_type(dset, name=str(structure.H5_FREE_SURFACE_DIMENSION_Y_ATTR['description']),
     location=structure.H5_FREE_SURFACE_DIMENSION_Y)
    l_y = dset[0]

    switch_free_surface = n_x > 0

    if switch_free_surface:
        logger.info('The computation of the wave elevation is enabled with '
         + str(n_x) + ' points in x direction of dimension' + str(l_x) 
         + ' and ' + str(n_y) + ' points in y direction of dimension ' + str(l_y))
    else:
        logger.info('The computation of the wave elevation is disabled')

    rad_case = [TCase() for x in range(n_radiation)]
    int_case = [TCase() for x in range(n_integration)]
    j_rad = 0
    j_int = 0

    for c in range(len(bodies)):
        body = bodies[c]
        # This has already been checked
        freedom_degree = body.get(structure.H5_FREEDOM_DEGREE)
        utility.check_array_ndim(freedom_degree, name='the freedom degree for body number ' + str(c),
            expected_ndim=2)
        utility.check_array_dim(logger, freedom_degree, name='the freedom degree for body number ' + str(c), expected_dim=7, dim_idx=1)

        m = freedom_degree.len()

        

        for i in range(m):
            case = TCase()
            case.i_case = freedom_degree[i, 0]
            case.direction = np.array(freedom_degree[i, 1:4])
            case.axis = np.array(freedom_degree[i, 4:7])
            case.body = c
            case.mode = i
            rad_case[j_rad + i] = case
        j_rad += m

        # This has already been checked
        generalised_forces = body.get(structure.H5_GENERALISED_FORCES)
        utility.check_array_ndim(generalised_forces, name='the generalised forces for body number ' + str(c),
            expected_ndim=2)


        m = generalised_forces.len()

        utility.check_array_dim(logger, generalised_forces, name='generalised forces for body number ' + str(c), expected_dim=7, dim_idx=1)

        for i in range(m):
            case = TCase()
            case.i_case = generalised_forces[i, 0]
            case.direction = np.array(generalised_forces[i, 1:4])
            case.axis = np.array(generalised_forces[i, 4:7])
            case.body = c
            case.mode = i
            int_case[j_int + i] = case

        j_int += m

    utility.log_and_print(logger, 'Summary of calculation')

    dset = hdf5_data.get(structure.H5_ENV_DEPTH)
    utility.check_dataset_type(dset, name='The depth of fluid (water)',
     location=structure.H5_ENV_DEPTH)
    depth = dset[0]
    if depth > 0:
        utility.log_and_print(logger, ' The water depth is ' + str(depth) + ' meter')
    else:
        utility.log_and_print(logger, ' The water depth is infinite ')

    utility.log_and_print(logger, '  -> ' + str(n_w) + ' wave frequencies from ' + str(w[0]) + ' to ' + str(w[n_w-1]))
    utility.log_and_print(logger, '  -> ' + str(n_beta) + str(' wave directions from  ')
                          + str(beta[0]) + ' to ' + str(beta[n_beta-1]))
    utility.log_and_print(logger, '  -> ' + str(n_radiation) + ' radiation problems')
    utility.log_and_print(logger, '  -> ' + str(n_integration) + ' forces')

    logger.info('Reading the mesh from the hdf5 file')
    mesh = read_mesh(hdf5_data, custom_config)
    write_mesh_l12(mesh, hdf5_data)
    write_mesh_l10(mesh, hdf5_data)

    mesh_tec_file = utility.get_setting(settings.MESH_TEC_FILE, custom_config, 'MESH_TEC_FILE')

    if mesh_tec_file:
        write_mesh_tec(mesh, mesh_tec_file)

    fnds = np.zeros((n_integration, mesh.n_panels*2**mesh.i_sym), settings.NEMOH_FLOAT)

    logger.info('Computing the integration nds')
    for j in range(n_integration):
        fnds[j, :] = compute_nds(mesh, int_case[j].body, int_case[j].i_case, int_case[j].direction, int_case[j].axis)


    dset = utility.require_dataset(hdf5_data, structure.H5_MESH_INTEGRATION, fnds.shape, dtype='f')
    utility.set_hdf5_attributes(dset, structure.H5_MESH_INTEGRATION_ATTR)
    dset[:, :] = fnds
    logger.info('The information about how the pressure has to be integrated ' 
        'over the body surfaces to obtain the requested forces has been saved in '
        + str(structure.H5_MESH_INTEGRATION) + '. Its characteristic is: '
        + str(dset))

    environment = utility.read_environment(hdf5_data)

    normal_velocity = np.zeros((mesh.n_panels*2**mesh.i_sym, (n_beta+n_radiation)*n_w), settings.NEMOH_COMPLEX)
    fk_force = np.zeros((n_w, n_beta, n_integration), settings.NEMOH_COMPLEX)

    logger.info('Computing the body conditions for each radiation and diffraction problem (normal velocities)'
        ' and the Froude Krylov forces for each of the diffraction problem (FK forces)')
    for i in range(n_w):
        for j in range(n_beta):

            result = compute_wave(mesh, w[i], beta[j], environment)
            pressure = result["pressure"]
            n_vel = result["n_vel"]
            normal_velocity[:, j+ i*(n_beta+n_radiation)] = n_vel
            # Calculate the corresponding FK forces
            for k in range(n_integration):
                #for c in range(mesh.n_panels*2**mesh.i_sym):
                    #fk_force[i, j, k] +=  pressure[c]*fnds[k, c]

                fk_force[i, j, k] = np.sum(pressure.flatten()*fnds[k, :].flatten())

        for j in range(n_radiation):
            n_vel = compute_radiation_condition(mesh, rad_case[j].body, rad_case[j].i_case, rad_case[j].direction,
                                        rad_case[j].axis)

            normal_velocity[:, j + n_beta + i*(n_beta+n_radiation)] = n_vel

    # Save body conditions
    n_problems = n_w*(n_radiation+n_beta)
    bc_omega = w.repeat(n_beta + n_radiation)
    dset = utility.require_dataset(hdf5_data, structure.H5_NORMAL_VELOCITY_W, bc_omega.shape, dtype='f', maxshape=(None))
    utility.set_hdf5_attributes(dset, structure.H5_NORMAL_VELOCITY_W_ATTR)
    dset[:] = bc_omega
    

    bc_switch_type = -np.ones(n_problems, dtype='f') # Set the whole array to -1
    # Set the first n_beta values to beta, skip the next n_radiation and so on until the end of the array
    for i in xrange(0, n_problems, n_beta + n_radiation):
        bc_switch_type[i:i + n_beta] = beta
    dset = utility.require_dataset(hdf5_data, structure.H5_NORMAL_VELOCITY_BETA, bc_switch_type.shape, dtype='f')
    utility.set_hdf5_attributes(dset, structure.H5_NORMAL_VELOCITY_BETA_ATTR)
    dset[:] = bc_switch_type

    logger.info('Saved the wave frequencies and directions for each of the ' + str(n_beta + n_radiation) 
        + ' radiation or integration cases leading to a total of ' + str(n_problems) + ' problems to be solved')


    temp = int(switch_potential)*np.ones(n_problems, dtype='i')
    dset = utility.require_dataset(hdf5_data, structure.H5_NORMAL_VELOCITY_SWITCH_POTENTIAL, temp.shape, dtype='i')
    utility.set_hdf5_attributes(dset, structure.H5_NORMAL_VELOCITY_SWITCH_POTENTIAL_ATTR)
    dset[:] = temp

    temp = int(switch_free_surface)*np.ones(n_problems, dtype='i')
    dset = utility.require_dataset(hdf5_data, structure.H5_NORMAL_VELOCITY_SWITCH_FREE_SURFACE, temp.shape, dtype='i')
    utility.set_hdf5_attributes(dset, structure.H5_NORMAL_VELOCITY_SWITCH_FREE_SURFACE_ATTR)
    dset[:] = temp

    temp = int(switch_kochin)*np.ones(n_problems, dtype='i')
    dset = utility.require_dataset(hdf5_data, structure.H5_NORMAL_VELOCITY_SWITCH_KOCHIN, temp.shape, dtype='i')
    utility.set_hdf5_attributes(dset, structure.H5_NORMAL_VELOCITY_SWITCH_KOCHIN_ATTR)
    dset[:] = temp

    logger.info('Saved wether or not to compute the potential, the free surface visualization coefficients '
        'and the Kochin function for each of the ' + str(n_problems) + ' problems to be solved; '
        'respectively at ' + structure.H5_NORMAL_VELOCITY_SWITCH_POTENTIAL + ', '
        + structure.H5_NORMAL_VELOCITY_SWITCH_FREE_SURFACE + ', '
        + structure.H5_NORMAL_VELOCITY_SWITCH_KOCHIN)

    dset = utility.require_dataset(hdf5_data, structure.H5_NORMAL_VELOCITY_VELOCITIES, normal_velocity.shape, dtype='F')
    utility.set_hdf5_attributes(dset, structure.H5_NORMAL_VELOCITY_VELOCITIES_ATTR)
    dset[:, :] = normal_velocity
    logger.info('Saved the normal velocities at ' + structure.H5_NORMAL_VELOCITY_VELOCITIES
        + ' with characteristics ' + str(dset))


    #fk_force_f = fk_force.flatten()
    #fk_force_o = np.vstack((np.abs(fk_force_f), np.arctan2(np.imag(fk_force_f), np.real(fk_force_f)))).transpose()
    logger.info('Computing the magnitude and angle of the FK forces')
    fk_force_o = np.zeros((n_integration*n_w, 2*n_beta+2*n_radiation), dtype='f')
    idx = 0
    for k in range(n_integration):
        for i in range(n_w):
            for c in range(n_beta):
                fk_force_o[idx, 2*c] = np.abs(fk_force[i, c, k])
                fk_force_o[idx, 2*c+1] = np.arctan2(np.imag(fk_force[i, c, k]), np.real(fk_force[i, c, k]))

            for c in range(2*n_radiation):
                fk_force_o[idx, 2*n_beta + c] = 0
            idx += 1


    dset = utility.require_dataset(hdf5_data, structure.H5_RESULTS_FK_FORCES, fk_force_o.shape, dtype='f')
    utility.set_hdf5_attributes(dset, structure.H5_RESULTS_FK_FORCES_ATTR)
    dset[:, :] = fk_force_o
    logger.info('Saved the magnitude and angle of the fk forces at' 
        + str(structure.H5_RESULTS_FK_FORCES) + ' with characteristics: ' + str(dset))

    dset = utility.require_dataset(hdf5_data, structure.H5_RESULTS_FK_FORCES_RAW, fk_force.shape, dtype='F')
    utility.set_hdf5_attributes(dset, structure.H5_RESULTS_FK_FORCES_RAW_ATTR)
    dset[:, :, :] = fk_force
    logger.info('Saved the raw imaginary fk forces numbers at' 
                + str(structure.H5_RESULTS_FK_FORCES_RAW) + ' with characteristics: ' + str(dset))



    fk_force_tec_file = utility.get_setting(settings.FK_FORCE_TEC_FILE, custom_config, 'FK_FORCE_TEC_FILE')
    if fk_force_tec_file:
        logger.info('Converting the FK forces to the Tecplot format at the file ' + str(fk_force_tec_file))
        write_fk_force_tec(int_case, fk_force, w, beta, fk_force_tec_file)

    #free_surface_v = [[-0.5*l_x+l_x*i/(n_x-1), -0.5*l_y+l_y*j/(n_y-1), 0.] for i in range(n_x) for j in range(
    #    n_y)]

    
    free_surface_v = np.zeros((3, n_x*n_y))
    logger.info('Computing the free surface coefficients matrix of shape ' + str(free_surface_v.shape))
    k = 0
    for i in range(n_x):
        for j in range(n_y):
            free_surface_v[0, k] = -0.5*l_x+l_x*i/(n_x-1)
            free_surface_v[1, k] = -0.5*l_y+l_y*j/(n_y-1)
            free_surface_v[2, k] = 0.
            k += 1

    #free_surface_v = np.array(free_surface_v)
    dset = utility.require_dataset(hdf5_data, structure.H5_MESH_FREE_SURFACE_VECTORS, free_surface_v.shape, dtype='f')
    utility.set_hdf5_attributes(dset, structure.H5_MESH_FREE_SURFACE_VECTORS_ATTR)
    dset[:, :] = free_surface_v
    logger.info('Saved the free surface coefficients at ' + str( structure.H5_MESH_FREE_SURFACE_VECTORS) 
                + ' with characteristics ' + str(dset))


    free_surface_v = np.zeros((0, 0))
    logger.info('Computing the free surface index matrix (identifying the corresponding problem) of shape ' + str(free_surface_v.shape))
    if (n_x-1) > 0 and (n_y-1) >0:
        #free_surface_v = [[j+i*n_y, j+1+i*n_y, j+1+(i+1)*n_y, j+(i+1)*n_y] for i in range(n_x-1) for j in
                        #range(n_y-1)]
        free_surface_v = np.zeros((4, (n_x-1)*(n_y-1)))
        k = 0
        for i in range(n_x-1):
            for j in range(n_y-1):
                free_surface_v[0, k] = j+i*n_y
                free_surface_v[1, k] = j+1+i*n_y
                free_surface_v[2, k] = j+1+(i+1)*n_y
                free_surface_v[3, k] = j+(i+1)*n_y
                k += 1
    #free_surface_v = np.array(free_surface_v)
    dset = utility.require_dataset(hdf5_data, structure.H5_MESH_FREE_SURFACE_INDEX, free_surface_v.shape, dtype='f')
    utility.set_hdf5_attributes(dset, structure.H5_MESH_FREE_SURFACE_INDEX_ATTR)
    dset[:, :] = free_surface_v
    logger.info('Saved the free surface index at ' + str(structure.H5_MESH_FREE_SURFACE_INDEX) 
        + ' with characteristics ' + str(dset))


    # Generate Kochin
    kochin = np.array([])
    if n_theta > 0:
        logger.info('Computing the ' + str(n_theta) + ' angles for which the Kochin function will be calculated equi-distantly '
        ' from ' + str(theta_min) + ' to ' + str(theta_max))
        if n_theta > 1:
            kochin = [(theta_min+(theta_max-theta_min)*j/(n_theta-1))*np.pi/180. for j in range(n_theta)]
        else:
            kochin = [theta_min*np.pi/180.]


    kochin = np.array(kochin)
    dset = utility.require_dataset(hdf5_data, structure.H5_MESH_KOCHIN, kochin.shape, dtype='f', maxshape=(None, ))
    utility.set_hdf5_attributes(dset, structure.H5_MESH_KOCHIN_ATTR)
    dset[:] = kochin
    logger.info('Saved the angle for which the Kochin function will be calculated at '
     + str(structure.H5_MESH_KOCHIN) + ' with characteristics ' + str(dset))

    # Save index of cases

    out = np.array([[k+1, int_case[k].body+1, int_case[k].mode+1] for k in range(n_integration)])
    dset = utility.require_dataset(hdf5_data, structure.H5_RESULTS_CASE_FORCE, out.shape, dtype='i')
    utility.set_hdf5_attributes(dset, structure.H5_RESULTS_CASE_FORCE_ATTR)
    dset[:, :] = out
    logger.info('Saved correlation between force ID number and body ID number at '
        + structure.H5_RESULTS_CASE_FORCE)

    out = np.array([[k+1, rad_case[k].body+1, rad_case[k].mode+1] for k in range(n_radiation)])
    dset = utility.require_dataset(hdf5_data, structure.H5_RESULTS_CASE_MOTION, out.shape, dtype='i')
    utility.set_hdf5_attributes(dset, structure.H5_RESULTS_CASE_MOTION_ATTR)
    dset[:, :] = out
    logger.info('Saved correlation between ID number of radiation problem,'
        ' ID of body which is actually moving and the number of degree of freedom '
        + structure.H5_RESULTS_CASE_MOTION)


    dset = utility.require_dataset(hdf5_data, structure.H5_RESULTS_CASE_BETA, beta.shape, dtype='f', maxshape=(None))
    utility.set_hdf5_attributes(dset, structure.H5_RESULTS_CASE_BETA_ATTR)
    dset[:] = beta
    logger.info('Saved the wave directions at ' + structure.H5_RESULTS_CASE_BETA)

    dset = utility.require_dataset(hdf5_data, structure.H5_RESULTS_CASE_W, w.shape, dtype='f', maxshape=(None))
    utility.set_hdf5_attributes(dset, structure.H5_RESULTS_CASE_W_ATTR)
    dset[:] = w
    logger.info('Saved the wave frequencies at ' + structure.H5_RESULTS_CASE_W)

    out = np.array([(theta_min+(theta_max-theta_min)*k/(n_theta-1))*np.pi/180. for k in range(n_theta)])
    dset = utility.require_dataset(hdf5_data, structure.H5_RESULTS_CASE_THETA, out.shape, dtype='f', maxshape=(None))
    utility.set_hdf5_attributes(dset, structure.H5_RESULTS_CASE_THETA_ATTR)
    dset[:] = out
    logger.info('Saved the kochin angles at ' + structure.H5_RESULTS_CASE_THETA)

    # Save radiation cases

    out = np.array([[rad_case[k].body+1, rad_case[k].i_case+1,  rad_case[k].direction[0], rad_case[k].direction[1], rad_case[k].direction[2],  rad_case[k].axis[0],  rad_case[k].axis[1] ,  rad_case[k].axis[2]] for k in range(n_radiation)])
    dset = utility.require_dataset(hdf5_data, structure.H5_RESULTS_CASE_RADIATION, out.shape, dtype='f')
    utility.set_hdf5_attributes(dset, structure.H5_RESULTS_CASE_RADIATION_ATTR)
    dset[:, :] = out
    logger.info('Saved the radiation cases (directions and axis) at ' + structure.H5_RESULTS_CASE_RADIATION)

    dset = utility.require_dataset(hdf5_data, structure.H5_RESULTS_CASE_BETA, beta.shape, dtype='f')
    utility.set_hdf5_attributes(dset, structure.H5_RESULTS_CASE_BETA_ATTR)
    dset[:] = beta
    logger.info('Saved the radiation wave directions at ' + structure.H5_RESULTS_CASE_BETA)


    # All the following are switches already logged
    switch_ode_influence = utility.get_setting(settings.USE_ODE_INFLUENCE_COEFFICIENTS, custom_config,
                                       'USE_ODE_INFLUENCE_COEFFICIENTS')

    use_higher_order = utility.get_setting(settings.USE_HIGHER_ORDER, custom_config,
                                       'USE_HIGHER_ORDER')

    num_panel_higher_order = utility.get_setting(settings.NUM_PANEL_HIGHER_ORDER, custom_config,
                                       'NUM_PANEL_HIGHER_ORDER')

    b_spline_order = utility.get_setting(settings.B_SPLINE_ORDER, custom_config,
                                       'B_SPLINE_ORDER')

    use_dipoles_implementation = utility.get_setting(settings.USE_DIPOLES_IMPLEMENTATION, custom_config,
                                       'USE_DIPOLES_IMPLEMENTATION')

    compute_yaw_moment = utility.get_setting(settings.COMPUTE_YAW_MOMENT, custom_config,
                                       'COMPUTE_YAW_MOMENT')

    compute_drift_forces = utility.get_setting(settings.COMPUTE_DRIFT_FORCES, custom_config,
                                       'COMPUTE_DRIFT_FORCES')

    thin_panels = utility.get_setting(settings.THIN_PANELS, custom_config,
                                       'THIN_PANELS')

    if num_panel_higher_order is not None and num_panel_higher_order > 0:
        dset = utility.require_dataset(hdf5_data, structure.H5_SOLVER_NUM_PANEL_HIGHER_ORDER, (1, ), dtype='i')
        utility.set_hdf5_attributes(dset, structure.H5_SOLVER_NUM_PANEL_HIGHER_ORDER_ATTR)
        dset[:] = int(num_panel_higher_order)
        logger.info('Saved the number of panel per patch in the higher order method at ' + structure.H5_SOLVER_NUM_PANEL_HIGHER_ORDER)

    if b_spline_order is not None and b_spline_order > 0:
        dset = utility.require_dataset(hdf5_data, structure.H5_SOLVER_B_SPLINE_ORDER, (1, ), dtype='i')
        utility.set_hdf5_attributes(dset, structure.H5_SOLVER_B_SPLINE_ORDER_ATTR)
        dset[:] = int(b_spline_order)
        logger.info('Saved The order of the B-Spline for the potential in the higher order panel method at' 
            + structure.H5_SOLVER_B_SPLINE_ORDER)

    if use_higher_order is not None:
        dset = utility.require_dataset(hdf5_data, structure.H5_SOLVER_USE_HIGHER_ORDER, (1, ), dtype='i')
        utility.set_hdf5_attributes(dset, structure.H5_SOLVER_USE_HIGHER_ORDER_ATTR)
        dset[:] = int(use_higher_order)
        logger.info('Saved the switch for using of the higher order panel method at ' +  structure.H5_SOLVER_USE_HIGHER_ORDER)


    if switch_ode_influence is not None:
        temp = int(switch_ode_influence)*np.ones(n_problems, dtype='i')
        dset = utility.require_dataset(hdf5_data, structure.H5_SOLVER_SWITCH_ODE_INFLUENCE, temp.shape, dtype='i')
        utility.set_hdf5_attributes(dset, structure.H5_SOLVER_SWITCH_ODE_INFLUENCE_ATTR)
        dset[:] = temp
        logger.info('Saved the switch for computing the influence coefficients using an Ordinary Differential equation at '
         +  structure.H5_SOLVER_SWITCH_ODE_INFLUENCE)

    if use_dipoles_implementation is not None:
        dset = utility.require_dataset(hdf5_data, structure.H5_SOLVER_USE_DIPOLES_IMPLEMENTATION, (1, ), dtype='i')
        utility.set_hdf5_attributes(dset, structure.H5_SOLVER_USE_DIPOLES_IMPLEMENTATION_ATTR)
        dset[:] = int(use_dipoles_implementation)
        logger.info('Saved the switch for using the dipoles Implementation for thin bodies at '
            + structure.H5_SOLVER_USE_DIPOLES_IMPLEMENTATION)

    if compute_yaw_moment is not None:
        dset = utility.require_dataset(hdf5_data, structure.H5_SOLVER_COMPUTE_YAW_MOMENT, (1, ), dtype='i')
        utility.set_hdf5_attributes(dset, structure.H5_SOLVER_COMPUTE_YAW_MOMENT_ATTR)
        dset[:] = int(compute_yaw_moment)
        logger.info('Saved the switch for computing the yaw moment at ' + structure.H5_SOLVER_COMPUTE_YAW_MOMENT)

    if compute_drift_forces is not None:
        dset = utility.require_dataset(hdf5_data, structure.H5_SOLVER_COMPUTE_DRIFT_FORCES, (1, ), dtype='i')
        utility.set_hdf5_attributes(dset, structure.H5_SOLVER_COMPUTE_DRIFT_FORCES_ATTR)
        dset[:] = int(compute_drift_forces)
        logger.info('Saved the switch for computing the drift forces at ' + structure.H5_SOLVER_COMPUTE_DRIFT_FORCES)

    if thin_panels is not None:
        temp = np.zeros(mesh.n_panels, dtype='i')
        for idx in thin_panels:
            if idx == -1:
                temp = np.ones(mesh.n_panels, dtype='i')
                break
            elif idx >= 0:
                temp[idx] = 1
        dset = utility.require_dataset(hdf5_data, structure.H5_SOLVER_THIN_PANELS, temp.shape, dtype='i')
        utility.set_hdf5_attributes(dset, structure.H5_SOLVER_THIN_PANELS_ATTR)
        dset[:] = temp
        logger.info('Saved the thin or not panels index at ' + structure.H5_SOLVER_THIN_PANELS)

    utility.log_exit(logger, signature, [None])


def preprocess(custom_config):
    """
    Configure and then run the preprocessor

    Args:
        custom_config, dict The custom configuration dictionary
    """
    signature = __name__ + '.preprocess(custom_config)'
    logger = logging.getLogger(__name__)
    utility.log_entrance(logging.getLogger(__name__), signature,
                        {'custom_config': custom_config})

    if not custom_config:
        custom_config = {}

    # Check if custom_config is a valid dict
    utility.check_type_value(custom_config, 'The input custom_config', dict, True)

    hdf5_file = utility.get_setting(settings.HDF5_FILE, custom_config, 'HDF5_FILE')
    nemoh_cal = utility.get_setting(settings.NEMOH_CALCULATIONS_FILE, custom_config, 'NEMOH_CALCULATIONS_FILE')
    input_file = utility.get_setting(settings.NEMOH_INPUT_FILE, custom_config, 'NEMOH_INPUT_FILE')


    # Check if hdf5_file is a string
    utility.check_str(hdf5_file, 'The path to the hdf5 file configured by HDF5_FILE')
    # Check wether or not nemoh_cal and input_file are valid string
    nemoh_cal_validation = utility.validate_str(nemoh_cal)
    input_file_validation = utility.validate_str(input_file)

    # If both input_file and nemoh_cal are not valid string then hdf5_file must be a valid file
    if not nemoh_cal_validation and not input_file_validation:
        utility.check_is_file(hdf5_file, 'The path to the hdf5 file configured by HDF5_FILE')
    else:
        # Otherwise we make sure that the directory to the hdf5_file exists and we will create the hdf5_file as part of the program
        utility.mkdir_p(os.path.abspath(os.path.dirname(hdf5_file)))
        utility.touch(hdf5_file) # This is to solve a bug in old version of hdf5: If file does not exist
        # when trying to append, it will create error

 

    with h5py.File(hdf5_file, "a") as hdf5_db:
        if nemoh_cal_validation:
            # If nemoh_cal is a valid string then it must be a valid file and we will convert its entries into the hdf5_file
            utility.check_is_file(nemoh_cal, 'The path to the nemoh calculations configured by NEMOH_CALCULATIONS_FILE')
            utility.convert_calculations(nemoh_cal, hdf5_db)

        if input_file_validation:
            # If input_file is a valid string then it must be a valid file and we will convert its entries into the hdf5_file
            utility.check_is_file(input_file, 'The path to the nemoh input configured by NEMOH_INPUT_FILE')
            utility.convert_input(input_file, hdf5_db)

        remove_irregular_frequencies = utility.get_setting(settings.REMOVE_IRREGULAR_FREQUENCIES, custom_config,
                                       'REMOVE_IRREGULAR_FREQUENCIES')
        if remove_irregular_frequencies is not None:
            dset = utility.require_dataset(hdf5_db, structure.H5_SOLVER_REMOVE_IRREGULAR_FREQUENCIES, (1, ), dtype='i')
            utility.set_hdf5_attributes(dset, structure.H5_SOLVER_REMOVE_IRREGULAR_FREQUENCIES_ATTR)
            dset[:] = int(remove_irregular_frequencies)
        else:
            settings.REMOVE_IRREGULAR_FREQUENCIES = hdf5_db.get(structure.H5_SOLVER_REMOVE_IRREGULAR_FREQUENCIES)[0]

        run(hdf5_db, custom_config)
        utility.log_and_print(logger, 'The preprocessing results are saved in the hdf5 file '
                              + utility.get_abs(hdf5_file))

    utility.log_exit(logging.getLogger(__name__), signature, [None])


def run_as_process(custom_config, queue):
    utility.setup_subprocess_logging(queue, logging.getLogger())
    return preprocess(custom_config)


if __name__ == '__main__':
    utility.setup_logging(default_conf_path=settings.LOGGING_CONFIGURATION_FILE, logging_path=settings.LOG_FILE)

    try:
        preprocess({})
        print('Pre processing successfully completed \n')
    except Exception as e:
        # exc_info=True means the stack trace will be printed automatically
        print('There was an error when running the application. Check the log file for more details \n')
        logging.getLogger(__name__).error('Program halted due to a fatal error whose detail is as follow: ',
                                          exc_info=True)