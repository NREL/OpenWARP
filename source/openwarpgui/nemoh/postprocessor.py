#!/usr/bin/env python
"""
This is the main program for the post processor.
It may be used to calculate RAOs and plot the free surface, wave elevation.
The aim of the postProcessor is to post process the results in order to provide the relevant quantities
(added mass, radiation damping, excitation force) in the usual format. It also provides
a framework to make relevant calculations.


Once run successfully, the following results files are created:

settings.RADIATION_COEFFICIENTS_TEC_FILE:
            This file contains the added mass and damping forces for the radiation problems.

settings.DIFFRACTION_FORCE_TEC_FILE: This file contains the diffraction force for the diffraction problems.

settings.EXCITATION_FORCE_TEC_FILE: This file contains the excitation force for the diffraction problems.

Changes in version 1.1:
    Added possibility to run the code with custom settings

Changes in version 1.2 (Implementation of Higher Order Panel Methods):
    Disable computation of wave elevation if higher panel method is used

Changes in version 1.3 (Implementation of Higher Order Panel Methods):
    Disable computation of wave elevation if dipoles implementation is used

Updated since version 1.4: (Hydrodynamic Data Exporter Assembly v1.0)
    Stored excitation forces, radiation damping coefficients, added mass
    in hdf5 file.

Changes in version 1.5 (Irregular Frequencies Assembly)
       Disable computation of wave elevation if  Irregular frequencies are
       removed.
       Applied some bug fixes to allow the shape of hdf5 file dataset 
       to be automatically resized.

Changes in version 1.6 (OpenWarp - Add Logging Functionality)
       Added support for logging.

Changes in version 1.7 (OPENWARP - FIX WAVE FREQUENCY AND DIRECTION CRASH BUG):
    1. Changed the way we do logging from this module when it is run
    as a child process.

"""

import utility
import numpy as np
import sys
import h5py
import structure
import settings
import os
from utility import cih
from utility import sih
import logging
import preprocessor


from models import TResult
from models import TIRF

__author__ = "yedtoss"
__copyright__ = "Copyright (C) 2014-2016 TopCoder Inc. All rights reserved."
__version__ = "1.7"


def get_irf(hdf5_data, result):
    """
    Gets the irf from the hdf5 file
    Args:
        hdf5_data: object, the hdf5 opened file
        result: object the hydrodynamic coefficients cases
    Returns:
        the irf
    """

    signature = __name__ + '.get_irf(hdf5_data, result)'
    logger = logging.getLogger(__name__)
    utility.log_entrance(logger, signature,
                         {"hdf5_data": str(hdf5_data), 'result': str(result)})

    dset = hdf5_data.get(structure.H5_COMPUTE_IRF)
    utility.check_dataset_type(dset,
                               name=str(structure.H5_COMPUTE_IRF_ATTR['description']),
                               location=structure.H5_COMPUTE_IRF)
    switch = dset[0]

    dset = hdf5_data.get(structure.H5_IRF_TIME_STEP)
    utility.check_dataset_type(dset,
                               name=str(structure.H5_IRF_TIME_STEP_ATTR['description']),
                               location=structure.H5_IRF_TIME_STEP)
    time_step = dset[0]


    dset = hdf5_data.get(structure.H5_IRF_DURATION)
    utility.check_dataset_type(dset,
                               name=str(structure.H5_IRF_DURATION_ATTR['description']),
                               location=structure.H5_IRF_DURATION)
    duration = dset[0]

    irf = TIRF()

    if switch == 1:
        irf = TIRF(int(duration/time_step), result.n_radiation, result.n_integration)
        for i in range(irf.n_time):
            irf.time[i] = i*time_step
    irf.switch = switch

    utility.log_exit(logger, signature, [str(irf)])

    return irf


def compute_irf(result, irf):
    """
    Computes the froude krylov forces for the given irf
    Args:
        result: object the hydrodynamic coefficients cases
        irf: object, the input irf
    Returns:
        the irf with the froude krylov forces computed.
    """
    signature = __name__ + '.compute_irf(result, irf)'
    logger = logging.getLogger(__name__)
    utility.log_entrance(logger, signature,
                         {"irf": str(irf), 'result': str(result)})

    logger.info('Computing the froude krylov forces for the given irf')
    for i in range(irf.n_time):
        for j in range(result.n_radiation):
            for k in range(result.n_integration):
                irf.k[i,j,k] = 0.
                for l in range(result.n_w -1):

                    irf.k[i, j, k] += - 0.5* (result.w[l+1]-result.w[l]) * (result.radiation_damping[l,
                                                                                                                   j,
                     k]*np.cos(result.w[l]*irf.time[i]) + result.radiation_damping[l+1, j, k]*np.cos(result.w[l+1]*irf.time[i]))

                irf.k[i, j, k] = (irf.k[i, j, k] * 2)/np.pi



    logger.info('Computing the added mass for the given irf')
    cm = np.zeros(result.n_w)
    for j in range(result.n_radiation):
        for k in range(result.n_integration):
            irf.added_mass[j, k] = 0
            for l in range(result.n_w):
                cm[l] = 0
                for i in range(irf.n_time-1):
                    cm[l] += 0.5*(irf.time[i+1]-irf.time[i])*(irf.k[i,j,k]*np.sin(result.w[l]*irf.time[i]) +
                                                                     irf.k[i+1, j, k]*np.sin(result.w[l]*irf.time[i+1]))
                cm[l]=(result.added_mass[l,j,k] + cm[l])/result.w[l]
                irf.added_mass[j,k]=irf.added_mass[j,k] +cm[l]
            irf.added_mass[j,k] = irf.added_mass[j,k]/result.n_w

    utility.log_exit(logger, signature, [str(irf)])

    return irf


def save_irf(irf, filename):
    """
    Saves the irf to a file in the tec format
    Args:
        irf: object, the irf
        filename: string, The path to the file where to save the irf
    """
    signature = __name__ + '.save_irf(irf, filename)'
    logger = logging.getLogger(__name__)
    utility.log_entrance(logger, signature,
                         {"irf": str(irf),
                          'filename': str(filename)})
    utility.mkdir_p(os.path.abspath(os.path.dirname(filename)))
    with open(filename, 'w') as inp:
        inp.write('VARIABLES="Time (s)"\n')
        for k in range(irf.n_integration):
            inp.write('"AddedMass '+ str(k+1) + '" "IRF ' + str(k+1) +'"\n')
        for j in range(irf.n_radiation):
            inp.write('Zone t="DoF ' + str(j+1) + '",I=\t' + str(irf.n_time) + ',F=POINT\n')
            for i in range(irf.n_time):
                inp.write(' ' + str(irf.time[i]) + ' ')
                s = ''
                for k in range(irf.n_integration):
                    s += ' ' + str(irf.added_mass[j, k]) + ' ' + ' ' + str(irf.k[i, j, k]) + ' '
                inp.write(s + '\n')

    utility.log_and_print(logger, utility.get_abs(filename) + ' contains the irf in tec format.')
    utility.log_exit(logger, signature, [None])


def read_results(hdf5_data):
    """
    Read the hydrodynamic coefficients cases from the hdf5 file
    Args:
        hdf5_data: object, the hdf5 opened file
    Returns:
        the hydrodynamic coefficients cases
    """
    signature = __name__ + '.read_results(hdf5_data)'
    logger = logging.getLogger(__name__)
    utility.log_entrance(logger, signature,
                         {"hdf5_data": str(hdf5_data)})

    idx_force = hdf5_data.get(structure.H5_RESULTS_CASE_FORCE)
    utility.check_dataset_type(idx_force,
                               name=str(structure.H5_RESULTS_CASE_FORCE_ATTR['description']),
                               location=structure.H5_RESULTS_CASE_FORCE)
    n_integration = idx_force.shape[0]

    idx_radiation = hdf5_data.get(structure.H5_RESULTS_CASE_MOTION)
    utility.check_dataset_type(idx_radiation,
                               name=str(structure.H5_RESULTS_CASE_MOTION_ATTR['description']),
                               location=structure.H5_RESULTS_CASE_MOTION)
    n_radiation = idx_radiation.shape[0]

    beta = hdf5_data.get(structure.H5_RESULTS_CASE_BETA)
    utility.check_dataset_type(beta,
                               name=str(structure.H5_RESULTS_CASE_BETA_ATTR['description']),
                               location=structure.H5_RESULTS_CASE_BETA)
    n_beta = beta.shape[0]

    w = hdf5_data.get(structure.H5_RESULTS_CASE_W)
    utility.check_dataset_type(beta,
                               name=str(structure.H5_RESULTS_CASE_W_ATTR['description']),
                               location=structure.H5_RESULTS_CASE_W)
    n_w = w.shape[0]

    theta = hdf5_data.get(structure.H5_RESULTS_CASE_THETA)
    n_theta = theta.shape[0]

    result = TResult(n_w, n_radiation, n_integration, n_theta, n_beta)
    result.idx_force = idx_force
    result.idx_radiation = idx_radiation
    result.beta = beta
    result.w = w
    result.theta = theta

    forces = hdf5_data.get(structure.H5_RESULTS_FORCES)
    for k in range(n_integration):
        c = 0
        for i in range(n_w):
            for j in range(n_beta):
                result.diffraction_force[i, j, k] = forces[k, c]*np.exp(complex(0, 1)*forces[k, c+1])
                c += 2
            for j in range(n_radiation):
                result.added_mass[i, j, k] = forces[k, c]
                result.radiation_damping[i, j, k] = forces[k, c+1]
                c += 2

    result.froudkrylov_force = hdf5_data.get(structure.H5_RESULTS_FK_FORCES_RAW)

    utility.log_exit(logger, signature, [str(result)])
    return result


def save_radiation_coefficients(result, filename):
    """
    Saves the radiation coefficient to a file in tec format
    Args:
        result: object, the hydrodynamic coefficients cases
        filename: The path to the file where to save
    """
    signature = __name__ + '.save_radiation_coefficients(result, filename)'
    logger = logging.getLogger(__name__)
    utility.log_entrance(logger, signature,
                         {"result": str(result),
                          'filename': str(filename)})

    utility.mkdir_p(os.path.abspath(os.path.dirname(filename)))
    with open(filename, 'w') as inp:

        inp.write('VARIABLES="w (rad/s)"\n')

        for k in range(result.n_integration):
            s = '"A\t' + str(result.idx_force[k, 1]) + '\t' + str(result.idx_force[k, 2])
            s = s + '" "B\t' + str(result.idx_force[k, 1]) + '\t' + str(result.idx_force[k, 2]) + '"\n'
            inp.write(s)

        for j in range(result.n_radiation):
            s = 'Zone t="Motion of body\t' + str(result.idx_radiation[j, 1])
            s += '\tin DoF\t' + str(result.idx_radiation[j, 2]) + '",I=\t' + str(result.n_w) + ',F=POINT\n'
            inp.write(s)
            for i in range(result.n_w):
                s = str(result.w[i]) + '\t'
                for k in range(result.n_integration):
                    s += str(result.added_mass[i, j, k]) + '\t' + str(result.radiation_damping[i, j, k]) + '\t'
                inp.write(s + '\n')

    utility.log_and_print(logger, utility.get_abs(filename) + ' contains the radiation coefficients in tec format.')
    utility.log_exit(logger, signature, [None])


def save_diffraction_force(result, filename):
    """
    Saves the diffraction forces to a file in tec format
    Args:
        result: object, the hydrodynamic coefficients cases
        filename: The path to the file where to save
    """
    signature = __name__ + '.save_diffraction_force(result, filename)'
    logger = logging.getLogger(__name__)
    utility.log_entrance(logger, signature,
                         {"result": str(result),
                          'filename': str(filename)})
    utility.mkdir_p(os.path.abspath(os.path.dirname(filename)))
    with open(filename, 'w') as inp:
        inp.write('VARIABLES="w (rad/s)"\n')
        for k in range(result.n_integration):
            s = '"abs(F\t' + str(result.idx_force[k, 1]) + '\t' + str(result.idx_force[k, 2])
            s += ')" "angle(F\t' + str(result.idx_force[k, 1]) + '\t' + str(result.idx_force[k, 2]) + ')"\n'
            inp.write(s)

        for j in range(result.n_beta):
            s = 'Zone t="Diffraction force - beta = ' + str(result.beta[j]*180./(4.* np.arctan(1.0)))
            s += ' deg",I=\t' + str(result.n_w) + ',F=POINT\n'
            inp.write(s)
            for i in range(result.n_w):
                s = str(result.w[i]) + '\t'
                for k in range(result.n_integration):
                    s += str(np.abs(result.diffraction_force[i, j, k])) + '\t'
                    s += str(np.arctan2(np.imag(result.diffraction_force[i, j, k]), np.real(result.diffraction_force[
                        i, j, k]))) + '\t'
                inp.write(s + '\n')

    utility.log_and_print(logger, utility.get_abs(filename) + ' contains the diffraction forces in tec format.')
    utility.log_exit(logger, signature, [None])


def save_excitation_force(result, filename):
    """
    Saves the excitation forces to a file in tec format
    Args:
        result: object, the hydrodynamic coefficients cases
        filename: The path to the file where to save
    """
    signature = __name__ + '.save_excitation_force(result, filename)'
    logger = logging.getLogger(__name__)
    utility.log_entrance(logger, signature,
                         {"result": str(result),
                          'filename': str(filename)})
    utility.mkdir_p(os.path.abspath(os.path.dirname(filename)))
    with open(filename, 'w') as inp:
        inp.write('VARIABLES="w (rad/s)"\n')

        for k in range(result.n_integration):
            s = '"abs(F\t' + str(result.idx_force[k, 1]) + '\t' + str(result.idx_force[k, 2])
            s += ')" "angle(F\t' + str(result.idx_force[k, 1]) + '\t' + str(result.idx_force[k, 2]) + ')"\n'
            inp.write(s)

        for j in range(result.n_beta):
            s = 'Zone t="Diffraction force - beta = ' + str(result.beta[j]*180./(4.*np.arctan(1.0)))
            s += ' deg",I=' + str(result.n_w) + ',F=POINT\n'
            inp.write(s)
            for i in range(result.n_w):
                s = str(result.w[i]) + '\t'
                for k in range(result.n_integration):
                    temp = result.diffraction_force[i, j, k] + result.froudkrylov_force[i, j, k]
                    s += str(np.abs(temp)) + '\t'
                    s += str(np.arctan2(np.imag(temp), np.real(temp))) + '\t'
                inp.write(s + '\n')

    utility.log_and_print(logger, utility.get_abs(filename) + ' contains the excitation forces in tec format.')
    utility.log_exit(logger, signature, [None])


def compute_raos(raos, result):
    """
    Computes the raos
    Args:
        raos: object, the input raos
        result: object, the hydrodynamic coefficients cases
    Returnss:
        the computed raos
    """
    return raos


def compute_wave_elevation(hdf5_data, environment, iw, ibeta, raos, result):
    """
    Computes the wave elevation
    Args:
        hdf5_data: object, the hdf5 opened file
        environment: object, the environment
        iw: int, the index of the wave frequency to use
        ibeta: int, the index of the wave direction to use
        raos: object, the raos
        result: the hydrodynamic coefficients cases

    Returns:
        A dictionary containing the wave elevation variables
    """

    signature = __name__ + '.compute_wave_elevation(hdf5_data, environment, iw, ibeta, raos, result)'
    logger = logging.getLogger(__name__)
    utility.log_entrance(logger, signature,
                         {"hdf5_data": str(hdf5_data),
                          "environment": str(environment),
                          "iw": str(iw),
                          "ibeta": str(ibeta),
                          "raos": str(raos),
                          "result": str(result)})

    dset = hdf5_data.get(structure.H5_FREE_SURFACE_POINTS_X)
    utility.check_dataset_type(dset, name=str(structure.H5_FREE_SURFACE_POINTS_X_ATTR['description']),
                               location=structure.H5_FREE_SURFACE_POINTS_X)
    nx = dset[0]

    dset = hdf5_data.get(structure.H5_FREE_SURFACE_POINTS_Y)
    utility.check_dataset_type(dset, name=str(structure.H5_FREE_SURFACE_POINTS_Y_ATTR['description']),
                               location=structure.H5_FREE_SURFACE_POINTS_Y)
    ny = dset[0]

    dset = hdf5_data.get(structure.H5_FREE_SURFACE_DIMENSION_X)
    utility.check_dataset_type(dset, name=str(structure.H5_FREE_SURFACE_DIMENSION_X_ATTR['description']),
                               location=structure.H5_FREE_SURFACE_DIMENSION_X)
    lx = dset[0]

    dset = hdf5_data.get(structure.H5_FREE_SURFACE_DIMENSION_Y)
    utility.check_dataset_type(dset, name=str(structure.H5_FREE_SURFACE_DIMENSION_Y_ATTR['description']),
     location=structure.H5_FREE_SURFACE_DIMENSION_Y)
    ly = dset[0]

    x = np.zeros(nx, dtype='f')
    y = np.zeros(ny, dtype='f')
    etai = np.zeros((nx, ny), dtype='F')
    etap = np.zeros((nx, ny), dtype='F')
    eta = np.zeros((nx, ny), dtype='F')

    for i in range(nx):
        x[i] = -0.5*lx+lx*(i)/(nx-1)

    for i in range(ny):
        y[i] = -0.5*ly+ly*(i)/(ny-1)

    w = result.w[iw]
    logger.info('Computing the wave number ...')
    kwave = utility.compute_wave_number(w, environment)
    logger.info('Wave number computed is ' + str(kwave))

    for i in range(nx):
        for j in range(ny):
            r = np.sqrt((x[i] - environment.x_eff)**2 + (y[j] - environment.y_eff)**2)
            theta = np.arctan2(y[j]-environment.y_eff, x[i]-environment.x_eff)
            k = 0
            while (k < result.n_theta -1) and (result.theta[k+1] < theta):
                k += 1
            if k == result.n_theta:
                raise ValueError(' Error: range of theta in Kochin coefficients is too small')

            coord = np.array([x[i], y[i], 0])
            one_wave = preprocessor.compute_one_wave(kwave, w, result.beta[ibeta], coord, environment)
            potential = one_wave["phi"]
            etai[i, j] = 1./environment.g*utility.II*w*one_wave["phi"]
            HKleft=0.
            HKright=0.
            for l in range(result.n_radiation):
                HKleft=HKleft+raos[l,iw,ibeta]*result.hkochin_radiation[iw,l,k]
                HKright=HKright+raos[l,iw,ibeta]*result.hkochin_radiation[iw,l,k+1]


            HKleft=HKleft+result.hkochin_diffraction[iw,ibeta,k]
            HKright=HKright+result.hkochin_diffraction[iw,ibeta,k+1]
            HKochin=HKleft+(HKright-HKleft)*(theta-result.theta[k])/(result.theta[k+1]-result.theta[k])
            if r > 0:
                potential = np.sqrt(kwave/(2.*np.pi*r))*cih(kwave,0.,environment.depth)* np.exp(utility.II*(kwave*r-0.25*np.pi))*HKochin
            else:
                potential = 0
            etap[i,j]=-utility.II*1./environment.g*utility.II*w*potential
            eta[i,j]=etai[i,j]+etap[i,j]

    rep = {"w": w, "x": x, "y": y, "eta": eta, "etai": etai, "etap": etap}
    utility.log_exit(logger, signature, [rep])
    return rep


def save_wave_elevation(w, etai, etap, eta, x, y, filename):
    """
    Save the wave elevation to a file in tec format
    Args:
        w: float, the wave frequency
        etai: 2D array, a wave elevation variable
        eta: 2D array, a wave elevation variable
        etap: 2D array, a wave elevation variable
        x: 1D array, a wave elevation variable
        y: 1D array, a wave elevation variable
        filename: string, the path to the file where to save
    """
    signature = __name__ + '.save_wave_elevation(w, etai, etap, eta, x, y, filename)'
    logger = logging.getLogger(__name__)
    utility.log_entrance(logger, signature,
                         {"w": w,
                          "etai": etai,
                          "etap": etap,
                          "eta": eta,
                          "x": x,
                          "y": y,
                          "filename": filename})

    nx = len(x)
    ny = len(y)
    utility.mkdir_p(os.path.abspath(os.path.dirname(filename)))
    with open(filename, 'w') as inp:
        s = 'VARIABLES="X" "Y" "etaI_C" "etaI_S" "etaP_C" "etaC_S" '
        s += '"etaI_C+etaP_C" "etaI_S+etaI_P" "|etaP|" "|etaI+etaP|"\n'
        inp.write(s)

        s = 'ZONE t="Wave frequency - w = '
        s += str(w)+ '",N=\t'+ str(nx*ny) + ', E=\t' + str((nx-1)*(nx-1)) + '\t , F=FEPOINT,ET=QUADRILATERAL\n'
        inp.write(s)
        for i in range(nx):
            for j in range(ny):
                s = str(x[i]) + '\t' + str(y[j]) + '\t' + str(np.real(etai[i, j])) + '\t'
                s += str(np.imag(etai[i,j])) + '\t' + str(np.real(etap[i, j])) + '\t'
                s += str(np.imag(etap[i,j])) + '\t' + str(np.real(etai[i,j]+etap[i,j])) + '\t'
                s += str(np.imag(etai[i,j]+etap[i,j])) + '\t' + str(np.abs(etap[i,j])) + '\t'
                s += str(np.abs(eta[i,j])) + '\n'
                inp.write(s)

        for i in range(nx-1):
            for j in range(ny-1):
                s = str(j+1+i*ny) + '\t' + str(j+1+(i+1)*ny) + '\t' + str(j+2+ i*ny) + '\n'
                inp.write(s)

    utility.log_and_print(logger, utility.get_abs(filename) + ' contains the wave elevation in tec format.')
    utility.log_exit(logger, signature, [None])


def run(hdf5_data, custom_config):
    """
    This function run the postprocessor
    Args:
        hdf5_data: object, the hdf5 opened file
        custom_config, dict The custom configuration dictionary
    """

    logger = logging.getLogger(__name__)
    signature = __name__ + '.run(hdf5_data, custom_config)'
    # No need to log the parameter of the method here as it will only be duplicate.
    # This function is never called directly by the user and always call from the postprocess function
    # which already logs the configuration.
    utility.log_entrance(logger, signature,
                         {})

    logger.info('Initialisation the post processing steps')

    logger.info('Reading environment data ...')
    environment = utility.read_environment(hdf5_data)
    logger.info('Read environment data' + str(environment))

    logger.info('Reading simulation results')
    result = read_results(hdf5_data)
    logger.info('Read solver result ' + str(result))

    logger.info('Post processing initialisation done !')

    # Saving to hdf5 file
    dset = utility.require_dataset(hdf5_data, structure.H5_RESULTS_ADDED_MASS, result.added_mass.shape, dtype='f')
    utility.set_hdf5_attributes(dset, structure.H5_RESULTS_ADDED_MASS_ATTR)
    dset[:, :, :] = result.added_mass
    logger.info('Saved ' + str(structure.H5_RESULTS_ADDED_MASS_ATTR['description']) +
                ' at ' + structure.H5_RESULTS_ADDED_MASS + ' with characteristics ' +
                str(dset))

    dset = utility.require_dataset(hdf5_data, structure.H5_RESULTS_RADIATION_DAMPING, result.radiation_damping.shape, dtype='f')
    utility.set_hdf5_attributes(dset, structure.H5_RESULTS_RADIATION_DAMPING_ATTR)
    dset[:, :, :] = result.radiation_damping
    logger.info('Saved ' + str(structure.H5_RESULTS_RADIATION_DAMPING_ATTR['description']) +
                ' at ' + structure.H5_RESULTS_RADIATION_DAMPING + ' with characteristics ' +
                str(dset))

    excitation_forces = result.diffraction_force + result.froudkrylov_force
    dset = utility.require_dataset(hdf5_data, structure.H5_RESULTS_EXCITATION_FORCES, excitation_forces.shape, dtype='F')
    utility.set_hdf5_attributes(dset, structure.H5_RESULTS_EXCITATION_FORCES_ATTR)
    dset[:, :, :] = excitation_forces
    logger.info('Saved ' + str(structure.H5_RESULTS_EXCITATION_FORCES_ATTR['description']) +
                ' at ' + structure.H5_RESULTS_EXCITATION_FORCES + ' with characteristics ' +
                str(dset))

    tec_file = utility.get_setting(settings.RADIATION_COEFFICIENTS_TEC_FILE, custom_config,
                                   'RADIATION_COEFFICIENTS_TEC_FILE')
    if tec_file:
        save_radiation_coefficients(result, tec_file)
        logger.info('Radiation coefficients successfully saved in tecplot format at ' +
                    str(tec_file))
    else:
        logger.info('Radiation coefficients tecplot format generation is disabled')

    tec_file = utility.get_setting(settings.DIFFRACTION_FORCE_TEC_FILE, custom_config,
                                   'DIFFRACTION_FORCE_TEC_FILE')

    if tec_file:
        save_diffraction_force(result, tec_file)
        logger.info('Diffraction forces successfully saved in tecplot format at ' +
                    str(tec_file))
    else:
        logger.info('Diffraction forces tecplot format generation is disabled')

    tec_file = utility.get_setting(settings.EXCITATION_FORCE_TEC_FILE, custom_config,
                                   'EXCITATION_FORCE_TEC_FILE')
    if tec_file:
        save_excitation_force(result, tec_file)
        logger.info('Excitation forces successfully saved in tecplot format at ' +
                    str(tec_file))
    else:
        logger.info('Excitation forces tecplot format generation is disabled')

    irf = get_irf(hdf5_data, result)

    if irf.switch == 1:
        irf = compute_irf(result, irf)
        # Saving to hdf5 file
        dset = utility.require_dataset(hdf5_data, structure.H5_RESULTS_ADDED_MASS_INFINITE, irf.added_mass.shape, dtype='f')
        utility.set_hdf5_attributes(dset, structure.H5_RESULTS_ADDED_MASS_INFINITE_ATTR)
        dset[:, :] = irf.added_mass

        tec_file = utility.get_setting(settings.IRF_TEC_FILE, custom_config,
                                       'IRF_TEC_FILE')
        if tec_file:
            save_irf(irf, tec_file)
            logger.info('IRF successfully saved in tecplot format at ' +
                        str(tec_file))
        else:
            logger.info('IRF tecplot format generation is disabled')
    else:
        logger.info('IRF computation is disabled')

    raos = np.zeros((result.n_integration, result.n_w, result.n_beta), dtype='F')
    raos = compute_raos(raos, result)

    tec_file = utility.get_setting(settings.WAVE_FIELD_TEC_FILE, custom_config,
                                   'WAVE_FIELD_TEC_FILE')

    dset = hdf5_data.get(structure.H5_SOLVER_USE_HIGHER_ORDER)
    utility.check_dataset_type(dset,
                               name=str(structure.H5_SOLVER_USE_HIGHER_ORDER_ATTR['description']),
                               location=structure.H5_SOLVER_USE_HIGHER_ORDER)
    use_higher_order = dset[0]

    dset = hdf5_data.get(structure.H5_SOLVER_USE_DIPOLES_IMPLEMENTATION)
    utility.check_dataset_type(dset,
                               name=str(structure.H5_SOLVER_USE_DIPOLES_IMPLEMENTATION_ATTR['description']),
                               location=structure.H5_SOLVER_USE_DIPOLES_IMPLEMENTATION)
    use_dipoles_implementation = dset[0]

    dset = hdf5_data.get(structure.H5_SOLVER_REMOVE_IRREGULAR_FREQUENCIES)
    utility.check_dataset_type(dset,
                               name=str(structure.H5_SOLVER_REMOVE_IRREGULAR_FREQUENCIES_ATTR['description']),
                               location=structure.H5_SOLVER_REMOVE_IRREGULAR_FREQUENCIES)
    remove_irregular_frequencies = dset[0]

    # Wave Elevation computation and tec generation
    if result.n_theta < 1:
        tec_file = None
        logger.info('Wave elevation tecplot format generation is disabled because there is no directions (Kochin)')


    if tec_file:
        if use_higher_order != 1 and use_dipoles_implementation != 1 and remove_irregular_frequencies != 1:
            res = compute_wave_elevation(hdf5_data, environment, 0, 0, raos, result)
            save_wave_elevation(res['w'], res['etai'], res["etap"], res["eta"], res["x"], res["y"],
                            tec_file)
            logger.info('Wave elevation successfully saved in tecplot format at ' +
                        str(tec_file))
        else:
            logger.info('Wave elevation computation is not supported when higher order panel, ' +
                        'used diplome implementation or remove irregular frequencies are enabled.' +
                        ' Disabling it.')
    else:
        logger.info('Wave elevation tecplot format generation is disabled')

    #print(' -> All results successfully saved.\n')


def postprocess(custom_config):
    """
    Configure and then run the postprocessor

    Args:
        custom_config, dict The custom configuration dictionary
    """
    signature = __name__ + '.postprocess(custom_config)'
    logger = logging.getLogger(__name__)
    utility.log_entrance(logger, signature,
                        {'custom_config': custom_config})

    if not custom_config:
        custom_config = {}

    hdf5_file = utility.get_setting(settings.HDF5_FILE, custom_config, 'HDF5_FILE')
    utility.check_is_file(hdf5_file, 'The path to the hdf5 file configured by HDF5_FILE')

    #utility.validate_file(hdf5_file, 'HDF5_FILE')
    with h5py.File(hdf5_file, "a") as hdf5_db:
        run(hdf5_db, custom_config)

    utility.log_and_print(logger, 'The post processing results are saved in the hdf5 file '
                          + utility.get_abs(hdf5_file))

    utility.log_exit(logging.getLogger(__name__), signature, [None])


def run_as_process(custom_config, queue):
    utility.setup_subprocess_logging(queue, logging.getLogger())
    return postprocess(custom_config)

if __name__ == '__main__':
    utility.setup_logging(default_conf_path=settings.LOGGING_CONFIGURATION_FILE, logging_path=settings.LOG_FILE)
    try:
        postprocess({})
        print('Post processing successfully completed.' + '\n')
    except Exception as e:
        print('There was an error when running the application. Check the log file for more details' + '\n')
        # exc_info=True means the stack trace will be printed automatically
        logging.getLogger(__name__).error('Program halted due to a fatal error whose detail is as follow: ',
                                          exc_info=True)