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

import preprocessor


from models import TResult
from models import TIRF

__author__ = "yedtoss, TCSASSEMBLER"
__copyright__ = "Copyright (C) 2014-2015 TopCoder Inc. All rights reserved."
__version__ = "1.5"


def get_irf(hdf5_data, result):
    """
    Gets the irf from the hdf5 file
    Args:
        hdf5_data: object, the hdf5 opened file
        result: object the hydrodynamic coefficients cases
    Returns:
        the irf
    """

    switch = hdf5_data.get(structure.H5_COMPUTE_IRF)[0]
    time_step = hdf5_data.get(structure.H5_IRF_TIME_STEP)[0]
    duration = hdf5_data.get(structure.H5_IRF_DURATION)[0]
    irf = TIRF()

    if switch == 1:
        irf = TIRF(int(duration/time_step), result.n_radiation, result.n_integration)
        for i in range(irf.n_time):
            irf.time[i] = i*time_step
    irf.switch = switch

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
    for i in range(irf.n_time):
        for j in range(result.n_radiation):
            for k in range(result.n_integration):
                irf.k[i,j,k] = 0.
                for l in range(result.n_w -1):

                    irf.k[i, j, k] += - 0.5* (result.w[l+1]-result.w[l]) * (result.radiation_damping[l,
                                                                                                                   j,
                     k]*np.cos(result.w[l]*irf.time[i]) + result.radiation_damping[l+1, j, k]*np.cos(result.w[l+1]*irf.time[i]))

                irf.k[i, j, k] = (irf.k[i, j, k] * 2)/np.pi



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

    return irf


def save_irf(irf, filename):
    """
    Saves the irf to a file in the tec format
    Args:
        irf: object, the irf
        filename: string, The path to the file where to save the irf
    """
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


def read_results(hdf5_data):
    """
    Read the hydrodynamic coefficients cases from the hdf5 file
    Args:
        hdf5_data: object, the hdf5 opened file
    Returns:
        the hydrodynamic coefficients cases
    """
    idx_force = hdf5_data.get(structure.H5_RESULTS_CASE_FORCE)
    n_integration = idx_force.shape[0]

    idx_radiation = hdf5_data.get(structure.H5_RESULTS_CASE_MOTION)
    n_radiation = idx_radiation.shape[0]

    beta = hdf5_data.get(structure.H5_RESULTS_CASE_BETA)
    n_beta = beta.shape[0]

    w = hdf5_data.get(structure.H5_RESULTS_CASE_W)
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

    return result


def save_radiation_coefficients(result, filename):
    """
    Saves the radiation coefficient to a file in tec format
    Args:
        result: object, the hydrodynamic coefficients cases
        filename: The path to the file where to save
    """
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


def save_diffraction_force(result, filename):
    """
    Saves the diffraction forces to a file in tec format
    Args:
        result: object, the hydrodynamic coefficients cases
        filename: The path to the file where to save
    """
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


def save_excitation_force(result, filename):
    """
    Saves the excitation forces to a file in tec format
    Args:
        result: object, the hydrodynamic coefficients cases
        filename: The path to the file where to save
    """
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

    nx = hdf5_data.get(structure.H5_FREE_SURFACE_POINTS_X)[0]
    ny = hdf5_data.get(structure.H5_FREE_SURFACE_POINTS_Y)[0]
    lx = hdf5_data.get(structure.H5_FREE_SURFACE_DIMENSION_X)[0]
    ly = hdf5_data.get(structure.H5_FREE_SURFACE_DIMENSION_Y)[0]

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
    kwave = utility.compute_wave_number(w, environment)

    for i in range(nx):
        for j in range(ny):
            r = np.sqrt((x[i] - environment.x_eff)**2 + (y[j] - environment.y_eff)**2)
            theta = np.arctan2(y[j]-environment.y_eff, x[i]-environment.x_eff)
            k = 0
            while (k < result.n_theta -1) and (result.theta[k+1] < theta):
                k += 1
            if k == result.n_theta:
                print(' Error: range of theta in Kochin coefficients is too small')
                sys.exit(1)

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

    return {"w": w, "x": x, "y": y, "eta": eta, "etai": etai, "etap": etap}


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


def run(hdf5_data, custom_config):
    """
    This function run the postprocessor
    Args:
        hdf5_data: object, the hdf5 opened file
        custom_config, dict The custom configuration dictionary
    """
    print('\n  -> Initialisation ...')

    try:
        environment = utility.read_environment(hdf5_data)
    except Exception as e:
        print('It looks like your hdf5 file is not correct. Please run ',
        'the preprocessor and the solver before running the postprocessor')
        sys.exit(1)

    result = read_results(hdf5_data)

    print('. Initialisation Done !\n')

    # Saving to hdf5 file
    dset = utility.require_dataset(hdf5_data, structure.H5_RESULTS_ADDED_MASS, result.added_mass.shape, dtype='f')
    utility.set_hdf5_attributes(dset, structure.H5_RESULTS_ADDED_MASS_ATTR)
    dset[:, :, :] = result.added_mass

    dset = utility.require_dataset(hdf5_data, structure.H5_RESULTS_RADIATION_DAMPING, result.radiation_damping.shape, dtype='f')
    utility.set_hdf5_attributes(dset, structure.H5_RESULTS_RADIATION_DAMPING_ATTR)
    dset[:, :, :] = result.radiation_damping

    excitation_forces = result.diffraction_force + result.froudkrylov_force
    dset = utility.require_dataset(hdf5_data, structure.H5_RESULTS_EXCITATION_FORCES, excitation_forces.shape, dtype='F')
    utility.set_hdf5_attributes(dset, structure.H5_RESULTS_EXCITATION_FORCES_ATTR)
    dset[:, :, :] = excitation_forces
    

    tec_file = utility.get_setting(settings.RADIATION_COEFFICIENTS_TEC_FILE, custom_config,
                                   'RADIATION_COEFFICIENTS_TEC_FILE')
    if tec_file:
        save_radiation_coefficients(result, tec_file)
        print('Radiation coefficients successfully saved.\n')

    tec_file = utility.get_setting(settings.DIFFRACTION_FORCE_TEC_FILE, custom_config,
                                   'DIFFRACTION_FORCE_TEC_FILE')
    if tec_file:
        save_diffraction_force(result, tec_file)
        print('Diffraction forces successfully saved.\n')

    tec_file = utility.get_setting(settings.EXCITATION_FORCE_TEC_FILE, custom_config,
                                   'EXCITATION_FORCE_TEC_FILE')
    if tec_file:
        save_excitation_force(result, tec_file)
        print('Excitation forces successfully saved.\n')

    
    irf = get_irf(hdf5_data, result)
    if not irf:
        print('It looks like your hdf5 file is not correct. Please run ',
        'the preprocessor and the solver before running the postprocessor')
        sys.exit(1)
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
            print('IRF successfully saved.\n')

    raos = np.zeros((result.n_integration, result.n_w, result.n_beta), dtype='F')
    raos = compute_raos(raos, result)

    

    tec_file = utility.get_setting(settings.WAVE_FIELD_TEC_FILE, custom_config,
                                   'WAVE_FIELD_TEC_FILE')
    if tec_file and hdf5_data.get(structure.H5_SOLVER_USE_HIGHER_ORDER)[0] != 1 and hdf5_data.get(structure.H5_SOLVER_USE_DIPOLES_IMPLEMENTATION)[0] != 1 and hdf5_data.get(structure.H5_SOLVER_REMOVE_IRREGULAR_FREQUENCIES)[0] != 1:
        res = compute_wave_elevation(hdf5_data, environment, 0, 0, raos, result)
        save_wave_elevation(res['w'], res['etai'], res["etap"], res["eta"], res["x"], res["y"],
                            tec_file)
        print('Wave elevation successfully saved.\n')

    print(' -> All results successfully saved.\n')


def postprocess(custom_config):
    """
    Configure and then run the postprocessor

    Args:
        custom_config, dict The custom configuration dictionary
    """

    if not custom_config:
        custom_config = {}

    hdf5_file = utility.get_setting(settings.HDF5_FILE, custom_config, 'HDF5_FILE')

    utility.validate_file(hdf5_file, 'HDF5_FILE')
    with h5py.File(hdf5_file, "a") as hdf5_db:
        run(hdf5_db, custom_config)

if __name__ == '__main__':
    postprocess({})