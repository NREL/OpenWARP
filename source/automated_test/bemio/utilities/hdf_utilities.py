from __future__ import division

import numpy as np

import h5py

from bemio.data_structures.bem import HydrodynamicData

def combine_h5(files):
    '''Function to combine df5 5 hydrodynamic files

    Parameters:
        files: list
            A list of files to be combined
    '''
    hdf5_data = {}
    for i,f in enumerate(files):
        hdf5_data[i] = h5py.File(f,'r')

    return hdf5_data

def check_parameters(hdf5_data):
    pass

def create_hydro_data(hdf5_data):
    body1 = hdf5_data[0]['body1']
    simulation_paramters = hdf5_data[0]['simulation_parameters']
    bem_data = hdf5_data[0]['bem_data']

    hydro_data = HydrodynamicData()
    hydro_data.T = simulation_paramters['T'].value
    hydro_data.g = simulation_paramters['g'].value
    hydro_data.rho = simulation_paramters['rho'].value
    hydro_data.scaled = simulation_paramters['scaled'].value
    hydro_data.w = simulation_paramters['w'].value
    hydro_data.water_depth = simulation_paramters['water_depth'].value
    hydro_data.wave_dir = simulation_paramters['wave_dir'].value
    hydro_data.bem_code = bem_data['code'].vlaue
    hydro_data.am.all = body_1['hydro_coeffs/added_mass/all'].value
    hydro_data.am.inf_freq = body_1['hydro_coeffs/added_mass/inf_freq'].value
    hydro_data.ex.im = body_1['hydro_coeffs/excitation/im'].value
    hydro_data.ex.mag = body_1['hydro_coeffs/excitation/mag'].value
    hydro_data.ex.phase = body_1['hydro_coeffs/excitation/phase'].value
    hydro_data.ex.re = body_1['hydro_coeffs/excitation/re'].value
    hydro_data.k = body_1['hydro_coeffs/linear_restoring_stiffness'].value



    return hydro_data
