"""
This is the cython wrapper which call the fotran module compute_nemoh
to run the nemoh solver

Changes in Version 1.1 (Drift forces and QTF Implementation of Nemoh)
                  Added parameter to store drift forces and yaw moment

Changes in version 1.2 (Implementation of Higher Order Panel Methods):
    Added logic to send handle additionnal Higher order panel variables to be sent to
    Nemoh Fortran

Changes in version 1.3 (Hydrodynamic Data Exporter Assembly v1.0)
       Send additional mesh properties like center of buoyancy, volume displacement to Nemoh Solver
       so that it can fill them.

Changes in version 1.4 (Irregular Frequencies Assembly)
       Added logic to send handle additionnal variables to be sent to Nemoh Fortran.
       Those additional variables determine whether or not Irregular frequencies should be
       removed and if so, the panels which are in the interior of the free surface.

Changes in version 1.5 (OpenWarp - Add Logging Functionality)
       Added support for logging.

Changes in version 1.6 (OPENWARP - FIX WAVE FREQUENCY AND DIRECTION CRASH BUG):
    1. Releasing the python GIL before running the fortran subroutine.
       Done to avoid polluting python in case there is segmentation fault from
       the fortran subroutine
"""

from numpy import ndarray
from numpy cimport ndarray
import numpy as np
cimport numpy as np
import sys
from ctypes import byref
from ctypes import c_char_p

__author__ = "yedtoss"
__copyright__ = "Copyright (C) 2014-2016 TopCoder Inc. All rights reserved."
__version__ = "1.6"


cdef extern:

    void compute_nemoh(float* rho, float * g, float * depth, float* xeff, float * yeff, float * zeff,
                        int* indiq_solver, int* max_iterations, int* restart_param, float* tol_gmres,
                        int** mesh_P, float** mesh_X, int* mesh_cPanel, float **mesh_XM,
                        float **mesh_N, float*  mesh_A,
                        int* mesh_Isym, int* Npoints, int* Npanels,
    int* Nbodies, int* Nproblems, int* NBCPanels, float complex ** bc_NormalVelocity, float* bc_Omega, int* bc_Switch_Potential,
    int*  bc_Switch_FreeSurface, int*  bc_Switch_Kochin, int*  bc_Switch_Type,
    float ** NDS, int* Nintegration, float * Theta,  int * Ntheta,
    int * NFSPoints, int * NFSPanels, float** meshfs_X, int** meshfs_P,
    float complex ** out_PHI, float complex ** out_PRESSURE, float complex ** out_HKochin, float **line,
    float**out_potential, int* n_potentials, int* n_tabulatedx, int* n_tabulatedz, int* n_points_simpson,
    float*** drift_forces, float** yaw_moment, int* FastInfluence_Switch, int* use_higher_order,
    int* n_beta, int* n_radiation, float** rad_case, float* beta, int* num_panel_higher_order, int* b_spline_order, 
    int* is_thin_body, int*use_dipoles_implementation, int* compute_yaw_moment, int* compute_drift_forces,
    float** center_buoyancy, float* displacement, float* waterplane_area, float***stifness, 
    int* is_interior_domain, int*remove_irregular_frequencies, int*log_level) nogil



def run_solver(data):
    """
    Python module to run the fortran nemoh solver
    Args:
        data: dict, dictionary containing the parameters to send to the fortran solver
    """
    cdef ndarray[float, mode="fortran", ndim=2] mesh_x = data["mesh_x"].astype(copy=False, dtype='f', order='F')
    cdef ndarray[int, mode="fortran", ndim=2] mesh_p = data["mesh_p"].astype(copy=False, dtype='i', order='F')
    cdef ndarray[int, mode="fortran", ndim=1] mesh_cpanel = data["mesh_cpanel"].astype(copy=False, dtype='i', order='F')
    cdef ndarray[float, mode="fortran", ndim=2] mesh_xm = data["mesh_xm"].astype(copy=False, dtype='f', order='F')
    cdef ndarray[float, mode="fortran", ndim=2] mesh_n = data["mesh_n"].astype(copy=False, dtype='f', order='F')
    cdef ndarray[float, mode="fortran", ndim=1] mesh_a = data["mesh_a"].astype(copy=False, dtype='f', order='F')
    cdef ndarray[float complex, mode="fortran", ndim=2] bc_normal_velocity = data["bc_normal_velocity"].astype(copy=False, dtype='F', order='F')
    cdef ndarray[float, mode="fortran", ndim=1] bc_omega = data["bc_omega"].astype(copy=False, dtype='f', order='F')
    cdef ndarray[int, mode="fortran", ndim=1] bc_switch_potential = data["bc_switch_potential"].astype(copy=False, dtype='i', order='F')
    cdef ndarray[int, mode="fortran", ndim=1] bc_switch_freesurface = data["bc_switch_freesurface"].astype(copy=False, dtype='i', order='F')
    cdef ndarray[int, mode="fortran", ndim=1] bc_switch_kochin = data["bc_switch_kochin"].astype(copy=False, dtype='i', order='F')
    cdef ndarray[int, mode="fortran", ndim=1] bc_switch_type = data["bc_switch_type"].astype(copy=False, dtype='i', order='F')
    cdef ndarray[float, mode="fortran", ndim=2] nds = data["nds"].astype(copy=False, dtype='f', order='F')
    cdef ndarray[float, mode="fortran", ndim=1] theta = data["theta"].astype(copy=False, dtype='f', order='F')
    cdef ndarray[float, mode="fortran", ndim=2] meshfs_x = data["meshfs_x"].astype(copy=False, dtype='f', order='F')
    cdef ndarray[int, mode="fortran", ndim=2] meshfs_p = data["meshfs_p"].astype(copy=False, dtype='i', order='F')
    cdef ndarray[float complex, mode="fortran", ndim=2] out_phi = data["out_phi"].astype(copy=False, dtype='F', order='F')
    cdef ndarray[float complex, mode="fortran", ndim=2] out_pressure = data["out_pressure"].astype(copy=False, dtype='F', order='F')
    cdef ndarray[float complex, mode="fortran", ndim=2] out_hkochin = data["out_hkochin"].astype(copy=False, dtype='F', order='F')
    cdef ndarray[float, mode="fortran", ndim=2] line = data["line"].astype(copy=False, dtype='f', order='F')
    cdef ndarray[float, mode="fortran", ndim=2] out_potential = data["out_potential"].astype(copy=False, dtype='f', order='F')
    cdef ndarray[float, mode="fortran", ndim=3] drift_forces = data["drift_forces"].astype(copy=False, dtype='f', order='F')
    cdef ndarray[float, mode="fortran", ndim=2] yaw_moment = data["yaw_moment"].astype(copy=False, dtype='f', order='F')
    cdef ndarray[int, mode="fortran", ndim=1] fast_influence_switch = data["fast_influence_switch"].astype(copy=False, dtype='i', order='F')
    cdef ndarray[float, mode="fortran", ndim=2] rad_case = data["rad_case"].astype(copy=False, dtype='f', order='F')
    cdef ndarray[float, mode="fortran", ndim=1] beta = data["beta"].astype(copy=False, dtype='f', order='F')
    cdef ndarray[int, mode="fortran", ndim=1] is_thin_body = data["is_thin_body"].astype(copy=False, dtype='i', order='F')

    cdef ndarray[float, mode="fortran", ndim=2] center_buoyancy = data["center_buoyancy"].astype(copy=False, dtype='f', order='F')
    cdef ndarray[float, mode="fortran", ndim=1] displacement = data["displacement"].astype(copy=False, dtype='f', order='F')
    cdef ndarray[float, mode="fortran", ndim=1] waterplane_area = data["waterplane_area"].astype(copy=False, dtype='f', order='F')
    cdef ndarray[float, mode="fortran", ndim=3] stifness = data["stifness"].astype(copy=False, dtype='f', order='F')
    cdef ndarray[int, mode="fortran", ndim=1] is_interior_domain = data["is_interior_domain"].astype(copy=False, dtype='i', order='F')

    cdef int nfs_panels = data["nfs_panels"]
    cdef int nfs_points = data["nfs_points"]
    cdef int n_theta = data["n_theta"]
    cdef int n_integration = data["n_integration"]
    cdef int nbc_panels = data["nbc_panels"]
    cdef int n_problems = data["n_problems"]
    cdef int n_bodies = data["n_bodies"]
    cdef int n_panels = data["n_panels"]
    cdef int n_points = data["n_points"]
    cdef int i_sym = data["i_sym"]
    cdef int n_potentials = data["n_potentials"]
    cdef float tol_gmres = data["tol_gmres"]
    cdef int restart_param = data["restart_param"]
    cdef int max_iterations = data["max_iterations"]
    cdef int indiq_solver = data["indiq_solver"]
    cdef float zeff = data["zeff"]
    cdef float yeff = data["yeff"]
    cdef float xeff = data["xeff"]
    cdef float depth = data["depth"]
    cdef float g = data["g"]
    cdef float rho = data["rho"]
    cdef int n_tabulatedx = data["n_tabulatedx"]
    cdef int n_tabulatedz = data["n_tabulatedz"]
    cdef int n_points_simpson = data["n_points_simpson"]
    cdef int use_higher_order = data["use_higher_order"]
    cdef int n_radiation = data["n_radiation"]
    cdef int n_beta = data["n_beta"]
    cdef int num_panel_higher_order = data["num_panel_higher_order"]
    cdef int b_spline_order = data["b_spline_order"]
    cdef int use_dipoles_implementation = data["use_dipoles_implementation"]
    cdef int compute_yaw_moment = data["compute_yaw_moment"]
    cdef int compute_drift_forces = data["compute_drift_forces"]
    cdef int remove_irregular_frequencies = data["remove_irregular_frequencies"]
    cdef int log_level = data["log_level"]
    
    

    with nogil:
        compute_nemoh(&rho, &g, &depth, &xeff, &yeff, &zeff,
        &indiq_solver, &max_iterations, &restart_param, &tol_gmres,
        <int**>mesh_p.data, <float**>mesh_x.data, <int*>mesh_cpanel.data, <float**>mesh_xm.data, <float**>mesh_n.data, <float*>mesh_a.data, &i_sym, &n_points,
        &n_panels, &n_bodies, &n_problems, &nbc_panels,
        <float complex**>(bc_normal_velocity.data), <float*>bc_omega.data, <int*>bc_switch_potential.data,
        <int*>bc_switch_freesurface.data, <int*>bc_switch_kochin.data, <int*>bc_switch_type.data, <float**>(nds.data), &n_integration, <float*>theta.data,
        &n_theta, &nfs_points, &n_panels, <float**>(meshfs_x.data), <int**>(meshfs_p.data), <float complex**>(out_phi.data),
        <float complex**>(out_pressure.data), <float complex**>(out_hkochin.data), <float**>line.data, <float**>out_potential.data,
        &n_potentials, &n_tabulatedx, &n_tabulatedz, &n_points_simpson, <float***>drift_forces.data, <float**>yaw_moment.data, <int*> fast_influence_switch.data,
        &use_higher_order, &n_beta, &n_radiation, <float**> rad_case.data, <float*> beta.data, &num_panel_higher_order, &b_spline_order, 
        <int*> is_thin_body.data, &use_dipoles_implementation, &compute_yaw_moment, &compute_drift_forces,
        <float**>center_buoyancy.data, <float*>displacement.data, <float*>waterplane_area.data,
        <float***>stifness.data, <int*>is_interior_domain.data, &remove_irregular_frequencies, &log_level)

    return 0