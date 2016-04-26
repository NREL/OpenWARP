#!/usr/bin/env python
"""
This module defines nemoh custom types or models.

Changes in version 1.1 (OpenWarp - Add Logging Functionality)
       Added support for logging.
"""

import numpy as np

__author__ = "yedtoss"
__copyright__ = "Copyright (C) 2014-2016 TopCoder Inc. All rights reserved."
__version__ = "1.1"


class TCase:
    """
    This class represents the Case type
    """
    def __init__(self):
        """
        Initialize a new instance of this type
        Args:
            self: The class itself
        """
        self.body = 0
        self.i_case = 0
        self.mode = 0
        self.direction = np.zeros(3)
        self.axis = np.zeros(3)


class TMesh:
    """
    This class represents the Mesh type
    """
    def __init__(self, n_points=0, n_panels=0, n_bodies=0):
        """
        Initialize a new instance of this type
        Args:
            self: The class itself
            n_points: int, the number of points
            n_panels: int, the number of panels
            n_bodies: int, the number of bodies
        """
        self.i_sym = 0
        self.n_points = int(n_points)
        self.n_panels = int(n_panels)
        self.n_bodies = int(n_bodies)
        self.x = np.zeros((3, n_points))
        self.n = np.zeros((3, n_panels))
        self.xm = np.zeros((3, n_panels))
        self.p = np.zeros((4, n_panels))
        self.c_panel = np.zeros(n_panels)
        self.a = np.zeros(n_panels)
        self.last_panel = np.zeros(n_bodies)
        self.cg = np.zeros((3, n_bodies))

    def __str__(self):
        sym_msg = 'no symmetry'
        if self.i_sym:
            sym_msg = 'symmetry'

        return ('Mesh of ' + str(self.n_points) + ' points and ' + str(self.n_panels) + ' panels '
                + ' with ' + sym_msg)


class TEnvironment:
    """
    This class represents the Environment type
    """
    def __init__(self):
        """
        Initialize a new instance of this type
        Args:
            self: The class itself
        """
        self.rho = 0
        self.g = 0
        self.depth = 0
        self.x_eff = 0
        self.y_eff = 0

    def __str__(self):
        depth_msg = 'infinite water depth'
        if self.depth > 0:
            depth_msg = 'water depth of ' + str(self.depth)
        return ('Environment with ' + depth_msg + ', gravity: ' + str(self.g)
                + ', sea water density: ' + str(self.rho)
                + ' wave measurements points coordinate: (' + str(self.x_eff)
                + ', ' + str(self.y_eff) + ')')


class TResult:
    """
    This class represents the Result type
    """
    def __init__(self, n_w=0, n_radiation=0, n_integration=0, n_theta=0, n_beta=0):
        """
        Initialize a new instance of this type
        Args:
            self: The class itself
            n_w: int, the length of omega
            n_radiation: int the number of radiations
            n_integration: int the number of integration
            n_theta: int, the length of theta
            n_beta: int, the length of beta
        """
        self.n_w = n_w
        self.n_radiation = n_radiation
        self.n_integration = n_integration
        self.n_theta = n_theta
        self.n_beta = n_beta
        self.w = np.zeros(n_w, dtype='f')
        self.idx_force = np.zeros((n_integration, 3), dtype='i')
        self.idx_radiation = np.zeros((n_integration, 3), dtype='i')
        self.beta = np.zeros(n_beta, dtype='f')
        self.diffraction_force = np.zeros((n_w, n_beta, n_integration), dtype='F')
        self.froudkrylov_force = np.zeros((n_w, n_beta, n_integration), dtype='F')
        self.added_mass = np.zeros((n_w, n_radiation, n_integration), dtype='f')
        self.radiation_damping = np.zeros((n_w, n_radiation, n_integration), dtype='f')
        self.theta = np.zeros((n_theta, ), dtype='f')
        self.hkochin_diffraction = np.zeros((n_w, n_beta, n_theta), dtype='F')
        self.hkochin_radiation = np.zeros((n_w, n_radiation, n_theta), dtype='F')

    def __str__(self):
        return ('Result with ' + str(self.n_w) + ' number of wave frequencies '
                + str(self.n_radiation) + ' radiations ' + str(self.n_integration) +
                ' integrations')


class TIRF:
    """
    This class represents the IRF type
    """
    def __init__(self, n_time=0, n_radiation=0, n_integration=0):
        """
        Initialize a new instance of this type
        Args:
            self: The class itself
            n_time:int, the length of time
            n_radiation: int, the number of radiation
            n_integration: int, the number of integration
        """
        self.switch = 0
        self.n_time = n_time
        self.n_radiation = n_radiation
        self.n_integration = n_integration
        self.time = np.zeros((n_time, ), dtype='f')
        self.k = np.zeros((n_time, n_radiation, n_integration), dtype='f')
        self.added_mass = np.zeros((n_radiation, n_integration), dtype='f')

    def __str__(self):
        return ('IRF with ' + str(self.n_time) + ' time steps'
                )