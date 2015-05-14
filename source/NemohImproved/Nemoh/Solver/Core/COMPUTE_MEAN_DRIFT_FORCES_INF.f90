!--------------------------------------------------------------------------------------
!
!Copyright (C) 2014 TopCoder Inc., All Rights Reserved.
!
!--------------------------------------------------------------------------------------

!   This module computes the means drift forces using the far-field momentum formulation for a infinite depth
!   case according to equation 4.1 and 4.2 of R1
!
!  R1 Computation of Higher-Order Hydrodynamic Forces on Ships and Offshore Structures in Waves
!             http://dspace.mit.edu/bitstream/handle/1721.1/79979/42664020.pdf?sequence=1
!
! Contest Drift forces and QTF Implementation of Nemoh
!
! Changes in version 1.1 (Implementation of Higher Order Panel Methods)
!       Added COMMON_TYPE module as dependency
!
! Changes in version 1.2 (Dipoles Implementation in NEMOH)
!       Switch from ParamsCommon to ParamsCommonInf
!
!   @author yedtoss
!   @version 1.2

module COMPUTE_MEAN_DRIFT_FORCES_INF

    USE COMMON_TYPE
    USE COM_VAR
    USE UTILITY
    USE KOCHIN

    implicit none

contains

    function driftintegral_fx(theta, param)
        ! computes abs(H(theta))**2*cos(theta)
        USE KOCHIN
        REAL, INTENT(IN) :: theta
        REAL :: driftintegral_fx
        TYPE(ParamsCommonInf) :: param
        COMPLEX :: HKochin

        CALL COMPUTE_KOCHIN(param%kwave, theta, HKochin, param%SolverVar)

        driftintegral_fx = (ABS(HKochin)**2)*COS(theta)

    end function driftintegral_fx

    function driftintegral_fy(theta, param)
        ! computes abs(H(theta))**2*sin(theta)

        USE KOCHIN
        REAL, INTENT(IN) :: theta
        REAL :: driftintegral_fy
        TYPE(ParamsCommonInf) :: param
        COMPLEX :: HKochin

        CALL COMPUTE_KOCHIN(param%kwave, theta, HKochin, param%SolverVar)

        driftintegral_fy = (ABS(HKochin)**2)*SIN(theta)

    end function driftintegral_fy


    subroutine compute_mean_drift_forces_infinite(param, drift_forces)
        ! This is the main subrouting computing the drift_forces
        ! It use romberg_trap to integrate the function driftintegral_fy and driftintegral_fx

        REAL :: PI
        TYPE(ParamsCommonInf) :: param
        REAL, DIMENSION(2) :: drift_forces !
        COMPLEX :: HKochin

        PI=4.*ATAN(1.)

        drift_forces(1) = (RHO* (param%kwave**2)/(8*PI))*romberg_trap(driftintegral_fx, 0., 2*PI, 1e-3, 20, param)
        drift_forces(2) = (RHO* (param%kwave**2)/(8*PI))*romberg_trap(driftintegral_fy, 0., 2*PI, 1e-3, 20, param)

        CALL COMPUTE_KOCHIN(param%kwave, PI+param%beta, HKochin, param%SolverVar)

        drift_forces(1) = drift_forces(1) + 0.5*RHO*param%w*param%AKH*cos(param%beta)*REAL(HKochin)
        drift_forces(2) = drift_forces(2) + 0.5*RHO*param%w*param%AKH*sin(param%beta)*REAL(HKochin)

    end subroutine compute_mean_drift_forces_infinite



end module COMPUTE_MEAN_DRIFT_FORCES_INF
