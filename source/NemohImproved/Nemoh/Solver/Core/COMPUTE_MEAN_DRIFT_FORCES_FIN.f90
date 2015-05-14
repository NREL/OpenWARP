!--------------------------------------------------------------------------------------
!
!Copyright (C) 2014 TopCoder Inc., All Rights Reserved.
!
!--------------------------------------------------------------------------------------

!   This module computes the means drift forces using the far-field momentum formulation for a finite depth
!   case according to equation 4.7 and 4.8 of R1
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

module COMPUTE_MEAN_DRIFT_FORCES_FIN

    USE COMMON_TYPE
    USE COM_VAR
    USE UTILITY
    USE KOCHIN

    implicit none

contains

    function driftintegralfin_fx(theta, param)
        ! computes abs(H(theta))**2*cos(theta)
        USE KOCHIN
        REAL, INTENT(IN) :: theta
        REAL :: driftintegralfin_fx
        TYPE(ParamsCommonInf) :: param
        COMPLEX :: HKochin

        CALL COMPUTE_KOCHIN(param%kwave, theta, HKochin, param%SolverVar)

        driftintegralfin_fx = (ABS(HKochin)**2)*COS(theta)

    end function driftintegralfin_fx

    function driftintegralfin_fy(theta, param)
        ! computes abs(H(theta))**2*sin(theta)

        USE KOCHIN
        REAL, INTENT(IN) :: theta
        REAL :: driftintegralfin_fy
        TYPE(ParamsCommonInf) :: param
        COMPLEX :: HKochin

        CALL COMPUTE_KOCHIN(param%kwave, theta, HKochin, param%SolverVar)

        driftintegralfin_fy = (ABS(HKochin)**2)*SIN(theta)

    end function driftintegralfin_fy


    subroutine compute_mean_drift_forces_finite(param, drift_forces)
        ! This is the main subrouting computing the drift_forces
        ! It use romberg_trap to integrate the function driftintegral_fy and driftintegral_fx

        REAL :: PI
        TYPE(ParamsCommonInf) :: param
        REAL, DIMENSION(2) :: drift_forces
        REAL :: const
        REAL :: m0
        REAL :: CgC
        REAL:: Z
        COMPLEX :: HKochin



        PI=4.*ATAN(1.)

        m0 = param%AMH

        Cgc = 0.5*(1 + 2*m0*Depth/sinh(2*m0*Depth))

        Z = (m0**2-param%kwave**2)/(param%kwave**2*Depth-m0**2*Depth-param%kwave)*cosh(m0*Depth)**2

        const = -(RHO*G*param%AKH**2/m0)*CgC

        drift_forces(1) =const*(Z**2/(4*PI)*romberg_trap(driftintegralfin_fx, 0., 2*PI, 1e-3, 20, param))
        drift_forces(2) = const*(Z**2/(4*PI)**romberg_trap(driftintegralfin_fy, 0., 2*PI, 1e-3, 20, param))

        CALL COMPUTE_KOCHIN(param%kwave, PI+param%beta, HKochin, param%SolverVar)

        drift_forces(1) = drift_forces(1) + (const)*Z*cos(param%beta)*aimag(HKochin)
        drift_forces(2) = drift_forces(2) + (const)*Z*sin(param%beta)*aimag(HKochin)

    end subroutine compute_mean_drift_forces_finite
end module COMPUTE_MEAN_DRIFT_FORCES_FIN
