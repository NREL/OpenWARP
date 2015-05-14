!--------------------------------------------------------------------------------------
!
!Copyright (C) 2014 TopCoder Inc., All Rights Reserved.
!
!--------------------------------------------------------------------------------------

!   This module computes the mean yaw moment using the far-field momentum formulation for a infinite depth
!   case according to equation 4.6 of R1
!
!  R1 Computation of Higher-Order Hydrodynamic Forces on Ships and Offshore Structures in Waves
!             http://dspace.mit.edu/bitstream/handle/1721.1/79979/42664020.pdf?sequence=1
!
! Contest Drift forces and QTF Implementation of Nemoh
!
!
! Changes in version 1.1 (Implementation of Higher Order Panel Methods)
!       Added COMMON_TYPE module as dependency
!
! Changes in version 1.2 (Dipoles Implementation in NEMOH)
!       Switch from ParamsCommon to ParamsCommonInf
!
!   @author yedtoss
!   @version 1.2

module COMPUTE_YAW_MOMENT_INF

    USE COMMON_TYPE
    USE COM_VAR
    USE UTILITY
    USE KOCHIN

    implicit none

contains

    complex function diff_func_yaw_infinite(theta, param)
        ! Simple wrapper to the kochin function
        USE KOCHIN

        REAL, INTENT(IN) :: theta
        COMPLEX:: fx
        TYPE(ParamsCommonInf) :: param

        CALL COMPUTE_KOCHIN(param%kwave, theta, fx, param%SolverVar)
        diff_func_yaw_infinite = fx

    end function diff_func_yaw_infinite

    complex function yawintegral_mz(theta, param)
        ! Derivative of the kochin function multiplied the the conjugate of the value of kochin function at theta

        USE KOCHIN
        REAL, INTENT(IN) :: theta
        TYPE(ParamsCommonInf) :: param
        COMPLEX :: HKochin

        CALL COMPUTE_KOCHIN(param%kwave, theta, HKochin, param%SolverVar)

        yawintegral_mz = CONJG(HKochin)*deriv_richardson_complex(diff_func_yaw_infinite,theta, 1, 1e-3, param)

    end function yawintegral_mz




    subroutine compute_yaw_moment_infinite(param, yaw_moment)
        ! This is the main subrouting computing the yaw moment
        ! It use romberg_trap_complex to integrate the function yawintegralfin_mz
        ! and deriv_richardson_complex to find the derivative of the diff_func_yaw_finite

        REAL :: PI
        TYPE(ParamsCommonInf) :: param ! input parameter
        REAL :: yaw_moment, tmp  ! yaw_moment is the output

        PI=4.*ATAN(1.)

        tmp = (-0.5*aimag(deriv_richardson_complex(diff_func_yaw_infinite,PI+param%beta, 1, 1e-3, param)))

        yaw_moment = (-RHO* (param%kwave)/(8*PI))*AIMAG(romberg_trap_complex(yawintegral_mz, 0., 2*PI, 1e-3, 2, param))

        yaw_moment = yaw_moment + RHO*param%w*param%AKH/param%kwave*tmp

    end subroutine compute_yaw_moment_infinite
end module COMPUTE_YAW_MOMENT_INF
