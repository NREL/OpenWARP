!--------------------------------------------------------------------------------------
!
!Copyright (C) 2014 TopCoder Inc., All Rights Reserved.
!
!--------------------------------------------------------------------------------------

!   This module computes the mean yaw moment using the far-field momentum formulation for a finite depth
!   case according to equation 4.9 of R1
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

module COMPUTE_YAW_MOMENT_FIN

    USE COMMON_TYPE
    USE COM_VAR
    USE UTILITY
    USE KOCHIN

    implicit none

contains

    complex function diff_func_yaw_finite(theta, param)
        ! Simple wrapper to the kochin function
        USE KOCHIN

        REAL, INTENT(IN) :: theta
        COMPLEX:: fx
        TYPE(ParamsCommonInf) :: param

        CALL COMPUTE_KOCHIN(param%kwave, theta, fx, param%SolverVar)
        diff_func_yaw_finite = fx

    end function diff_func_yaw_finite

    complex function yawintegralfin_mz(theta, param)
        ! Derivative of the kochin function multiplied the the conjugate of the value of kochin function at theta

        USE KOCHIN
        REAL, INTENT(IN) :: theta
        COMPLEX :: yawintegral_mz
        TYPE(ParamsCommonInf) :: param
        COMPLEX :: HKochin

        CALL COMPUTE_KOCHIN(param%kwave, theta, HKochin, param%SolverVar)

        yawintegralfin_mz = CONJG(HKochin)*deriv_richardson_complex(diff_func_yaw_finite,theta, 1, 1e-3, param)

    end function yawintegralfin_mz




    subroutine compute_yaw_moment_finite(param, yaw_moment)
        ! This is the main subrouting computing the yaw moment
        ! It use romberg_trap_complex to integrate the function yawintegralfin_mz
        ! and deriv_richardson_complex to find the derivative of the diff_func_yaw_finite

        REAL :: PI
        TYPE(ParamsCommonInf) :: param ! input parameter
        REAL :: yaw_moment ! output
        REAL :: const
        REAL :: m0
        REAL :: CgC
        REAL:: Z
        REAL:: tmp

        PI=4.*ATAN(1.)

        m0 = param%AMH

        Cgc = 0.5*(1 + 2*m0*Depth/sinh(2*m0*Depth))

        Z = (m0**2-param%kwave**2)/(param%kwave**2*Depth-m0**2*Depth-param%kwave)*cosh(m0*Depth)**2

        const = -(RHO*G*param%AKH**2/m0)*CgC

        yaw_moment = const*(Z**2/(4*PI)*AIMAG(romberg_trap_complex(yawintegralfin_mz, 0., 2*PI, 1e-3, 2, param)))

        tmp = const*(-0.5*REAL(deriv_richardson_complex(diff_func_yaw_finite,PI+param%beta, 1, 1e-3, param)))

        yaw_moment = yaw_moment + tmp

    end subroutine compute_yaw_moment_finite
end module COMPUTE_YAW_MOMENT_FIN
