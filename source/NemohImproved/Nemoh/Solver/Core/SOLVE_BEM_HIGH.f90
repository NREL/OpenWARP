!--------------------------------------------------------------------------------------
!
!Copyright (C) 2014 TopCoder Inc., All Rights Reserved.
!
!--------------------------------------------------------------------------------------

!   This module solves a boundary element method for a finite or infinite depth using
!   a three dimensional higher order method. The higher order method used is based on
!   B-Spline as described in the references given below.
!   The influence matrix is setup is solved using an SVD (Singular Value Decomposition) solver.
!   More specific description can be read in section 15.5 of reference R2.
!   Here the green function and derivative are computed using a series approximation and a tabulation
!   as previously done for the lower order panel method in Nemoh.
!   After solving the influence matrix, the potential at the boundary and at the given points in the fluid domain is
!   computed. The pressure at the centroid of each patch are computed and the forces are derived.
!   The main subroutine of this module is SOLVE_POTENTIAL_HIGH
!
!                       REFERENCES
!   R1:  A three dimensional higher order panel method based on
!        B-splines. http://dspace.mit.edu/bitstream/handle/1721.1/11127/34279481.pdf?sequence=1
!
!   R2 Wamit User Manual http://www.wamit.com/manualupdate/V70_manual.pdf
!
!   R3 A higher order panel method for large-amplitude simulation of bodies in waves.
!      Ph. D. Thesis http://dspace.mit.edu/bitstream/handle/1721.1/9708/42663926.pdf?sequence=1
!
!   R4 Fast Multipole Boundary Element Method: Theory and Applications in Engineering
!   https://books.google.bj/books?id=1p4SCnM5UIYC
!   http://www.amazon.com/Fast-Multipole-Boundary-Element-Method-ebook/dp/B002TRJ076
!
!   Contest Implementation of Higher Order Panel Methods version 1.0
!
!   @author yedtoss
!   @version 1.0

MODULE SOLVE_BEM_HIGH

    USE COM_VAR
    USE ELEMENTARY_FNS
    USE UTILITY
    USE gaussm3
    USE COMMON_TYPE
    USE GREEN_FUNCTION
    USE OUTPUT
    USE M_SOLVER
    implicit none

contains

    complex FUNCTION influence_dik(u, v, param)
        ! This function computes the function to be integrated according to equation 15.32 of reference R2
        ! It does not compute the integral. So to compute the formula of equation 15.32, it is expected to call
        ! a subrouting which will integrate this function. For example:
        ! qgss2d_simple(influence_dik, -1., 1., -1., 1., 3, param)
        ! Arguments
        ! u,v coordinate in parametric space of the integral
        ! param: common variables that contains useful variables
        !        param%j is the index of the source patch (corresponding to uf)
        !        param%k is the index of the basis function in the v direction for the source patch
        !        param%i is the index of the basis function in the u direction for the source patch
        !        param%jj is the index of the field patch (corresponding to u)
        !        param%kk is the index of the basis function in the v direction for the field patch
        !        param%ii is the index of the basis function in the u direction for the field patch

        REAL, INTENT(IN) :: u, v
        TYPE(ParamsCommonInf) :: param

        influence_dik = cmplx(GET_BSPLINE(u, param%SolverVar%KNOTS_U(param%j), param%i),0)

        influence_dik = influence_dik * GET_BSPLINE(v, param%SolverVar%KNOTS_V(param%j), param%k)

        IF(param%jj /= -1 .and. param%ii /= -1 .and. param%kk /= -1) THEN
            influence_dik = influence_dik * GET_BSPLINE(u, param%SolverVar%KNOTS_U(param%jj), param%ii)
            influence_dik = influence_dik * GET_BSPLINE(v, param%SolverVar%KNOTS_V(param%jj), param%kk)
        END IF


    END FUNCTION influence_dik


    complex FUNCTION influence_si(uf, vf, param)
        ! This function computes the function of the outer integral (including the inner integral) in equation 15.34
        ! of reference R2
        ! It does not compute the integral. So to compute the whole integral of equation 15.34, it is expected to call
        ! a subrouting which will integrate this function. For example:
        ! qgss2d_simple(influence_si, -1., 1., -1., 1., 3, param)
        ! Arguments
        ! uf,vf coordinate in parametric space of the outer integral (source point)
        ! param: common variables that contains useful variables
        !        param%j is the index of the source patch (corresponding to uf)
        !        param%k is the index of the basis function in the v direction for the source patch
        !        param%i is the index of the basis function in the u direction for the source patch
        !        param%jj is the index of the field patch (corresponding to u)
        !        param%kk is the index of the basis function in the v direction for the field patch
        !        param%ii is the index of the basis function in the u direction for the field patch

        REAL, INTENT(IN) :: uf, vf
        TYPE(ParamsCommonInf) :: param

        param%uf = uf
        param%vf = vf

        influence_si = cmplx(GET_BSPLINE(uf, param%SolverVar%KNOTS_U(param%j), param%i),0)

        influence_si = influence_si * GET_BSPLINE(vf, param%SolverVar%KNOTS_V(param%j), param%k)

        influence_si = influence_si * qgss2d_simple(influence_ugreen, -1., 1., -1., 1., 3, param)


    END FUNCTION influence_si


    complex FUNCTION influence_bigdik(uf, vf, param)
        ! This function computes the function of the outer integral (including the inner integral) in equation 15.33
        ! of reference R2
        ! It does not compute the integral. So to compute the whole integral of equation 15.33, it is expected to call
        ! a subrouting which will integrate this function. For example:
        ! qgss2d_simple(influence_bigdik, -1., 1., -1., 1., 3, param)
        ! Arguments
        ! uf,vf coordinate in parametric space of the outer integral (source point)
        ! param: common variables that contains useful variables
        !        param%j is the index of the source patch (corresponding to uf)
        !        param%k is the index of the basis function in the v direction for the source patch
        !        param%i is the index of the basis function in the u direction for the source patch
        !        param%jj is the index of the field patch (corresponding to u)
        !        param%kk is the index of the basis function in the v direction for the field patch
        !        param%ii is the index of the basis function in the u direction for the field patch

        REAL, INTENT(IN) :: uf, vf
        TYPE(ParamsCommonInf) :: param

        param%uf = uf
        param%vf = vf
        influence_bigdik = qgss2d_simple(influence_udgreen, -1., 1., -1., 1., 3, param)
        influence_bigdik = influence_bigdik * GET_BSPLINE(uf, param%SolverVar%KNOTS_U(param%j), param%i)
        influence_bigdik = influence_bigdik * GET_BSPLINE(vf, param%SolverVar%KNOTS_V(param%j), param%k)

    END FUNCTION influence_bigdik


    complex FUNCTION influence_ugreen(u, v, param)
        ! This function computes the function of the inner integral in equation 15.34 of reference R2
        ! It does not compute the integral. So to compute the whole integral of equation 15.34, it is expected to call
        ! a subrouting which will integrate this function. For example:
        ! qgss2d_simple(influence_ugreen, -1., 1., -1., 1., 3, param)
        ! Note that this function simultaneously computes the function of the integral at the right hand side
        ! of equation 15.29 of reference R2
        ! Arguments
        ! u,v coordinate in parametric space of the inner integral (field point)
        ! param: common variables that contains useful variables
        !        param%j is the index of the source patch (corresponding to uf)
        !        param%k is the index of the basis function in the v direction for the source patch
        !        param%i is the index of the basis function in the u direction for the source patch
        !        param%jj is the index of the field patch (corresponding to u)
        !        param%kk is the index of the basis function in the v direction for the field patch
        !        param%ii is the index of the basis function in the u direction for the field patch

        REAL, INTENT(IN) :: u, v
        REAL::  uf, vf
        INTEGER:: panel_idxf
        TYPE(ParamsCommonInf) :: param
        complex :: tmp
        REAL:: xx, yy, zz, XPO, YPO, ZPO

        panel_idxf = param%j

        IF(param%use_cartesian_coordinate == 1) THEN
            XPO = param%xx
            YPO = param%yy
            ZPO = param%zz
        ELSE
            uf = param%uf
            vf = param%vf
            CALL transform_uv_to_cartesian(panel_idxf, uf ,vf, XPO, YPO, ZPO)

        END IF

        CALL transform_uv_to_cartesian(param%jj, u ,v, xx, yy, zz)

        influence_ugreen = cmplx(compute_jacobian(param%jj, u, v), 0)
        tmp = compute_green_wave_part(param, u, v, XPO, YPO, ZPO, param%jj)
        tmp = tmp + compute_green_rankine_part(param, u, v, XPO, YPO, ZPO, param%jj)

        influence_ugreen = influence_ugreen * tmp * COMPUTE_NORMAL_POTENTIAL(param, param%jj, xx, yy, zz)
    end FUNCTION influence_ugreen





    complex FUNCTION influence_udgreen(u, v, param)
        ! This function computes the function of the inner integral in equation 15.33 of reference R2
        ! It does not compute the integral. So to compute the whole integral of equation 15.33, it is expected to call
        ! a subrouting which will integrate this function. For example:
        ! qgss2d_simple(influence_udgreen, -1., 1., -1., 1., 3, param)
        ! Note that this function simultaneously computes the function of the integral of the last term
        ! at the left hand side of equation 15.29 of reference R2 if param%compute_potential is not 0
        ! Arguments
        ! u,v coordinate in parametric space of the inner integral (field point)
        ! param: common variables that contains useful variables
        !        param%j is the index of the source patch (corresponding to uf)
        !        param%k is the index of the basis function in the v direction for the source patch
        !        param%i is the index of the basis function in the u direction for the source patch
        !        param%jj is the index of the field patch (corresponding to u)
        !        param%kk is the index of the basis function in the v direction for the field patch
        !        param%ii is the index of the basis function in the u direction for the field patch
        REAL, INTENT(IN) :: u, v
        REAL::  uf, vf, XPO, YPO, ZPO
        INTEGER:: panel_idxf
        TYPE(ParamsCommonInf) :: param
        complex :: tmp


        panel_idxf = param%j

        IF(param%use_cartesian_coordinate == 1) THEN
            XPO = param%xx
            YPO = param%yy
            ZPO = param%zz
        ELSE
            uf = param%uf
            vf = param%vf
            CALL transform_uv_to_cartesian(panel_idxf, uf ,vf, XPO, YPO, ZPO)

        END IF

        IF(param%compute_potential == 0) THEN
            influence_udgreen = cmplx(GET_BSPLINE(u, param%SolverVar%KNOTS_U(param%jj), param%ii), 0)
        ELSE
            influence_udgreen = COMPUTE_POTENTIAL_BOUNDARY_PATCH(param%jj, u, v, param%B, param)
        END IF
        influence_udgreen = influence_udgreen * GET_BSPLINE(v, param%SolverVar%KNOTS_V(param%jj), param%kk)
        influence_udgreen = influence_udgreen * compute_jacobian(param%jj, u, v)
        tmp = compute_dgreen_wave_part(param, u, v, XPO, YPO, ZPO, param%jj)
        tmp = tmp + compute_dgreen_rankine_part(param, u, v, XPO, YPO, ZPO, param%jj)

        influence_udgreen = influence_udgreen * tmp
    END FUNCTION



    SUBROUTINE transform_cartesian_to_uv(panel_idx, xx, yy, zz, u, v)
        ! This subroutine converts a cartesian coordinate to a parametric coordinate for a panel
        ! Currently it is a mock and expected to be change in later contest if needed.
        ! The current implementation does not use it. This function could be useful if we want
        ! to compute the potential at the boundary exactly without doing any integration (as currently).
        REAL :: u,v,xx,yy,zz, r, theta, rmin, rmax, tmp
        INTEGER:: panel_idx, M(4), I
        REAL, PARAMETER :: PI=4.*ATAN(1.)

        ! Getting index of the coordinate of the 4 vertices
        M(1) = M1(panel_idx)
        M(2) = M2(panel_idx)
        M(3) = M3 (panel_idx)
        M(4) = M4 (panel_idx)
        DO I = 1, 4

            tmp = sqrt(X(M(I))**2 + Y(M(I))**2)

            IF(I == 1 .or. tmp < rmin) THEN
                rmin = tmp
            END IF

            IF(I == 1 .or. tmp > rmax) THEN
                rmax = tmp
            END IF
        END DO

        r = sqrt(xx**2 + yy**2)
        theta = atan(yy/xx)
        v = theta/PI
        u = (r -rmin/2.0 -rmax/2.0)/(-rmin/2.0 + rmax/2.0)
    END SUBROUTINE transform_cartesian_to_uv


    SUBROUTINE transform_uv_to_cartesian(panel_idx, u ,v, xx, yy, zz)
        ! This subroutine converts a parametric coordinate to a cartesian coordinate
        ! It performs it using equations 2.49 and 2.46 (with alpha=4) of reference R4
        ! The generated u,v parametric coordinate are in the range -1, 1
        !  Arguments
        ! panel_idx (input):  The id of the panel where the coordinate belong to
        ! xx, yy, zz (input): The cartesian coordinates
        !  u,v (output): The parametric coordinate
        REAL :: u,v,xx,yy,zz, r, theta, rmin, rmax, tmp, N1, N2,   N3, N4
        INTEGER:: panel_idx, M(4), I
        REAL, PARAMETER :: PI=4.*ATAN(1.)

        ! Getting index of the coordinate of the 4 vertices
        M(1) = M1(panel_idx)
        M(2) = M2(panel_idx)
        M(3) = M3 (panel_idx)
        M(4) = M4 (panel_idx)

        N1 = (1-u)*(1-v)
        N2 = (1+u)*(1-v)
        N3 = (1+u)*(1+v)
        N4 = (1-u)*(1+v)


        xx = 0.25*N1*X(M(1)) + 0.25*N2*X(M(2)) + 0.25*N3*X(M(3)) + 0.25*N4*X(M(4))
        yy = 0.25*N1*Y(M(1)) + 0.25*N2*Y(M(2)) + 0.25*N3*Y(M(3)) + 0.25*N4*Y(M(4))
        zz=  0.25*N1*Z(M(1)) + 0.25*N2*Z(M(2)) + 0.25*N3*Z(M(3)) + 0.25*N4*Z(M(4))
    END SUBROUTINE transform_uv_to_cartesian



    REAL FUNCTION compute_jacobian(panel_idx, u, v)
        ! This function computes the Jacobian of the tranformation performed by the subroutine transform_uv_to_cartesian
        ! It performs it using equations 2.53 of reference R4
        !  Arguments
        ! panel_idx (input):  The id of the panel where the coordinate belong to
        !  u,v (input): The parametric coordinate in the range -1, 1
        REAL :: dxdu, dxdv, dydu, dydv, dzdu, dzdv ! Derivative of x,y,z with respect to u and so on
        REAL :: u,v,xx,yy,zz, r, theta, rmin, rmax, dN1du, dN1dv, dN2du, dN2dv, dN3du, dN3dv, dN4du, dN4dv
        REAL :: jx , jy, jz
        INTEGER:: panel_idx, M(4), I
        REAL, PARAMETER :: PI=4.*ATAN(1.)

        ! Getting index of the coordinate of the 4 vertices
        M(1) = M1(panel_idx)
        M(2) = M2(panel_idx)
        M(3) = M3 (panel_idx)
        M(4) = M4 (panel_idx)


        dN1du = -(1-v)
        dN1dv = -(1-u)
        dN2du = 1-v
        dN2dv = -(1+u)
        dN3du = (1+v)
        dN3dv = (1+u)
        dN4du = -(1+v)
        dN4dv = (1-u)

        dxdu = 0.25*dN1du*X(M(1)) + 0.25*dN2du*X(M(2)) + 0.25*dN3du*X(M(3)) + 0.25*dN4du*X(M(4))
        dxdv = 0.25*dN1dv*X(M(1)) + 0.25*dN2dv*X(M(2)) + 0.25*dN3dv*X(M(3)) + 0.25*dN4dv*X(M(4))

        dydu = 0.25*dN1du*Y(M(1)) + 0.25*dN2du*Y(M(2)) + 0.25*dN3du*Y(M(3)) + 0.25*dN4du*Y(M(4))
        dydv = 0.25*dN1dv*Y(M(1)) + 0.25*dN2dv*Y(M(2)) + 0.25*dN3dv*Y(M(3)) + 0.25*dN4dv*Y(M(4))

        dzdu = 0.25*dN1du*Z(M(1)) + 0.25*dN2du*Z(M(2)) + 0.25*dN3du*Z(M(3)) + 0.25*dN4du*Z(M(4))
        dzdv = 0.25*dN1dv*Z(M(1)) + 0.25*dN2dv*Z(M(2)) + 0.25*dN3dv*Z(M(3)) + 0.25*dN4dv*Z(M(4))

        jx = dydu*dzdv-dydv*dzdu
        jy = dxdv*dzdu-dxdu*dzdv
        jz=  dxdu*dydv-dxdv*dydu

        compute_jacobian = sqrt( jx**2 + jy**2 + jz**2)
    END FUNCTION compute_jacobian



    SUBROUTINE compute_green_derivative_wave_part(param, u, v, XPO, YPO, ZPO, panel_idx, green, dgreen)
        ! This subroutine computes the value of of wave part of the green function and it's derivative
        ! for the finite and infinite depth case
        ! Arguments
        ! param (input): common variables that contains useful variables
        !        See param definition for explanation of the member used here
        ! u, v (inputs): coordinate in parametric space of the field point
        ! XPO, YPO, ZPO (inputs): coordinate in cartesian of the source point
        ! panel_idx (input):  The field point panel
        ! green (output) The value of the wave part of the green function
        ! dgreen (output) The value of the derivative of the wave part of the green function
        complex :: green, dgreen
        REAL:: u,v,uf,vf
        INTEGER:: panel_idx
        INTEGER:: J
        REAL:: XPO,YPO,ZPO, xx, yy, zz
        TYPE(ParamsCommonInf) :: param
        REAL:: FS1(NFA,2),FS2(NFA,2)
        REAL:: VSX1(NFA,2),VSY1(NFA,2),VSZ1(NFA,2)
        REAL:: VSX2(NFA,2),VSY2(NFA,2),VSZ2(NFA,2)

        J = panel_idx

        CALL transform_uv_to_cartesian(panel_idx, u ,v, xx, yy, zz)

        IF((Depth .EQ. 0.) .OR. (param%AMH .GE. 20))  THEN ! infinite case
            CALL VVV(2, 2,J,XPO,YPO,ZPO, xx, yy, zz, param%SolverVar) ! finite part
            CALL VNV(2, J,XPO,YPO,ZPO, xx, yy, zz, param%SolverVar, green, dgreen)

        ELSE
            CALL VVV(2, 2,J,XPO,YPO,ZPO, xx, yy, zz, param%SolverVar)
            CALL VNVF(2, param%AM0,param%AMH,param%NEXP,J,XPO,YPO,ZPO, xx, yy, zz, param%SolverVar, green , dgreen)

        END IF
    END SUBROUTINE  compute_green_derivative_wave_part


    COMPLEX FUNCTION compute_green_wave_part(param, u, v, XPO, YPO, ZPO, panel_idx)
        ! This function computes the value of wave part of the green function
        ! for the finite and infinite depth case
        ! Arguments
        ! param (input): common variables that contains useful variables
        !        See param definition for explanation of the member used here
        ! u, v (inputs): coordinate in parametric space of the field point
        ! XPO, YPO, ZPO (inputs): coordinate in cartesian of the source point
        ! panel_idx (input):  The field point panel
        complex :: green, dgreen
        REAL:: u,v,uf,vf, XPO,YPO,ZPO
        INTEGER :: panel_idx
        TYPE(ParamsCommonInf) :: param


        CALL compute_green_derivative_wave_part(param, u, v, XPO, YPO, ZPO, panel_idx, green, dgreen)

        compute_green_wave_part = green

    END FUNCTION compute_green_wave_part


    COMPLEX FUNCTION compute_dgreen_wave_part(param, u, v, XPO, YPO, ZPO, panel_idx)
        ! This function computes the value of the derivative of the wave part of the green function
        ! for the finite and infinite depth case
        ! Arguments
        ! param (input): common variables that contains useful variables
        !        See param definition for explanation of the member used here
        ! u, v (inputs): coordinate in parametric space of the field point
        ! XPO, YPO, ZPO (inputs): coordinate in cartesian of the source point
        ! panel_idx (input):  The field point panel
        complex :: green, dgreen
        REAL:: u, v, uf, vf, XPO, YPO, ZPO
        INTEGER :: panel_idx
        TYPE(ParamsCommonInf) :: param

        CALL compute_green_derivative_wave_part(param, u, v, XPO, YPO, ZPO, panel_idx, green, dgreen)

        compute_dgreen_wave_part = dgreen

    END FUNCTION compute_dgreen_wave_part


    COMPLEX FUNCTION compute_dgreen_rankine_part(param, u, v, XPO, YPO, ZPO, panel_idx)
        ! This function computes the value of the derivative of the rankine part of the green function
        ! for the finite and infinite depth case
        ! Arguments
        ! param (input): common variables that contains useful variables
        !        See param definition for explanation of the member used here
        ! u, v (inputs): coordinate in parametric space of the field point
        ! XPO, YPO, ZPO (inputs): coordinate in cartesian of the source point
        ! panel_idx (input):  The field point panel
        complex :: green, dgreen
        REAL:: u,v,uf,vf
        REAL:: XPO,YPO,ZPO, xx, yy, zz, tmp
        INTEGER:: J
        INTEGER :: panel_idx
        TYPE(ParamsCommonInf) :: param

        J = panel_idx

        CALL transform_uv_to_cartesian(panel_idx, u ,v, xx, yy, zz)

        IF((Depth .EQ. 0.) .OR. (param%AMH .GE. 20))  THEN ! infinite case
            CALL VVV(2, 1,J,XPO,YPO,ZPO, xx, yy, zz, param%SolverVar) ! finite part

        ELSE
            CALL VVV(2, 2,J,XPO,YPO,ZPO, xx, yy, zz, param%SolverVar)

        END IF



        tmp = XN(J)*param%SolverVar%VSXP + YN(J)*param%SolverVar%VSYP + ZN(J)*param%SolverVar%VSZP

        compute_dgreen_rankine_part = cmplx(tmp, 0)

    END FUNCTION compute_dgreen_rankine_part


    COMPLEX FUNCTION compute_green_rankine_part(param, u, v, XPO, YPO, ZPO, panel_idx)
        ! This function computes the value of the rankine part of the green function
        ! for the finite and infinite depth case
        ! Arguments
        ! param (input): common variables that contains useful variables
        !        See param definition for explanation of the member used here
        ! u, v (inputs): coordinate in parametric space of the field point
        ! XPO, YPO, ZPO (inputs): coordinate in cartesian of the source point
        ! panel_idx (input):  The field point panel
        complex :: green, dgreen
        REAL:: u,v,uf,vf
        REAL:: XPO,YPO,ZPO, xx, yy, zz, tmp
        INTEGER:: J
        INTEGER :: panel_idx
        TYPE(ParamsCommonInf) :: param

        J = panel_idx

        CALL transform_uv_to_cartesian(panel_idx, u ,v, xx, yy, zz)

        IF((Depth .EQ. 0.) .OR. (param%AMH .GE. 20))  THEN ! infinite case
            CALL VVV(2, 1,J,XPO,YPO,ZPO, xx, yy, zz, param%SolverVar) ! finite part

        ELSE
            CALL VVV(2, 2,J,XPO,YPO,ZPO, xx, yy, zz, param%SolverVar)

        END IF

        compute_green_rankine_part = cmplx(param%SolverVar%FSP, 0)

    END FUNCTION compute_green_rankine_part


    COMPLEX FUNCTION COMPUTE_NORMAL_POTENTIAL(param, panel_idx, xx, yy, zz)
        ! Function to compute the normal potential at any point (xx, yy, zz) of a panel identified by panel_idx
        ! It is expected to replace the old NVEL variable for higher order panel method
        ! The old NVEL only contain the normal potential at the centroid of the panel

        REAL :: xx,yy,zz ! The coordinate of the point
        INTEGER :: panel_idx ! The panel id
        INTEGER :: i
        COMPLEX :: Phi,p,Vx,Vy,Vz
        TYPE(ParamsCommonInf) :: param ! The variable containing useful parameters
        REAL :: direction(3), axis(3)

        direction(1) = param%rad_case(param%rad_case_idx, 3)
        direction(2) = param%rad_case(param%rad_case_idx, 4)
        direction(3) = param%rad_case(param%rad_case_idx, 5)

        axis(1) = param%rad_case(param%rad_case_idx, 6)
        axis(2) = param%rad_case(param%rad_case_idx, 7)
        axis(3) = param%rad_case(param%rad_case_idx, 8)

        i = panel_idx

        IF (param%is_wave == 1) THEN
            CALL Compute_Wave(param%kwave,param%w,param%beta,xx,yy,zz,Phi,p,Vx,Vy,Vz,param%Environment)
            COMPUTE_NORMAL_POTENTIAL = Vx*param%Mesh%N(1,i)+Vy*param%Mesh%N(2,i)+Vz*param%Mesh%N(3,i)
        ELSE
            IF (int(param%rad_case(param%rad_case_idx, 2)) == 2 ) THEN
                IF(param%Mesh%cPanel(i) == int(param%rad_case(param%rad_case_idx, 1))) THEN

                    Vx = direction(2)*(zz - axis(3)) - direction(3)*(yy - axis(2))
                    Vy = direction(3)*(xx - axis(1)) - direction(1)*(zz - axis(3))
                    Vz = direction(1)*(yy - axis(2)) - direction(2)*(xx - axis(1))
                    COMPUTE_NORMAL_POTENTIAL = Vx*param%Mesh%N(1,i)+Vy*param%Mesh%N(2,i)+Vz*param%Mesh%N(3,i)

                ELSE
                    COMPUTE_NORMAL_POTENTIAL = param%NVEL(i)  ! cmplx(0, 0)
                END IF
            ELSE
                COMPUTE_NORMAL_POTENTIAL = param%NVEL(i)
            END IF
        END IF

    END FUNCTION COMPUTE_NORMAL_POTENTIAL




    SUBROUTINE SOLVE_POTENTIAL_HIGH(param, PRESSURE, PHI, MeshFS, Switch_FS,  Switch_Potential, NDS, Force, &
        & out_potential)
        !This subroutine solves a boundary element method for a finite or infinite depth using
        !   a three dimensional higher order method. The higher order method used is based on
        !   B-Spline.
        !   The influence matrix is setup using equation 15.30 of reference R2
        !   The corresponding linear system is solved using an SVD (Singular Value Decomposition) solver.
        !   More specific description can be read in section 15.5 of reference R2.
        !   Here the green function and derivative are computed using a series approximation and a tabulation
        !   as previously done for the lower order panel method in Nemoh.
        !   After solving the influence matrix, the potential at the boundary and at the given points in the fluid
        !   domain is computed. The pressure at the centroid of each patch are computed and the forces are derived.
        !  Arguments
        !  param (input/temporary vars): common variables that contains useful variables
        !        See param definition for explanation of the member used here
        !  PRESSURE (output) The array to contain the pressure at the centroid of each patch
        !  PHI (output) The array to contain the potential in the fluid domain of the given point
        !  MeshFS (input) Mesh structure used to get the points for which to compute the potential in the fluid
        !                 domain
        !  Switch_FS (input): Whether or not to compute the potential at the fluid domain
        !  Switch_Potential (input) : Whether or not to save the pressure
        !  NDS (input):  Array containing the nds at each center of a patch
        !  Force (output): Array to contain the computed dynamic forces
        !  out_potential (output): Array to contain the computed potential/pressure/centroid values

        COMPLEX, DIMENSION(:, :), ALLOCATABLE :: ZIJ
        COMPLEX, DIMENSION(:, :), ALLOCATABLE :: A
        COMPLEX, DIMENSION(:, :), ALLOCATABLE :: B
        complex :: tmp, tmp2
        COMPLEX, DIMENSION(:), POINTER :: ZPB,ZPS
        INTEGER :: I, J,  K,  II, JJ, KK, num, curx, cury
        TYPE(TempVar), POINTER :: SolverVar
        REAL AMH
        REAL :: u,v, tmpr
        COMPLEX,DIMENSION(:),ALLOCATABLE :: NVEL
        INTEGER :: NEXP
        TYPE(TMesh) :: MeshFS
        COMPLEX :: PRESSURE(:)
        COMPLEX :: PHI(:)
        INTEGER :: Switch_FS, Switch_Potential
        REAL,DIMENSION(:,:) :: NDS
        COMPLEX,DIMENSION(:) :: Force
        COMPLEX,PARAMETER :: CII=CMPLX(0.,1.)
        REAL :: out_potential(:)
        REAL, PARAMETER :: PI=4.*ATAN(1.)
        INTEGER, DIMENSION(:), ALLOCATABLE:: IPIV
        REAL, DIMENSION(:), ALLOCATABLE:: S, RWORK
        INTEGER:: INFO, RANK, LWORK
        COMPLEX, DIMENSION(:), ALLOCATABLE:: WORK



        TYPE(ParamsCommonInf), TARGET :: param

        SolverVar => param%SolverVar

        CALL COMPUTE_NEXP(param)

        ZPB => SolverVar%ZPB
        ZPS => SolverVar%ZPS

        CALL INITIALIZE_BEM_HIGH(param)

        LWORK = 2*min(SolverVar%Num_Unknows,SolverVar%Num_Unknows) + max(SolverVar%Num_Unknows, &
            & SolverVar%Num_Unknows,1) + 5

        ALLOCATE(WORK(LWORK))
        ALLOCATE(RWORK(5*min(SolverVar%Num_Unknows,  SolverVar%Num_Unknows)))
        ALLOCATE(ZIJ(SolverVar%Num_Unknows, SolverVar%Num_Unknows +1))
        ALLOCATE(IPIV(SolverVar%Num_Unknows))
        ALLOCATE(S(SolverVar%Num_Unknows))
        ALLOCATE(A(SolverVar%Num_Unknows, SolverVar%Num_Unknows))
        ALLOCATE(B(SolverVar%Num_Unknows, 1))

        ZIJ = cmplx(0., 0.)

        num = 0

        param%SolverVar = SolverVar

        param%use_cartesian_coordinate = 0
        param%compute_potential = 0



        ! Computing influence Matrix
        DO J=1,SolverVar%Num_Patches


            DO I=1,SolverVar%Num_Basis_U(J)


                DO K=1,SolverVar%Num_Basis_V(J)

                    curx = SolverVar%Panel_Idx(J) + I*K
                    param%j = J
                    param%i = I
                    param%k = K


                    ! dik Computing integral at equation 15.32 of reference R2
                    DO II = 1,SolverVar%Num_Basis_U(J)

                        DO KK=1,SolverVar%Num_Basis_V(J)

                            param%ii = II
                            param%jj = J
                            param%kk = KK

                            ! Using tmp to avoid long formula
                            tmp = qgss2d_simple(influence_dik, -1., 1., -1., 1., 3, param)
                            ZIJ(curx, SolverVar%Panel_Idx(J) + II*KK) = 2*PI*tmp

                        END DO

                    END DO




                    !! DIK Computing integral at equation 15.33 of reference R2
                    DO JJ=1,SolverVar%Num_Patches

                        DO II = 1,SolverVar%Num_Basis_U(JJ)

                            DO KK=1,SolverVar%Num_Basis_V(JJ)



                                cury = SolverVar%Panel_Idx(JJ) + II*KK
                                param%ii = II
                                param%jj = JJ
                                param%kk = KK

                                ! Using tmp to avoid long formulas
                                tmp = qgss2d_simple(influence_bigdik, -1., 1., -1., 1., 3, param)
                                ZIJ(curx, cury) = ZIJ(curx, cury) + tmp

                            END DO

                        END DO

                    END DO



                    !!  SI Computing integral at equation 15.34 of reference R2
                    cury = SolverVar%Num_Unknows + 1

                    param%jj = J

                    ZIJ(curx, cury) = qgss2d_simple(influence_si, -1., 1., -1., 1., 3, param)



                END DO




            END DO



        END DO

        ! Setting up matrix A and B of the linear equations Ax=B
        DO I = 1, SolverVar%Num_Unknows

            DO J=1, SolverVar%Num_Unknows

                A(I, J) = ZIJ(I,J)

            END DO

            B(I, 1) = ZIJ(I, SolverVar%Num_Unknows+1)

        END DO
        ! Computing influence matrix ended


        ! Solve the Matrix

        !------------------------------------------------!
        ! Solve this using a SVD decomposition, solving equation 15.30 of reference R2
        CALL CGELSS( SolverVar%Num_Unknows, SolverVar%Num_Unknows, 1, A, SolverVar%Num_Unknows, B, &
            & SolverVar%Num_Unknows, S, 1e-20, RANK, WORK, LWORK, RWORK, INFO )

        IF( INFO.GT.0 ) THEN
            WRITE(*,*)'The algorithm for computing the SVD failed to converge'
            WRITE(*,*)INFO, ' off-diagonal elements of an intermediate'
            WRITE(*,*)'bidiagonal form did not converge to zero.'
            STOP
        END IF

        !------------------------------------------------!




        ! Computing the boundary potential at the centroid of a panel/patch
        ! We can re-use the B-Spline formulation of the potential but we are using the general
        ! equation 15.29 of reference R2 to get it.
        ! In the higher order panel method we can get the potential at any point of a patch
        ! We are returning here only the potential at the centroid of a patch
        DO J=1,SolverVar%Num_Patches

            ! Commented, once the subroutine transform_cartesian_to_uv is not anymore mocked it could be use
            ! if needed
            !CALL transform_cartesian_to_uv(J, XG(J), YG(J), ZG(J), u, v) to re-enable here
            !ZPB(J) = COMPUTE_POTENTIAL_BOUNDARY_PATCH(J, u, v, B, param)
            ZPB(J) = COMPUTE_POTENTIAL_DOMAIN_HIGH(XG(J), YG(J), ZG(J), param, B)
            ZPS(J) = cmplx(0., 0.)

        END DO

        !       Assemble pressure at patches centroid
        DO J=1,SolverVar%Num_Patches
            PRESSURE(J)=RHO*CII*param%W*ZPB(J) !*AIRE(i)
        END DO

         !       Save free surface elevation
        IF (Switch_FS.EQ.1) THEN
            DO j=1,MeshFS%Npoints
                PHI(j) = COMPUTE_POTENTIAL_DOMAIN_HIGH(MeshFS%X(1,j),MeshFS%X(2,j),0.,param, B)
            END DO
        END IF

        !       Save output
        IF (Switch_Potential.EQ.1) THEN
            CALL SAVE_POTENTIAL(PRESSURE, out_potential)
        END IF

        !       Calculate force coefficients at centroid of panel
        DO i=1,param%n_integration
            DO j=1,param%NpanIsy
                Force(i)=Force(i)-PRESSURE(j)*NDS(i,j)
            END DO
        END DO


        ! Deallocate allocated arrays
        if (allocated(ZIJ)) deallocate(ZIJ)
        if (allocated(IPIV)) deallocate(IPIV)
        if (allocated(A)) deallocate(A)
        if (allocated(B)) deallocate(B)
        if (allocated(RWORK)) deallocate(RWORK)
        if (allocated(WORK)) deallocate(WORK)
        CALL DEALLOCATE_BEM_HIGH(param)

    END SUBROUTINE SOLVE_POTENTIAL_HIGH

    COMPLEX FUNCTION COMPUTE_POTENTIAL_BOUNDARY_PATCH(panel_idx, u, v, B, param)
        ! Compute potential at a given point in the boundary for a patch using equation 7.4 of reference R2
        ! It is only valid for the higher panel method
        ! Arguments (inputs)
        !panel_idx: The id of the panel of the given point
        ! u, v: The parametric coordinate
        ! B: Array, solution of the linear system. Contain the phi coefficients
        ! param: common variables that contains useful variables
        !        See param definition for explanation of the member used here
        INTEGER :: panel_idx, J,   I,   K ! The id of the patch/panels
        REAL :: xx, yy, zz ! The coordinate of a point in the domain of the given patch
        REAL:: tmpr
        INTEGER:: curx
        COMPLEX ::   tmp
        REAL :: u,v
        TYPE(ParamsCommonInf) :: param
        COMPLEX, DIMENSION(:, :), ALLOCATABLE :: B

        J = panel_idx

        tmp = cmplx(0., 0.)

        DO I=1, param%SolverVar%Num_Basis_U(J)

            DO K=1, param%SolverVar%Num_Basis_V(J)

                curx = param%SolverVar%Panel_Idx(J) + I*K

                tmpr = GET_BSPLINE(u, param%SolverVar%KNOTS_U(J), I)*GET_BSPLINE(v, param%SolverVar%KNOTS_U(J), K)

                tmp = tmp + B(curx, 1)*tmpr

            END DO


        END DO

        COMPUTE_POTENTIAL_BOUNDARY_PATCH = tmp

    END FUNCTION COMPUTE_POTENTIAL_BOUNDARY_PATCH


    COMPLEX FUNCTION COMPUTE_POTENTIAL_DOMAIN_HIGH(xx, yy, zz, param, B)
        ! Compute the potential at a given point in the fluid domain or in the boundary
        ! for the higher panel method
        ! Arguments
        ! xx,yy,zz: Cartesian coordinate of the given point
        ! B: Array, solution of the linear system. Contain the phi coefficients
        ! param: common variables that contains useful variables
        !        See param definition for explanation of the member used here
        TYPE(ParamsCommonInf) :: param
        COMPLEX :: tmp
        REAL :: xx, yy, zz ! The coordinate of a point in the domain of the given patch
        REAL :: uf,vf
        INTEGER:: J
        COMPLEX, DIMENSION(:, :), ALLOCATABLE :: B ! The solution to the linear system

        COMPUTE_POTENTIAL_DOMAIN_HIGH = cmplx(0, 0)

        param%use_cartesian_coordinate = 1
        param%compute_potential = 1
        param%B = B

        param%xx = xx
        param%yy = yy
        param%zz = zz

        DO J=1,param%SolverVar%Num_Patches

            param%jj = J

            tmp = qgss2d_simple(influence_ugreen, -1., 1., -1., 1., 3, param)

            tmp = tmp - qgss2d_simple(influence_udgreen, -1., 1., -1., 1., 3, param)

            COMPUTE_POTENTIAL_DOMAIN_HIGH = COMPUTE_POTENTIAL_DOMAIN_HIGH + tmp
        END DO

    END FUNCTION COMPUTE_POTENTIAL_DOMAIN_HIGH

    SUBROUTINE INITIALIZE_BEM_HIGH(param)
        ! This initializes the higher order panel method. It setups the common parameters.
        ! It computes the knot vectors, their length, the number of basis per patch
        ! param(input/output): common variables that contains useful variables
        !        See param definition for explanation of the member used here
        TYPE(ParamsCommonInf), TARGET :: param
        TYPE(TempVar), POINTER :: SolverVar
        INTEGER :: I,J, num
        REAL :: total

        SolverVar => param%SolverVar


        SolverVar%Num_Patches = IMX

        ALLOCATE(SolverVar%Num_Basis_U(IMX))
        ALLOCATE(SolverVar%Num_Basis_V(IMX))
        ALLOCATE(SolverVar%Panel_Idx(IMX))
        ALLOCATE(SolverVar%KNOTS_U(IMX))
        ALLOCATE(SolverVar%KNOTS_V(IMX))

        DO I = 1, SolverVar%Num_Patches
            SolverVar%Num_Basis_U(I) = param%NSPLIN + param%KSPLIN -1
            SolverVar%Num_Basis_V(I) = param%NSPLIN + param%KSPLIN -1
        END DO

        num = 0

        DO I = 1, SolverVar%Num_Patches

            SolverVar%Panel_Idx(I) = num

            num = num + SolverVar%Num_Basis_U(I)*SolverVar%Num_Basis_V(I)

        END DO

        SolverVar%Num_Unknows = num

        DO I = 1, SolverVar%Num_Patches

            SolverVar%KNOTS_U(I)%POLY_ORDER = param%KSPLIN
            SolverVar%KNOTS_V(I)%POLY_ORDER = param%KSPLIN

            SolverVar%KNOTS_V(I)%LENGTH = SolverVar%Num_Basis_V(I) + param%KSPLIN
            SolverVar%KNOTS_U(I)%LENGTH = SolverVar%Num_Basis_U(I) + param%KSPLIN

            ALLOCATE(SolverVar%KNOTS_V(I)%KNOTS(0:SolverVar%KNOTS_V(I)%LENGTH))
            ALLOCATE(SolverVar%KNOTS_U(I)%KNOTS(0:SolverVar%KNOTS_U(I)%LENGTH))

        END DO


        ! Computing knot vector as suggested in R1 page 195
        DO I = 1, SolverVar%Num_Patches

            total = SolverVar%KNOTS_V(I)%LENGTH + 1 -2*SolverVar%KNOTS_V(I)%POLY_ORDER

            DO J=0, SolverVar%KNOTS_V(I)%LENGTH

                IF (J < SolverVar%KNOTS_V(I)%POLY_ORDER) THEN
                    SolverVar%KNOTS_V(I)%KNOTS(J) = -1
                ELSE IF (SolverVar%KNOTS_V(I)%LENGTH-J < SolverVar%KNOTS_V(I)%POLY_ORDER) THEN
                    SolverVar%KNOTS_V(I)%KNOTS(J) = 1
                ELSE
                    SolverVar%KNOTS_V(I)%KNOTS(J) = -1+ 2.*real(J-SolverVar%KNOTS_V(I)%POLY_ORDER)/total

                END IF

            END DO

            total = SolverVar%KNOTS_U(I)%LENGTH + 1 -2*SolverVar%KNOTS_U(I)%POLY_ORDER

            DO J=0, SolverVar%KNOTS_U(I)%LENGTH

                IF (J < SolverVar%KNOTS_U(I)%POLY_ORDER) THEN
                    SolverVar%KNOTS_U(I)%KNOTS(J) = -1
                ELSE IF (SolverVar%KNOTS_U(I)%LENGTH-J < SolverVar%KNOTS_U(I)%POLY_ORDER) THEN
                    SolverVar%KNOTS_U(I)%KNOTS(J) = 1
                ELSE
                    SolverVar%KNOTS_U(I)%KNOTS(J) = -1+ 2.*real(J-SolverVar%KNOTS_U(I)%POLY_ORDER)/total

                END IF

            END DO


        END DO
    END SUBROUTINE INITIALIZE_BEM_HIGH


    SUBROUTINE DEALLOCATE_BEM_HIGH(param)
        ! This deallocates allocated data done in INITIALIZE_BEM_HIGH
        TYPE(ParamsCommonInf) :: param
        TYPE(TempVar), TARGET :: SolverVar
        INTEGER :: I,J

        SolverVar = param%SolverVar

        DEALLOCATE(SolverVar%Num_Basis_U)
        DEALLOCATE(SolverVar%Num_Basis_V)
        DEALLOCATE(SolverVar%Panel_Idx)

        DO I = 1, SolverVar%Num_Patches
            DEALLOCATE(SolverVar%KNOTS_V(I)%KNOTS)
            DEALLOCATE(SolverVar%KNOTS_U(I)%KNOTS)

        END DO

        DEALLOCATE(SolverVar%KNOTS_U)
        DEALLOCATE(SolverVar%KNOTS_V)

    END SUBROUTINE DEALLOCATE_BEM_HIGH


end MODULE SOLVE_BEM_HIGH
