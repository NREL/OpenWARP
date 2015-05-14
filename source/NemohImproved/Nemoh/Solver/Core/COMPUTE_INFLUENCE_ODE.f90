!--------------------------------------------------------------------------------------
!
!Copyright (C) 2014 TopCoder Inc., All Rights Reserved.
!
!--------------------------------------------------------------------------------------

!   This module computes the influence coefficient according to R2 for the regular part and
!  to R1 for the rankine part. It only perform the computation for the infinite depth case
!
!  R1 Distributions of sources and normal dipoles
!over a quadrilateral panel (Available in another contest (Dipoles, request permission) document paper 2
!http://community.topcoder.com/tc?
!module=DownloadDocument&docid=27516374 )
!
!  R2  A second order Ordinary Differential Equation for the frequency
!domain Green function http://www.iwwwfb.org/Abstracts/iwwwfb28/iwwwfb28_12.pdf
!
! Contest Code Acceleration of the Calculation of Influence Coefficients of Nemoh
!
!   @author yedtoss
!   @version 1.0

module COMPUTE_INFLUENCE_ODE

    USE ODE
    USE COM_VAR
    USE COMMON_TYPE
    implicit none


contains

    complex function influence_infinite_source_wave_part(w,y, param)
        ! This is the function representing the expression of the ODE of the green function
        ! See equation 11 of reference R2
        real r, w, Z
        complex y(0:10)
        TYPE(ParamsCommonInf) :: param
        complex tmp

        r = param%r
        Z = param%Z

        tmp = 2.*(1. + Z*w**2)/(sqrt(r**2 + Z**2)) + w*(w**2*Z + 3./4.)*y(1)

        tmp = tmp - (w**4*(r**2 + Z**2) + w**2*Z +1 )* y(0)
        if(abs(w) < tolerance) then
            tmp = 0
        else
            tmp = tmp/(w**2/4.)
        end if

        influence_infinite_source_wave_part = tmp
    end function


    complex function influence_infinite_dipoles_wave_part(w,y, param)
        ! This is the function representing the expression of the ODE of the derivative of the green function with
        ! respect to r
        ! See corresponding equation in Deployment guide

        real r, w, Z
        complex y(0:10)
        TYPE(ParamsCommonInf) :: param
        complex tmp

        r = param%r
        Z = param%Z

        tmp = -6.*r*(1. + Z*w**2)/((r**2 + Z**2)**(3./2.))
        tmp = tmp + w*(w**2*Z + 3./4.)*y(1)

        tmp = tmp - (w**4*(r**2 + Z**2) + 3.*w**2*Z + 3 )* y(0)

        if(abs(w) < tolerance) then
            tmp = 0
        else
            tmp = tmp/(w**2/4.)
        end if

        influence_infinite_dipoles_wave_part = tmp

    end function


    complex function influence_infinite_dipoles_wave_part_z(w,y, param)
        ! This is the function representing the expression of the ODE of the derivative of the green function with
        ! respect to z
        ! See corresponding equation in Deployment guide

        real r, w, z
        complex y(0:10)
        TYPE(ParamsCommonInf) :: param
        complex tmp

        r = param%r
        Z = param%Z

        tmp = (-8.*Z -(w**2)*(4.*Z**2 - 2.*r**2))/((r**2 + Z**2)**(3./2.))
        tmp = tmp + w*(w**2*Z + 3./4.)*y(1)

        tmp = tmp - (w**4*(r**2 + Z**2) + 3.*w**2*Z + 4 )* y(0)

        if(abs(w) < tolerance) then
            tmp = 0
        else
            tmp = tmp/(w**2/4.)
        end if

        influence_infinite_dipoles_wave_part_z = tmp

    end function

    real function compute_rankine_dipoles(ISP, IFP)

        ! Computing formula 2.14 with 2.15 as insight of
        ! Distributions of sources and normal dipoles over a quadrilateral panel

        INTEGER :: ISP ! source point

        INTEGER :: IFP !Field point

        INTEGER:: I, M(4)

        REAL:: s1, c1, s2, c2, nu_n, nu_n1, psi_n, psi_n1, delta_n, delta_n1, rn, rn_1, xx, yy, zz

        REAL ::  s3, c3

        REAL, PARAMETER :: PI=4.*ATAN(1.)

        !!$OMP CRITICAL
        !We are safe without a critical or atomic operation. Although we are
        ! using it's value, a value of 0 is safe and when it is 1 it is also safe
        IF(is_rankine_dipoles_computed == 1) THEN
            compute_rankine_dipoles = rankine_dipoles_cache(ISP, IFP)
            RETURN
        END IF
        !!$OMP END CRITICAL


        ! Getting index of the coordinate of the 4 vertices
        M(1) = M1(IFP)
        M(2) = M2(IFP)
        M(3) = M3 (IFP)
        M(4) = M4 (IFP)

        compute_rankine_dipoles = 0

        ! Getting centroid coordinate
        xx = XG(ISP) ! Source point x
        yy = YG(ISP) ! Source point y
        zz = ZG(ISP) ! Source point z


        DO I =1,4

            nu_n = Y(M(I)) ! y coordinate of vertice n
            nu_n1 = Y(modulo(M(I),4) +1) ! y coordinate of vertice n +1
            psi_n = X(M(I)) ! x coordinate of vertice n
            psi_n1 = X(modulo(M(I),4) +1) ! x coordinate of vertice n +1

            delta_n = nu_n1 - nu_n
            delta_n1 = Y(modulo(M(I+1),4) + 1) - nu_n1

            ! If both vertices are the same, it means we have a triangle instead of quadrilateral.
            ! We can thus skip this point
            if (abs(delta_n) < tolerance .and. abs(psi_n1 - psi_n) < tolerance) then
                CYCLE
            end if


            ! Spherical coordinate http://mathworld.wolfram.com/SphericalCoordinates.html
            ! Radial distance between vertice n and the field point
            rn = sqrt((xx - psi_n)**2 + (yy - nu_n)**2 + zz**2)
            ! Radial distance between vertice n+1 and the field point
            rn_1 = sqrt((xx - psi_n1)**2 + (yy - nu_n1)**2 + zz**2)

            s1 = delta_n* ((xx -psi_n)**2 + zz**2) - delta_n* (xx -psi_n) *(yy -nu_n)
            c1 = rn*zz* delta_n

            s2 = delta_n* ((xx -psi_n1)**2 + zz**2) - delta_n* (xx -psi_n1) *(yy -nu_n1)
            c2 = rn_1*zz* delta_n


            s3 = s1*c2 - s2*c1
            c3 = c1*c2 + s1*s2

            ! If c3 is near 0 the angle tends to pi/2
            if (abs(c3) < tolerance) then
                compute_rankine_dipoles = compute_rankine_dipoles + PI/2.
            else
                compute_rankine_dipoles = compute_rankine_dipoles + atan(s3/c3)
            end if

        END DO

        !!$OMP CRITICAL
        !We are safe without a critical or atomic operation. Although we are
        ! updating it's value, any thread updating it simulatenous will do so with the same value
        rankine_dipoles_cache(ISP, IFP) = compute_rankine_dipoles
        !!$OMP END CRITICAL


    end function


    real function compute_rankine_sources(ISP, IFP)


        ! Computing formula 3.10
        ! Distributions of sources and normal dipoles over a quadrilateral panel

        INTEGER :: ISP ! source point

        INTEGER :: IFP !Field point

        REAL:: s1, c1, s2, c2, nu_n, nu_n1, psi_n, psi_n1, delta_n, delta_n1, rn, rn_1, xx, yy, zz

        REAL ::  s3, c3,s ,c, theta_n, theta_n1, sn

        INTEGER:: I, M(4)

        REAL, PARAMETER :: PI=4.*ATAN(1.)

        !!$OMP CRITICAL
        !We are safe without a critical or atomic operation. Although we are
        ! using it's value, a value of 0 is safe and when it is 1 it is also safe
        IF(is_rankine_sources_computed == 1) THEN
            compute_rankine_sources = rankine_sources_cache(ISP, IFP)
            RETURN
        END IF
        !!$OMP END CRITICAL

        ! Getting index of the coordinate of the 4 vertices
        M(1) = M1(IFP)
        M(2) = M2(IFP)
        M(3) = M3 (IFP)
        M(4) = M4 (IFP)

        compute_rankine_sources = 0

        ! Getting centroid coordinate
        xx = XG(ISP) ! Source point x
        yy = YG(ISP) ! Source point y
        zz = ZG(ISP) ! Source point z


        DO I =1,4

            nu_n = Y(M(I)) ! y coordinate of vertice n
            nu_n1 = Y(modulo(M(I),4) +1) ! y coordinate of vertice n +1
            psi_n = X(M(I)) ! x coordinate of vertice n
            psi_n1 = X(modulo(M(I),4) +1) ! x coordinate of vertice n +1

            delta_n = nu_n1 - nu_n
            delta_n1 = Y(modulo(M(I+1),4) + 1) - nu_n1

            ! If both vertices are the same, it means we have a triangle instead of quadrilateral.
            ! We can thus skip this point
            if (abs(delta_n) < tolerance .and. abs(psi_n1 - psi_n) < tolerance) then
                CYCLE
            end if

            ! If delta_n is near 0 the angle tends to pi/2
            if (abs(delta_n) < tolerance) then
                theta_n = PI/2.
            else
                !http://stackoverflow.com/questions/2676719/calculating-the-angle-between-the-line-defined-by-two-points
                theta_n = atan(delta_n/(psi_n1 - psi_n))
            endif

            ! sn is the length between vertices n and vertices n+1. Note that z=0 because it is a flat panel
            sn = sqrt((psi_n1 - psi_n)**2 + (delta_n1- delta_n)**2)



            ! Spherical coordinate http://mathworld.wolfram.com/SphericalCoordinates.html
            ! Radial distance between vertice n and the field point
            rn = sqrt((xx - psi_n)**2 + (yy - nu_n)**2 + zz**2)
            ! Radial distance between vertice n+1 and the field point
            rn_1 = sqrt((xx - psi_n1)**2 + (yy - nu_n1)**2 + zz**2)

            ! Using intermediate variable s to avoid long expression or complicated multi line expressions
            s = (xx- psi_n)* sin(theta_n) - (yy - nu_n)*cos(theta_n)*log((rn + rn_1 + sn)/(rn + rn_1 - sn))

            compute_rankine_sources = compute_rankine_sources + s



        END DO

        compute_rankine_sources = compute_rankine_sources - zz*compute_rankine_dipoles(ISP, IFP)
        !!$OMP CRITICAL
        !We are safe without a critical or atomic operation. Although we are
        ! updating it's value, any thread updating it simulatenous will do so with the same value
        rankine_sources_cache(ISP, IFP) = compute_rankine_sources
        !!$OMP END CRITICAL
    end function

    subroutine compute_influence_infinite_sources(w, ZIJ)

        INTEGER :: ISP ! source point

        INTEGER :: IFP !Field point

        REAL:: r, G0, dG0, Z, w

        complex yi(0:10),t(50)
        COMPLEX, DIMENSION(:, :):: ZIJ ! influence coefficients
        COMPLEX:: Gn(IMX, IMX) ! Cache of green function

        TYPE(ParamsCommonInf) :: param

        DO ISP=1,IMX

            DO IFP=1,IMX

                Z = ZG(ISP) + ZG(IFP)
                r = sqrt((XG(ISP) - XG(IFP))**2 + (YG(ISP) - YG(IFP))**2 )

                G0 = cmplx(2./(sqrt(r**2 + Z**2)), 0.)
                dG0 = cmplx(0 , 0)

                yi(0) = G0
                yi(1) = dG0

                param%w = w
                param%Z = Z
                param%r = r


                !  Compute Wave part/ regular part of green function
                if(ISP <= IFP) then
                    ! A 3x4 runge kutta is enough for convergence
                    t = Equadifnc(influence_infinite_source_wave_part, 0., w, yi, 1, 2, 3, param )
                    Gn(ISP, IFP) = t(2)
                    Gn(IFP, ISP) = t(2)
                else
                    t(2) = Gn(IFP, ISP)
                end if

                !  Double integral done by multiplying by the area
                ZIJ(ISP, IFP) = t(2)* AIRE(IFP) +  compute_rankine_sources(ISP, IFP)

            END DO

        END DO
    end subroutine compute_influence_infinite_sources



    subroutine compute_influence_infinite_dipoles(w, ZIJ)
        ! This compute the influence coefficients containing the derivative of the green function
        !

        INTEGER :: ISP ! source point correspond to the collocation point

        INTEGER :: IFP !Field point correspond to the flat panel

        REAL:: r, G0, dG0, Z, w

        complex yi(0:10),t(50), yi1(0:10),t1(50), tmp
        COMPLEX, DIMENSION(:, :):: ZIJ ! influence coefficients
        COMPLEX, DIMENSION(IMX, IMX):: GR, GZ ! Cache of green function derivative with respect to r and z respectively

        TYPE(ParamsCommonInf) :: param

        DO ISP=1,IMX
            DO IFP=1,IMX

                Z = ZG(ISP) + ZG(IFP)
                r = sqrt((XG(ISP) - XG(IFP))**2 + (YG(ISP) - YG(IFP))**2 )
                if(abs(r) < tolerance .or. ISP == IFP) then
                    ZIJ(ISP, IFP) = 0.5
                    CYCLE
                end if

                G0 = cmplx(-2.*r/((r**2 + Z**2)**(3./2.)), 0.)
                dG0 = cmplx(0, 0)

                yi(0) = G0
                yi(1) = dG0

                param%w = w
                param%Z = Z
                param%r = r

                !  Compute Wave part/ regular part of green function derivative with respect to r
                if(ISP <= IFP) then
                    ! A 3x4 runge kutta is enough for convergence
                    t = Equadifnc(influence_infinite_dipoles_wave_part, 0., w, yi, 1, 2, 3, param )
                    GR(ISP, IFP) = t(2)
                    GR(IFP, ISP) = t(2)
                else
                    t(2) = GR(IFP, ISP)
                end if

                ! Compute Wave part/ regular part of green function derivative with respect to Z
                yi1(0) = cmplx((-2.*Z)/((r**2 + Z**2)**(3./2.)), 0.)
                yi1(1) = cmplx(0, 0)
                if (abs(ZN(IFP)) < tolerance) then
                    t1(2) = 0
                else
                    if(ISP <= IFP) then
                        ! A 3x4 runge kutta is enough for convergence
                        t1 = Equadifnc(influence_infinite_dipoles_wave_part_z, 0., w, yi1, 1, 2, 3, param )
                        GZ(ISP, IFP) = t1(2)
                        GZ(IFP, ISP) = t1(2)
                    else
                        t1(2) = GZ(IFP, ISP)
                    end if
                end if

                ! Compute coefficient transforming derivative in nM' to derivative in r and z
                tmp = cmplx(-XN(IFP)*(XG(ISP) - XG(IFP))/r -YN(IFP)*(YG(ISP) - YG(IFP))/r, 0)
                ! Double integral done by multiplying by the area
                ZIJ(ISP, IFP) = (t(2)*tmp + ZN(IFP)*t1(2))* AIRE(IFP) +  compute_rankine_dipoles(ISP, IFP)

            END DO
        END DO
    end subroutine compute_influence_infinite_dipoles
end module COMPUTE_INFLUENCE_ODE
