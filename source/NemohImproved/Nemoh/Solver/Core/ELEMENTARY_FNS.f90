!--------------------------------------------------------------------------------------
!
!Copyright (C) 2014 TopCoder Inc., All Rights Reserved.
!
!--------------------------------------------------------------------------------------

!--------------------------------------------------------------------------------------
!
!   Copyright 2014 Ecole Centrale de Nantes, 1 rue de la No�, 44300 Nantes, France
!
!   Licensed under the Apache License, Version 2.0 (the "License");
!   you may not use this file except in compliance with the License.
!   You may obtain a copy of the License at
!
!       http://www.apache.org/licenses/LICENSE-2.0
!
!   Unless required by applicable law or agreed to in writing, software
!   distributed under the License is distributed on an "AS IS" BASIS,
!   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!   See the License for the specific language governing permissions and
!   limitations under the License. 
!
!   Contributors list:
!   - G. Delhommeau
!   - P. Guevel
!   - J.C. Daubisse
!   - J. Singh  
!
!--------------------------------------------------------------------------------------


!   This module contains utility mathematics function used throughout the code
!
! Changes in version 1.5 (Implementation of Higher Order Panel Methods)
!       Added GET_BSPLINE function and compute_wave subroutine
!
!   @author yedtoss
!   @version 1.2
MODULE ELEMENTARY_FNS
    USE COMMON_TYPE

    IMPLICIT NONE
 
CONTAINS
    !-----------------------------------------------------------------------

    COMPLEX FUNCTION GG(Z,CEX)
        COMPLEX:: Z,CEX,Y
        REAL::T

        IF(REAL(Z)+16. <= 0) THEN       ! case 1 pp 368
            Y=1./Z
            GG=Y*(1.+Y*(-1.+Y*(2.+Y*(-6.+Y*(24.+Y*(-120.))))))
            RETURN
        END IF

        T=AIMAG(Z)
        IF(ABS(T)-10. > 0) THEN          !case 3 pp 368
            GG=0.711093/(Z+0.415775)+0.278518/(Z+2.29428)+0.010389/(Z+6.2900)
            IF(T < 0) THEN
                GG=GG-(0.,3.14159265)*CEX
                RETURN

            ELSE
                GG=GG+(0.,3.14159265)*CEX
                RETURN

            END IF
        END IF

        IF(REAL(Z)+0.5 <= 0) THEN !case 2 pp 368
            IF(T < 0) THEN

                    !   Case 4 pp. 369
                GG=((((((( (1.000000,1.3935496E-06)*Z+ (15.82958,-20.14222))  &
                    &*Z+ (-70.52863,-227.9511))*Z+ (-985.4221,-226.6272))*Z        &
                    &+ (-1202.318,1580.907))*Z+ (953.2441,1342.447))*Z             &
                    &+ (417.3716,-196.6665))*Z+ (-9.881266,-24.24952))/            &
                    &(((((((( (1.000000,0.0000000E+00)*Z+ (16.83184,-20.14481))*Z  &
                    &+ (-55.66969,-248.1167))*Z+ (-1068.640,-434.4707))*Z          &
                    &+ (-2082.250,1522.471))*Z+ (383.3455,2730.378))*Z             &
                    &+ (1216.791,351.7189))*Z+ (115.3926,-161.2647))*Z             &
                    &+ (-3.777369,-4.510900))-(0.,3.14159265)*CEX
                RETURN

            ELSE

                    !   Case 4 pp. 369
                GG=((((((( (1.000000,-1.3935496E-06)*Z+ (15.82958,20.14222))  &
                    &*Z+ (-70.52863,227.9511))*Z+ (-985.4221,226.6272))*Z          &
                    &+ (-1202.318,-1580.907))*Z+ (953.2441,-1342.447))*Z           &
                    &+ (417.3716,196.6665))*Z+ (-9.881266,24.24952))/              &
                    &(((((((( (1.000000,0.0000000E+00)*Z+ (16.83184,20.14481))*Z   &
                    &+ (-55.66969,248.1167))*Z+ (-1068.640,434.4707))*Z            &
                    &+ (-2082.250,-1522.471))*Z+ (383.3455,-2730.378))*Z           &
                    &+ (1216.791,-351.7189))*Z+ (115.3926,161.2647))*Z             &
                    &+ (-3.777369,4.510900))+(0.,3.14159265)*CEX
                RETURN

            END IF

        END IF


        GG=-(CLOG(Z)*(.1E+01+Z*(0.23721365E+00+Z*(0.206543E-01+Z*(0.763297E-03+&
            & Z*0.97087007E-05))))+0.5772156649E+00*(0.99999207E+00+&
            & Z*(-0.149545886E+01+Z*(0.41806426E-01+Z*(-0.3000591E-01+&
            & Z*(0.19387339E-02+Z*(-0.51801555E-03)))))))/(0.1E+01+&
            & Z*(-0.76273617E+00+Z*(0.28388363E+00+Z*(-0.66786033E-01+Z*(0.12982719E-01&
            & +Z*(-0.8700861E-03+Z*0.2989204E-03))))))

        IF(T < 0) THEN
            GG=GG-(0.,3.14159265)*CEX
            RETURN

        ELSE
            GG=GG+(0.,3.14159265)*CEX
            RETURN

        END IF

    END FUNCTION
    !-------------------------------------------------------------------------------!

    REAL FUNCTION CIH(AK,Z,H)
        REAL:: AK,Z,H

        IF((AK*H.LE.20).AND.(AK*H.GT.0.))THEN
            CIH=COSH(AK*(Z+H))/COSH(AK*H)
        ELSE
            CIH=EXP(AK*Z)
        END IF
        RETURN
    END FUNCTION
    !-------------------------------------------------------------------------------!

    REAL FUNCTION SIH(AK,Z,H)
        REAL:: AK,Z,H

        IF((AK*H.LE.20).AND.(AK*H.GT.0.))THEN
            SIH=SINH(AK*(Z+H))/COSH(AK*H)
        ELSE
            SIH=EXP(AK*Z)
        END IF
        RETURN
    END FUNCTION
!-------------------------------------------------------------------------------!



!------------------------------------------------------------------!
    !!  Evaluate the indicated B-spline function at given point.      !!
    !------------------------------------------------------------------!
    !!  INPUT  : PT - Point, KVEC - A valid knot vector,              !!
    !!           INDX - Index of B-spline functions starting from 0   !!
    !!  OUTPUT : "REAL(8)" B-spline function value                    !!
    !------------------------------------------------------------------!

    REAL FUNCTION GET_BSPLINE(PT, KVEC, IN_INDX)

        REAL, INTENT(IN) :: PT
        TYPE(KNOT_VECTOR), INTENT(IN) :: KVEC
        INTEGER, INTENT(IN) :: IN_INDX

        REAL :: N(0:KVEC%POLY_ORDER), SAVED, TEMP
        REAL,PARAMETER :: tolerance = 1e-8
        INTEGER :: I, J, K, L, INDX

        ! The given index is 1 based and we need 0 based
        INDX = IN_INDX -1

        ! Check if the given point is between min. of knot values and max. of knot values.
        IF (INDX==0 .AND. ABS(PT-KVEC%KNOTS(0))<TOLERANCE) THEN
            GET_BSPLINE = 1.0D0
        ELSE IF (INDX==(KVEC%LENGTH-KVEC%POLY_ORDER-1) .AND. ABS(PT-KVEC%KNOTS(KVEC%LENGTH))<TOLERANCE) THEN
            GET_BSPLINE = 1.0D0
        ELSE IF ((PT<KVEC%KNOTS(INDX)-TOLERANCE) .OR. (PT>=KVEC%KNOTS(INDX+KVEC%POLY_ORDER+1)+TOLERANCE)) THEN
            GET_BSPLINE = 0.0D0
        ELSE
            ! Check if the given point is in the corresponding knot span.
            DO J = 0, KVEC%POLY_ORDER
                IF (PT>=KVEC%KNOTS(INDX+J) .AND. PT<KVEC%KNOTS(INDX+J+1)) THEN
                    N(J) = 1.0D0
                ELSE
                    N(J) = 0.0D0
                ENDIF
            ENDDO
            ! Compute the B-spline function corresponding to the INDX using the recursive formular
            DO K = 1, KVEC%POLY_ORDER
                IF (ABS(N(0))<TOLERANCE) THEN
                    SAVED = 0.0D0
                ELSE
                    SAVED = ((PT-KVEC%KNOTS(INDX))*N(0)) / (KVEC%KNOTS(INDX+K)-KVEC%KNOTS(INDX))
                ENDIF
                DO J = 0, KVEC%POLY_ORDER-K
                    IF (ABS(N(J+1))<TOLERANCE) THEN
                        N(J) = SAVED
                        SAVED = 0.0D0
                    ELSE
                        TEMP = N(J+1) / (KVEC%KNOTS(INDX+J+K+1)-KVEC%KNOTS(INDX+J+1))
                        N(J) = SAVED + (KVEC%KNOTS(INDX+J+K+1)-PT)*TEMP
                        SAVED = (PT-KVEC%KNOTS(INDX+J+1))*TEMP
                    ENDIF
                ENDDO
            ENDDO
            GET_BSPLINE = N(0)
        ENDIF

    END FUNCTION GET_BSPLINE


    !   Calculate the complex potential, pressure and fluid velocities for a regular wave eta=sin(k*wbar-wt)
    SUBROUTINE Compute_Wave(k,w,beta,x,y,z,Phi,p,Vx,Vy,Vz,Environment)
    IMPLICIT NONE
    REAL :: k,w,beta,x,y,z
    COMPLEX :: Phi,p,Vx,Vy,Vz
    TYPE(TEnvironment) :: Environment
    REAL :: wbar
    COMPLEX,PARAMETER :: II=CMPLX(0.,1.)
    wbar=(x-Environment%XEFF)*COS(Beta)+(y-Environment%YEFF)*SIN(Beta)
    Phi=-Environment%g/w*CIH(k,z,Environment%Depth)*CEXP(II*k*wbar)
    p=-Environment%rho*Environment%g*II*CIH(k,z,Environment%Depth)*CEXP(II*k*wbar)
    Vx=-Environment%g/w*II*k*COS(beta)*CIH(k,z,Environment%Depth)*CEXP(II*k*wbar)
    Vy=-Environment%g/w*II*k*SIN(beta)*CIH(k,z,Environment%Depth)*CEXP(II*k*wbar)
    Vz=-Environment%g/w*k*SIH(k,z,Environment%Depth)*CEXP(II*k*wbar)
    END SUBROUTINE Compute_wave
END MODULE
