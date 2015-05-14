!--------------------------------------------------------------------------------------
!
!Copyright (C) 2014 TopCoder Inc., All Rights Reserved.
!
!--------------------------------------------------------------------------------------

!   This module solves a boundary element method for bodies with thin submerged elements.
!   It works with a finite or infinite fluid depth. This module is using
!   the lower order panel method. The bodies is assumed to be composed of two parts
!   which can be empty. One part consists of conventional element S_s, the other is a
!   surface consisting of thin 'dipoles' elements. To avoid singularity of the
!   original potential formulation of the problem when the thickness of part of the body
!   tends to zero, this method used an alternative integral equation. The
!   alternative form is equation 15.46 and 15.47 of reference R2 (or equation 6.2 and 6.3
!   of reference R1).
!   It is a pair of two simultaneous integral equation which is solved with a technique
!   similar to the one describe in 15.3 of reference R2 or to Chapter 6 of reference R1
!
!                       REFERENCES
!   R1:  A Rankine panel method as a tool for the hydrodynamic design of complex marine vehicles.
!   http://dspace.mit.edu/bitstream/handle/1721.1/50364/39696248.pdf?sequence=1
!
!   R2 Wamit User Manual http://www.wamit.com/manualupdate/V70_manual.pdf
!
!   Contest Dipoles Implementation in NEMOH version 1.0
!
!   @author yedtoss
!   @version 1.0
module SOLVE_BEM_DIPOLES_THIN

    USE COMPUTE_GREEN_FD
    USE COMPUTE_GREEN_INFD
    USE COM_VAR
    USE COMMON_TYPE
    USE OUTPUT
    USE UTILITY
    USE COMPUTE_GREEN_FREESURFACE

    implicit none

contains


    complex function compute_green_dipoles_wrapper(var, param)
        ! This function will compute the  of the derivative of the green function with respect to n_xi with n_xi the
        ! normal at field point coordinate. It is expected to be derived (again) with respect to n_x with
        ! and n_x the normal at source point coordinate.
        ! Arguments
        !   var(input) represent the variable with respect to which we are computing and we will
        !           derive this function with respect to n_x. It can represent either the x or y
        !           or z coordinate.
        !   param (input/temporary vars): common variables that contains useful variables
        !        See param definition for explanation of the member used here
        !

        REAL, INTENT(IN) :: var
        TYPE(ParamsCommonInf) :: param
        REAL:: xx, yy, zz
        INTEGER:: ISP, IFP
        REAL:: PCOS,PSIN
        REAL, PARAMETER :: PI=4.*ATAN(1.)

        ISP = param%ISP
        IFP = param%IFP
        xx = param%xx
        yy = param%yy
        zz = param%zz

        if (ISP > 0 .and. ISP <= IMX) THEN
            ! Getting centroid coordinate
            xx = XG(ISP) ! Source point x
            yy = YG(ISP) ! Source point y
            zz = ZG(ISP) ! Source point z
        end if

        if(param%coordinate == 0) THEN

            xx = var
        ELSE IF(param%coordinate == 1) THEN
            yy = var

        ELSE
            zz = var

        END IF

        IF(param%is_infinite == 1)  THEN ! infinite case
            call VAVINFD(2, 1,xx, yy, zz,ISP,IFP, param%SolverVar)  !1/r+1/r1

            call VNSINFD(2, 1,ISP,IFP,xx,yy,zz, param%SolverVar)

        ELSE
            call VAVFD(2, 2,xx,yy,zz,ISP,IFP, param%SolverVar)  !1/r+1/r1
            call VNSFD(2, param%AM0,param%AMH,param%NEXP,ISP,IFP,xx,yy,zz,  param%SolverVar)

        END IF

        PCOS=param%SolverVar%VSXM1*XN(IFP)+param%SolverVar%VSYM1*YN(IFP)+param%SolverVar%VSZM1*ZN(IFP)
        PSIN=param%SolverVar%VSXM2*XN(IFP)+param%SolverVar%VSYM2*YN(IFP)+param%SolverVar%VSZM2*ZN(IFP)

        compute_green_dipoles_wrapper = CMPLX(PCOS,PSIN) * (-4*PI)


    end function


    complex function compute_green_double_dipoles(xx, yy, zz, ISP, IFP, param)

        ! This function will compute the double derivative of the green function with respect to n_xi and n_x with n_xi
        ! the normal at field point coordinate andn_x the normal at source point coordinate.
        ! Arguments
        !   ISP (input) index of the source point panel.
        !   IFP (input) index of the field point panel
        !   xx, yy, zz(input). Coordinate of source point to be used if ISP is invalid
        !   param (input/temporary vars): common variables that contains useful variables
        !        See param definition for explanation of the member used here
        !

        INTEGER :: ISP ! source point

        INTEGER :: IFP !Field point

        REAL:: xx, yy, zz

        complex:: dxx, dyy, dzz

        TYPE(ParamsCommonInf) :: param

        REAL, PARAMETER :: PI=4.*ATAN(1.)

        IF(ISP > 0 .and. ISP <= IMX) THEN
            ! Getting centroid coordinate
            xx = XG(ISP) ! Source point x
            yy = YG(ISP) ! Source point y
            zz = ZG(ISP) ! Source point z

        END IF

        param%ISP = ISP
        param%IFP = IFP
        param%xx = xx
        param%yy = yy
        param%zz = zz

        ! Derivate with respect to x
        param%coordinate = 0
        dxx = deriv_richardson_complex(compute_green_dipoles_wrapper,xx, 1, 1e-3, param)

        ! Derivate with respect to y
        param%coordinate = 1
        dyy = deriv_richardson_complex(compute_green_dipoles_wrapper,yy, 1, 1e-3, param)

        ! Derivate with respect to z
        param%coordinate = 2
        dzz = deriv_richardson_complex(compute_green_dipoles_wrapper,zz, 1, 1e-3, param)

        compute_green_double_dipoles = XN(ISP)*dxx + YN(ISP)*dyy + ZN(ISP)*dzz



    end function compute_green_double_dipoles

    SUBROUTINE SOLVE_POTENTIAL_DIPOLES_THIN(param, PRESSURE, PHI, MeshFS, Switch_FS,  Switch_Potential, NDS, Force, &
        & out_potential)
        !  This subroutine solves a boundary element method for a finite or infinite depth using
        !   a three dimensional dipoles implementation.
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

        REAL:: xx, yy, zz
        REAL:: PCOS,PSIN
        INTEGER:: ISP, IFP
        COMPLEX, DIMENSION(:, :), ALLOCATABLE :: ZIJ
        COMPLEX, DIMENSION(:), POINTER :: ZPB,ZPS
        TYPE(TMesh) :: MeshFS
        COMPLEX :: PRESSURE(:)
        COMPLEX :: PHI(:)
        INTEGER :: Switch_FS, Switch_Potential
        TYPE(TempVar), POINTER :: SolverVar
        REAL, POINTER :: SP1,SM1,SP2,SM2
        REAL, POINTER :: VSXP1,VSXP2,VSYP1,VSYP2,VSZP1,VSZP2
        REAL, POINTER :: VSXM1,VSXM2,VSYM1,VSYM2,VSZM1,VSZM2
        REAL,DIMENSION(:,:) :: NDS
        COMPLEX,DIMENSION(:) :: Force
        COMPLEX,PARAMETER :: CII=CMPLX(0.,1.)
        COMPLEX,PARAMETER :: CRR=CMPLX(0.,1.)
        REAL :: out_potential(:)
        REAL, PARAMETER :: PI=4.*ATAN(1.)
        TYPE(ParamsCommonInf), TARGET :: param
        COMPLEX,DIMENSION(:),ALLOCATABLE :: NVEL ! normal velocity at centroid of a patch
        INTEGER:: I, J
        INTEGER, DIMENSION(:), ALLOCATABLE:: IPIV
        REAL, DIMENSION(:), ALLOCATABLE:: S, RWORK
        INTEGER:: INFO, RANK, LWORK
        COMPLEX, DIMENSION(:), ALLOCATABLE:: WORK
        COMPLEX, DIMENSION(:, :), ALLOCATABLE :: A
        COMPLEX, DIMENSION(:, :), ALLOCATABLE :: B

        SolverVar => param%SolverVar

        ZPB => SolverVar%ZPB
        ZPS => SolverVar%ZPS

        SP1 => SolverVar%SP1
        SM1 => SolverVar%SM1
        SP2 => SolverVar%SP2
        SM2 => SolverVar%SM2
        VSXP1 => SolverVar%VSXP1
        VSXP2 => SolverVar%VSXP2
        VSYP1 => SolverVar%VSYP1
        VSYP2 => SolverVar%VSYP2
        VSZP1 => SolverVar%VSZP1
        VSZP2 => SolverVar%VSZP2
        VSXM1 => SolverVar%VSXM1
        VSXM2 => SolverVar%VSXM2
        VSYM1 => SolverVar%VSYM1
        VSYM2 => SolverVar%VSYM2
        VSZM1 => SolverVar%VSZM1
        VSZM2 => SolverVar%VSZM2

        ALLOCATE(ZIJ(IMX, IMX+1))
        ALLOCATE(IPIV(IMX))
        ALLOCATE(A(IMX, IMX))
        ALLOCATE(B(IMX, 1))


        ZIJ = cmplx(0., 0.)

        NVEL = param%NVEL

        xx = 0.
        yy = 0.
        zz = 0.

        CALL COMPUTE_NEXP(param)



        DO ISP=1,IMX

            DO IFP=1,IMX

                IF(param%is_thin_body(ISP) == 1) THEN

                    ! Using left hand side of equation (15.47)  Reference
                    ! double derivative of green function multiply by Area of IFP

                    ZIJ(ISP, IFP) =  ZIJ(ISP, IFP) + compute_green_double_dipoles(xx, yy, zz, ISP, IFP, param)


                    ! Using right hand side of equation (15.47)  Reference second term
                    IF(param%is_infinite == 1)  THEN ! infinite case
                        call VAVINFD(1, 1,XG(ISP),YG(ISP),ZG(ISP),ISP,IFP, SolverVar)  !1/r+1/r1

                        call VNSINFD(1, 1,ISP,IFP,XG(ISP),YG(ISP),ZG(ISP), SolverVar)

                    ELSE
                        call VAVFD(1, 2,XG(ISP),YG(ISP),ZG(ISP),ISP,IFP, param%SolverVar)  !1/r+1/r1
                        call VNSFD(1, param%AM0,param%AMH,param%NEXP,ISP,IFP,XG(ISP),YG(ISP),ZG(ISP),  param%SolverVar)

                    END IF

                    PCOS=VSXM1*XN(ISP)+VSYM1*YN(ISP)+VSZM1*ZN(ISP)
                    PSIN=VSXM2*XN(ISP)+VSYM2*YN(ISP)+VSZM2*ZN(ISP)
                    ZIJ(ISP, IMX + 1 ) = ZIJ(ISP, IMX + 1 ) + NVEL(IFP) * CMPLX(PCOS,PSIN) * (-4*PI)


                ELSE
                    ! Using left hand side of equation (15.46)  Reference
                    ! integral of green function derivative with respect to normal of field point
                    IF(param%is_infinite == 1)  THEN ! infinite case
                        call VAVINFD(2, 1,XG(ISP),YG(ISP),ZG(ISP),ISP,IFP, SolverVar)  !1/r+1/r1

                        call VNSINFD(2, 1,ISP,IFP,XG(ISP),YG(ISP),ZG(ISP), SolverVar)

                    ELSE
                        call VAVFD(2, 2,XG(ISP),YG(ISP),ZG(ISP),ISP,IFP, param%SolverVar)  !1/r+1/r1
                        call VNSFD(2, param%AM0,param%AMH,param%NEXP,ISP,IFP,XG(ISP),YG(ISP),ZG(ISP),  param%SolverVar)

                    END IF

                    PCOS=VSXM1*XN(IFP)+VSYM1*YN(IFP)+VSZM1*ZN(IFP)
                    PSIN=VSXM2*XN(IFP)+VSYM2*YN(IFP)+VSZM2*ZN(IFP)
                    ZIJ(ISP, IFP) =  ZIJ(ISP, IFP) +  CMPLX(PCOS,PSIN) * (-4*PI)

                    IF(ISP == IFP) THEN

                        ! Using left hand side of equation (15.46)  Reference, first term
                        ZIJ(ISP, IFP) =  ZIJ(ISP, IFP) + 2*PI*CRR

                    END IF

                    ! Using right hand side of equation (15.46)  Reference
                    ! No need to recompute green function here. It has already beeen computed
                    ! in previous instructions.

                    ZIJ(ISP, IMX + 1 ) = ZIJ(ISP, IMX + 1 ) + NVEL(IFP) * CMPLX(SP1+SM1,SP2+SM2) * (-4*PI)

                END IF

            END DO

            IF(param%is_thin_body(ISP) == 1) THEN

                ! Using right hand side of equation (15.47)  Reference
                ! -4pi*nvel + intregal nvel green function double derivative

                ZIJ(ISP, IMX + 1 ) = ZIJ(ISP, IMX + 1 ) -4*PI*NVEL(ISP)

            END IF

        END DO




        ! Setting up matrix A and B of the linear equations Ax=B
        DO I = 1, IMX

            DO J=1, IMX

                A(I, J) = ZIJ(I,J)

            END DO

            B(I, 1) = ZIJ(I, IMX+1)

        END DO


        !------------------------------------------------!
        ! Using lapack linear solver with partial pivoting
        CALL CGESV( IMX, 1, A, IMX, IPIV, B, IMX, INFO )

        IF( INFO.GT.0 ) THEN
            WRITE(*,*)'The diagonal element of the triangular factor of A,'
            WRITE(*,*)'U(',INFO,',',INFO,') is zero, so that'
            WRITE(*,*)'A is singular; the solution could not be computed.'
            STOP
        END IF
        !------------------------------------------------!



        DO I=1,IMX
            ZPS(I) = cmplx(0., 0.)
            ZPB(I) = B(I, 1)
        END DO

           !       Assemble pressure at patches centroid
        DO J=1,IMX
            PRESSURE(J)=RHO*CII*param%W*ZPB(J) !*AIRE(i)
        END DO



          !       Save free surface elevation
        IF (Switch_FS.EQ.1) THEN
            DO j=1,MeshFS%Npoints
                PHI(j) = COMPUTE_POTENTIAL_DOMAIN_THIN(MeshFS%X(1,j),MeshFS%X(2,j),0.,param, B)
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

    END SUBROUTINE SOLVE_POTENTIAL_DIPOLES_THIN


    COMPLEX FUNCTION COMPUTE_POTENTIAL_DOMAIN_THIN(xx, yy, zz, param, B)
        ! Compute the potential at a given point in the fluid domain or in the boundary
        ! for the thin dipoles formulation
        ! Arguments
        ! xx,yy,zz: Cartesian coordinate of the given point
        ! B: Array, solution of the linear system. Contain the phi coefficients
        ! param: common variables that contains useful variables
        !        See param definition for explanation of the member used here
        TYPE(ParamsCommonInf), TARGET :: param
        COMPLEX :: tmp
        REAL :: xx, yy, zz ! The coordinate of a point in the domain of the given patch
        REAL :: uf,vf
        INTEGER:: IFP, ISP
        COMPLEX, DIMENSION(:, :), ALLOCATABLE :: B ! The solution to the linear system
        TYPE(TempVar), POINTER :: SolverVar
        REAL, POINTER :: SP1,SM1,SP2,SM2
        REAL, POINTER :: VSXP1,VSXP2,VSYP1,VSYP2,VSZP1,VSZP2
        REAL, POINTER :: VSXM1,VSXM2,VSYM1,VSYM2,VSZM1,VSZM2
        REAL:: PCOS,PSIN
        REAL, PARAMETER :: PI=4.*ATAN(1.)

        SolverVar => param%SolverVar

        SP1 => SolverVar%SP1
        SM1 => SolverVar%SM1
        SP2 => SolverVar%SP2
        SM2 => SolverVar%SM2
        VSXP1 => SolverVar%VSXP1
        VSXP2 => SolverVar%VSXP2
        VSYP1 => SolverVar%VSYP1
        VSYP2 => SolverVar%VSYP2
        VSZP1 => SolverVar%VSZP1
        VSZP2 => SolverVar%VSZP2
        VSXM1 => SolverVar%VSXM1
        VSXM2 => SolverVar%VSXM2
        VSYM1 => SolverVar%VSYM1
        VSYM2 => SolverVar%VSYM2
        VSZM1 => SolverVar%VSZM1
        VSZM2 => SolverVar%VSZM2

        COMPUTE_POTENTIAL_DOMAIN_THIN = cmplx(0, 0)

        ISP = 1  ! Disable ISP panel

        DO IFP=1,IMX

            SP1=0.
            SP2=0.
            SM1=0.
            SM2=0.
            VSXM1 = 0.
            VSYM1 = 0.
            VSZM1 = 0.
            VSXM2 = 0.
            VSYM2 = 0.
            VSZM2 = 0.
            ! Using right hand side of equation (15.46)  Reference
            IF(param%is_infinite == 1)  THEN ! infinite case
                CALL VVV(2, 1,IFP,xx,yy,zz, SolverVar) ! finite part
                CALL VNV(2, IFP,xx,yy,zz, SolverVar)  !infinite part


            ELSE
                CALL VVV(2, 2,IFP,xx,yy,zz, SolverVar)
                CALL VNVF(2, param%AM0,param%AMH,param%NEXP,IFP,xx,yy,zz, SolverVar)

            END IF

            tmp =  param%NVEL(IFP) * CMPLX(SP1+SM1,SP2+SM2) * (-4*PI)

            ! Using left hand side of equation (15.46)  Reference
            ! integral of green function derivative with respect to normal of field point
            ! No need to recompute green function here. It has already beeen computed
            ! in previous instructions.

            PCOS=VSXM1*XN(IFP)+VSYM1*YN(IFP)+VSZM1*ZN(IFP)
            PSIN=VSXM2*XN(IFP)+VSYM2*YN(IFP)+VSZM2*ZN(IFP)
            tmp =  tmp -  B(IFP, 1)* CMPLX(PCOS,PSIN) * (-4*PI)

            COMPUTE_POTENTIAL_DOMAIN_THIN = COMPUTE_POTENTIAL_DOMAIN_THIN + tmp
        END DO

        COMPUTE_POTENTIAL_DOMAIN_THIN = COMPUTE_POTENTIAL_DOMAIN_THIN/(2.*PI)



    END FUNCTION COMPUTE_POTENTIAL_DOMAIN_THIN



end module SOLVE_BEM_DIPOLES_THIN
