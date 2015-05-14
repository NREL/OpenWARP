!--------------------------------------------------------------------------------------
!
!Copyright (C) 2015 TopCoder Inc., All Rights Reserved.
!
!--------------------------------------------------------------------------------------

!   This module solves a boundary element method for bodies using the potential formulation,
!   and by removing irregular frequencies.
!   It works with a finite or infinite fluid depth. This module is using
!   the lower order panel method. It assumes a discretization of both the body and the
!   interior of the free surface. The contour of the interior free surface is determined using
!   waterline segments. Waterline is the line where the hull of a ship meets the surface of the water,
!   in concept or reality (See http://en.wikipedia.org/wiki/Waterline)
!   Panels from both discretisation are used to derive the a system
!   of two simultaneous linear equation to solve. The main reference for the implementation is R1
!
!                       REFERENCES
!   R1:  Irregular frequency removal from the boundary integral equation for the wave-body problem
!        http://dspace.mit.edu/bitstream/handle/1721.1/11691/32279180.pdf?sequence=1
!
!   Contest Irregular Frequencies Assembly
!
!   @author TCSASSEMBLER
!   @version 1.0


MODULE SOLVE_BEM_IRREGULAR

    USE COMPUTE_GREEN_FD
    USE COMPUTE_GREEN_INFD
    USE COM_VAR
    USE COMMON_TYPE
    USE OUTPUT
    USE UTILITY
    USE COMPUTE_GREEN_FREESURFACE

    IMPLICIT NONE

    CONTAINS

    SUBROUTINE SOLVE_POTENTIAL_IRREGULAR(param, PRESSURE, PHI, MeshFS, Switch_FS,  Switch_Potential, NDS, Force, &
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

        ALLOCATE(ZIJ(IMX, IMX + 1))
        ALLOCATE(IPIV(IMX))
        ALLOCATE(A(IMX, IMX))
        ALLOCATE(B(IMX, 1))


        ZIJ = cmplx(0., 0.)

        NVEL = param%NVEL

        xx = 0.
        yy = 0.
        zz = 0.

        CALL COMPUTE_NEXP(param)

        Force = 0.



        DO ISP=1,IMX

            DO IFP=1,IMX

                ! Using left hand side of equation (5.33) or (5.34)  Reference
                    ! integral of green function derivative with respect to normal of field point
                    IF(param%is_infinite == 1)  THEN ! infinite case
                        call VAVINFD(2, 1,XG(ISP),YG(ISP),ZG(ISP),ISP,IFP, SolverVar)  !1/r+1/r1

                        call VNSINFD(2, 1,ISP,IFP,XG(ISP),YG(ISP),ZG(ISP), SolverVar)

                    ELSE
                        call VAVFD(2, 2,XG(ISP),YG(ISP),ZG(ISP),ISP,IFP, param%SolverVar)  !1/r+1/r1
                        call VNSFD(2, param%AM0,param%AMH,param%NEXP,ISP,IFP,XG(ISP),YG(ISP),ZG(ISP),  param%SolverVar)

                    END IF

                    PCOS = VSXM1*XN(IFP) + VSYM1*YN(IFP) + VSZM1*ZN(IFP)
                    PSIN = VSXM2*XN(IFP) + VSYM2*YN(IFP) + VSZM2*ZN(IFP)

                    ! Using left hand side of equation (5.33) or (5.34)  Reference second or third term
                    ZIJ(ISP, IFP) =  ZIJ(ISP, IFP) +  CMPLX(PCOS,PSIN) * (-4*PI)

                    IF(ISP == IFP) THEN

                        IF(param%is_interior_domain(ISP) == 1) THEN

                            ! Using left hand side of equation (5.34)  Reference, first term
                        ZIJ(ISP, IFP) =  ZIJ(ISP, IFP) - 4*PI*CRR
                        ELSE

                        ! Using left hand side of equation (5.33)  Reference, first term
                        ZIJ(ISP, IFP) =  ZIJ(ISP, IFP) + 2*PI*CRR

                        ENDIF

                    END IF

                    ! Using right hand side of equation (5.33) or (5.34)  Reference
                    ! No need to recompute green function here. It has already beeen computed
                    ! in previous instructions.

                    ZIJ(ISP, IMX + 1 ) = ZIJ(ISP, IMX + 1 ) + NVEL(IFP)*CMPLX(SP1 + SM1, SP2 + SM2)*(-4*PI)

                    IF(ZIJ(ISP, IFP) /= ZIJ(ISP, IFP)) THEN

                        ZIJ(ISP, IFP) = 1e-4

                    ENDIF

                    IF(ZIJ(ISP, IMX + 1 ) /= ZIJ(ISP, IMX + 1 )) THEN

                        ZIJ(ISP, IMX + 1 ) = 1e-4

                    ENDIF



            END DO

        END DO

        ! Setting up matrix A and B of the linear equations Ax=B
        DO I = 1,IMX

            DO J=1,IMX

                A(I, J) = ZIJ(I,J)

            END DO

            B(I, 1) = ZIJ(I, IMX + 1)

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
            PRESSURE(J) = RHO*CII*param%W*ZPB(J) !*AIRE(i)
        END DO



          !       Save free surface elevation
        IF (Switch_FS.EQ.1) THEN
            DO j=1,MeshFS%Npoints
                PHI(j) = COMPUTE_POTENTIAL_DOMAIN_IRREGULAR(MeshFS%X(1,j),MeshFS%X(2,j),0.,param, B)
            END DO
        END IF




         !       Save output
        IF (Switch_Potential.EQ.1) THEN
            CALL SAVE_POTENTIAL(PRESSURE, out_potential)
        END IF

        !       Calculate force coefficients at centroid of panel
        DO i=1,param%n_integration
            DO j=1,param%NpanIsy
                Force(i) = Force(i) - PRESSURE(j)*NDS(i,j)
            END DO
        END DO

        ! Deallocate allocated arrays
        if (allocated(ZIJ)) deallocate(ZIJ)
        if (allocated(IPIV)) deallocate(IPIV)
        if (allocated(A)) deallocate(A)
        if (allocated(B)) deallocate(B)

    END SUBROUTINE SOLVE_POTENTIAL_IRREGULAR


    COMPLEX FUNCTION COMPUTE_POTENTIAL_DOMAIN_IRREGULAR(xx, yy, zz, param, B)
        ! Compute the potential at a given point in the fluid domain or in the boundary
        ! for the potential formulation with irregular frequency
        ! Arguments
        ! xx,yy,zz: Cartesian coordinate of the given point
        ! B: Array, solution of the linear system. Contain the phi coefficients
        ! param: common variables that contains useful variables
        !        See param definition for explanation of the member used here
        TYPE(ParamsCommonInf), TARGET :: param
        COMPLEX :: tmp
        REAL :: xx, yy, zz ! The coordinate of a point in the domain of the given patch. Point is on the free surface
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

        COMPUTE_POTENTIAL_DOMAIN_IRREGULAR = cmplx(0, 0)

        ISP = 1  ! Disable ISP panel

        DO IFP=1,IMX

            SP1 = 0.
            SP2 = 0.
            SM1 = 0.
            SM2 = 0.
            VSXM1 = 0.
            VSYM1 = 0.
            VSZM1 = 0.
            VSXM2 = 0.
            VSYM2 = 0.
            VSZM2 = 0.
            ! Using right hand side of equation (5.34)  Reference
            IF(param%is_infinite == 1)  THEN ! infinite case
                CALL VVV(2, 1, IFP, xx, yy, zz, SolverVar) ! finite part
                CALL VNV(2, IFP, xx, yy, zz, SolverVar)  !infinite part


            ELSE
                CALL VVV(2, 2, IFP, xx, yy, zz, SolverVar)
                CALL VNVF(2, param%AM0, param%AMH, param%NEXP, IFP, xx, yy, zz, SolverVar)

            END IF

            tmp =  param%NVEL(IFP)*CMPLX(SP1 + SM1, SP2 + SM2)*(-4*PI)

            ! Using left hand side of equation (5.34)  Reference, second or third term
            ! integral of green function derivative with respect to normal of field point
            ! No need to recompute green function here. It has already beeen computed
            ! in previous instructions.

            PCOS = VSXM1*XN(IFP) + VSYM1*YN(IFP) + VSZM1*ZN(IFP)
            PSIN = VSXM2*XN(IFP) + VSYM2*YN(IFP) + VSZM2*ZN(IFP)
            tmp =  tmp - B(IFP, 1)*CMPLX(PCOS, PSIN)*(-4*PI)

            COMPUTE_POTENTIAL_DOMAIN_IRREGULAR = COMPUTE_POTENTIAL_DOMAIN_IRREGULAR + tmp
        END DO

        ! Dividing by constant term in front of first term in equation (5.34)
        COMPUTE_POTENTIAL_DOMAIN_IRREGULAR = COMPUTE_POTENTIAL_DOMAIN_IRREGULAR/(-4.*PI)



    END FUNCTION COMPUTE_POTENTIAL_DOMAIN_IRREGULAR


END MODULE SOLVE_BEM_IRREGULAR
