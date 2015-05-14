!--------------------------------------------------------------------------------------
!
!Copyright (C) 2014 TopCoder Inc., All Rights Reserved.
!
!--------------------------------------------------------------------------------------

!  Module to contains common types
!  All common types should be moved from COM_VAR to COMMON_TYPE later on.
!
! Contest Code Acceleration of the Calculation of Influence Coefficients of Nemoh
!
! Changes in version 1.1 (Implementation of Higher Order Panel Methods)
!       Added additional types to handle higher order panel methods and move some types from COM_VAR
!
! Changes in version 1.2 (Irregular Frequencies Assembly)
!      Added is_infinite member to ParamsCommonInf type indicating
!      whether or not to remove irregular frequencies.
!
!   @author yedtoss
!   @version 1.2

module COMMON_TYPE

    USE MMesh
    implicit none

    ! Definition of TYPE Environment
TYPE TEnvironment
    REAL :: RHO         ! Sea water density
    REAL :: G           ! Gravity constant
    REAL :: Depth       ! Water depth
    REAL :: XEFF,YEFF   ! Coordinates of point where the incident wave is measured
END TYPE TEnvironment

    !!  Knot vector including information about the length of knot vector and degree of the B-spline
    TYPE KNOT_VECTOR
        INTEGER :: LENGTH = -1, POLY_ORDER = -1
        REAL, DIMENSION(:), ALLOCATABLE :: KNOTS   ! Should start from 0 index
    END TYPE KNOT_VECTOR


    TYPE TempVar

        ! --- Boundary value problem ---------------------------------------
        ! normal velocity array as input
        !    REAL, DIMENSION(:), ALLOCATABLE :: NVEL
        ! period
        REAL :: T
        ! Computed and computed potential on original (B) and symmetric boundary(S)
        COMPLEX, DIMENSION(:), ALLOCATABLE :: ZPB,ZPS
        ! Source ditribution
        COMPLEX, DIMENSION(:), ALLOCATABLE :: ZIGB,ZIGS

        ! --- Solver ------------------------------------------------------

        ! Linear complex matrix to be solved
        COMPLEX, DIMENSION(:, :), ALLOCATABLE :: ZIJ
        REAL :: FSP,FSM,VSXP,VSYP,VSZP,VSXM,VSYM,VSZM
        ! Variable for storage of Greens function
        REAL :: SP1,SM1,SP2,SM2
        ! Variable for storage of Greens function
        REAL :: VSXP1,VSXP2,VSYP1,VSYP2,VSZP1,VSZP2
        REAL :: VSXM1,VSXM2,VSYM1,VSYM2,VSZM1,VSZM2

        INTEGER:: NQ
        REAL:: CQ(101),QQ(101),AMBDA(31),AR(31)

        INTEGER :: ProblemNumber, Switch_FastInfluence


        ! High order specific parameter
        INTEGER :: Num_Patches ! Number of patches
        INTEGER , DIMENSION(:), ALLOCATABLE:: Num_Basis_U, Num_Basis_V ! Number of spline basis function in (u, v)

        INTEGER , DIMENSION(:), ALLOCATABLE:: Panel_Idx ! Panel_Idx(i) is the index of the first panel for patch i

        TYPE(KNOT_VECTOR) , DIMENSION(:), ALLOCATABLE:: KNOTS_U ! knot vectors of U size of num_patches

        TYPE(KNOT_VECTOR) , DIMENSION(:), ALLOCATABLE:: KNOTS_V ! knot vectors of V size of num_patches

        INTEGER :: Num_Unknows  ! number total of unknows

    END TYPE TempVar

    ! Type to keep parameters to pass to function to integrate
    ! Currently only for drift forces and yaw moment
    TYPE ParamsCommonInf

        REAL :: kwave  ! Wave number
        REAL :: THeta
        REAL  :: AMH  ! m0
        REAL  :: AM0  ! m0
        REAL::AKH  ! A
        REAL:: beta  ! Wave incident angle
        REAL:: w ! the frequency
        REAL ::  r

        INTEGER :: i,j,k, ii, jj, kk ! variables to keep hold of index data
        INTEGER:: ISP ! variable to keep hold of source index
        INTEGER:: IFP ! variable to keep hold of field index
        REAL:: Z ! the Z coordinate
        REAL:: c ! constant

        ! High Panel parameter
        REAL :: uf, vf ! to keep hold of temporary parametric space coordinate. Usually for the field point
        REAL :: xx, yy, zz ! To keep hold of temporary cartesian coordinate. Usually for the field point
        INTEGER :: use_cartesian_coordinate ! 1 to use cartesian coordinate directly, 0 to convert u,v to cartesian
        INTEGER :: compute_potential ! 1 to compute the potential, 0 to use the B-Spline basis function
        INTEGER:: panel_idxf, panel_idx ! Keep hold of index of field point panel and source point panel respectively
        TYPE(TempVar) :: SolverVar ! Variable to keep the boundary value problem temporary variables
        COMPLEX,DIMENSION(:),ALLOCATABLE :: NVEL ! normal velocity at centroid of a patch
        INTEGER :: NEXP ! temporary variable
        COMPLEX, DIMENSION(:, :), ALLOCATABLE :: B ! The solution to the linear system

        INTEGER:: is_wave ! 1 to indicate that this is a wave problem, 0 for other problem

        TYPE(TEnvironment) :: Environment ! The environment of the problem
        TYPE(TMesh) :: Mesh ! The mesh to solve

        REAL, ALLOCATABLE:: rad_case(:, :) ! The radiation case

        INTEGER :: ProblemNumber ! The index of current problem to solve
        INTEGER:: rad_case_idx  ! The radiation case index to use for current problem

        INTEGER:: n_integration, NpanIsy ! Number of integration and number of patches
        INTEGER :: KSPLIN ! order of B-Spline
        INTEGER:: NSPLIN  ! Number of B-Spline per patches


        INTEGER, ALLOCATABLE:: is_thin_body(:)  ! whether a panel is a thin body or not
        INTEGER, ALLOCATABLE:: is_interior_domain(:) ! Array indicating whether or not a panel is in the body
                                                     ! or the interior free surface domain. 1 to indicate the interior, 0 for the body.

        INTEGER:: coordinate ! represent the coordinate to use

        INTEGER:: is_infinite ! 1 if this is an infinite problem, 0 otherwise




    END TYPE ParamsCommonInf





end module COMMON_TYPE
