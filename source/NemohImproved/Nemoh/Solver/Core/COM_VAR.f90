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
!   - J. Singh
!   - A. Babarit
!
!--------------------------------------------------------------------------------------


!   This module contains read only variables used by all problems. It also contains writeable variable.
!   These writeable variables are temporary variables used by each problem for computation. There are
!   grouped under the type TempVar and passed to each subroutine needing it.
!
!   Changes in version 1.2 (Drift forces and QTF Implementation of Nemoh)
!    Added ParamsCommon derived type. It should be pass to function for which we want the integrate or
!    derivative. Use to store additional paramaters
!
!   Changes in version 1.3 (Code Acceleration of the Calculation of Influence Coefficients of Nemoh)
!       Added variables for caching rankine values.
!
! Changes in version 1.4 (Implementation of Higher Order Panel Methods)
!       MOVE some types to COMMON_TYPE
!
! Changes in version 1.5 (Dipoles Implementation in NEMOH)
!       Remove ParamsCommon structure
!
!   @author yedtosss
!   @version 1.5
MODULE COM_VAR

    USE COMMON_TYPE

    IMPLICIT NONE

   !---- Tabulation ----------
   INTEGER:: TABULATION_JZ  ! Represents Number of points in z direction of tabulated data
   INTEGER:: TABULATION_IR  ! Represents Number of points in x direction of tabulated data
   INTEGER:: TABULATION_NPINTE ! Represents Number of sub intervals used
                               ! to approximate the green function integral using simpson rule

    ! --- Environment --------------------------------------------------
    ! Volumique mass of the fluid in KG/M**3
    REAL :: RHO
    ! Gravity
    REAL :: G
    ! depth of the domain
    REAL :: Depth
    ! Coordinates where the waves are measured
    REAL :: XEFF,YEFF,ZEFF

    ! --- Geometrical quantities ---------------------------------------
    ! Mesh file
    CHARACTER(80) :: MESHFILE
    INTEGER :: LFILE
    ! No. of points in the surface mesh
    INTEGER :: NP
    ! No of total panels in the surface mesh
    INTEGER :: NFA
    !
    INTEGER :: IMX
    INTEGER :: IXX
    !!!!!!!!!
    !Symmetry: 0 for no symmetry and1 for symmetry
    INTEGER:: NSYMY
    !
    REAL :: ZER
    ! DIST(I) is the maximum distance (projection on XY plane!) of a panel I from other points
    REAL, DIMENSION(:), ALLOCATABLE :: DIST,TDIS
    !vertices of the panel, size NFA(no of facettes)
    INTEGER, DIMENSION(:), ALLOCATABLE :: M1,M2,M3,M4
    ! Array of cordinates of points of size NP (no of points)
    REAL, DIMENSION(:), ALLOCATABLE :: X,Y,Z
    ! Normal vector
    REAL, DIMENSION(:), ALLOCATABLE :: XN,YN,ZN
    ! Array of centre of gravity for each panel
    REAL, DIMENSION(:), ALLOCATABLE :: XG,YG,ZG
    ! Array for surface area of the panels
    REAL, DIMENSION(:), ALLOCATABLE :: AIRE

    !Values used in interpolation of infinite part of the Greens function
    !REAL :: XR(328),XZ(46),APD1X(328,46),APD1Z(328,46),APD2X(328,46),APD2Z(328,46)
    REAL, DIMENSION(:), ALLOCATABLE :: XR, XZ
    REAL, DIMENSION(:, :), ALLOCATABLE :: APD1X, APD1Z, APD2X, APD2Z

    ! --- Reading and writing units -----------------------------------
    ! File for visualization of potential
    INTEGER :: Sav_potential

    ! Which solver: (0) direct solver, (1): GMRES
    INTEGER:: Indiq_solver

    INTEGER :: maxIterations, restartParam
    REAL :: TOL_GMRES

    ! Variable to avoid recomputing influence matrix when period does not change
    INTEGER, DIMENSION(:), ALLOCATABLE :: ProblemSavedAt, ProblemToSaveLocation
    COMPLEX, DIMENSION(:,:,:,:), ALLOCATABLE:: InfluenceMatrixCache

    COMPLEX, DIMENSION(:, :, : , :), ALLOCATABLE :: TSafCache
    INTEGER, DIMENSION(:, :, :), ALLOCATABLE :: TSipivCache
    REAL, DIMENSION(:, :, :), ALLOCATABLE :: TSrCache,TScCache, SM1Cache, SM2Cache, SP1Cache, SP2Cache
    REAL, DIMENSION(:, :), ALLOCATABLE :: CQCache, QQCache, AMBDACache, ARCache
    INTEGER, DIMENSION(:), ALLOCATABLE :: NEXPCache
    CHARACTER(1), DIMENSION(:,:),ALLOCATABLE :: TSequedCache


    REAL, DIMENSION(:, :), ALLOCATABLE :: rankine_sources_cache  ! Cache for rankine source integral
    REAL, DIMENSION(:, :), ALLOCATABLE :: rankine_dipoles_cache  ! Cache for rankine dipoles integral
    INTEGER:: is_rankine_sources_computed ! Variable to indicate if rankine source has been computed once
    INTEGER:: is_rankine_dipoles_computed ! Variable to indicate if rankine dipole has been computed once




    REAL,PARAMETER :: tolerance = 1e-8


CONTAINS


    SUBROUTINE ALLOCATE_TEMP(SolverVar)
        TYPE(TempVar) :: SolverVar
        ALLOCATE(SolverVar%ZIGB(NFA),SolverVar%ZIGS(NFA)) !sources
        ALLOCATE(SolverVar%ZPB(NFA),SolverVar%ZPS(NFA))   !potential
        !Linear complex matrix to be solved by direct solver
        !For GMRES workspace is calculated inside the main subroutine
        IF(Indiq_solver == 0) THEN
            ALLOCATE(SolverVar%ZIJ(NFA,NFA+1))
        ELSE
            ALLOCATE(SolverVar%ZIJ(NFA,NFA))
        END IF

    END SUBROUTINE ALLOCATE_TEMP


    SUBROUTINE DEALLOCATE_TEMP(SolverVar)
        TYPE(TempVar) :: SolverVar
        !if(Indiq_solver .eq. 0)DEALLOCATE(SolverVar%ZIJ)
        DEALLOCATE(SolverVar%ZIJ)
        DEALLOCATE(SolverVar%ZPB,SolverVar%ZPS)
        DEALLOCATE(SolverVar%ZIGB,SolverVar%ZIGS)

    END SUBROUTINE DEALLOCATE_TEMP
       
END MODULE COM_VAR
