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
!
!--------------------------------------------------------------------------------------


! This module will intialize common (read only for problems) variables shared by all problems
!  Changes in version 1.2:
!        Make sure this module takes it's input from variable instead of reading it
!
! Changes in version 1.3 (Code Acceleration of the Calculation of Influence Coefficients of Nemoh)
!       Remove initialization of tabulated green function from here. It will be perform in Nemoh.f90
!       when not using ODE for influence coefficients
!
!   @author yedtoss
!   @version 1.3
MODULE INITIALIZATION

IMPLICIT NONE

CONTAINS
!---------------------------------------------------------------------------
    SUBROUTINE INITIALIZE(ID,NF,NSYM,XF,YF,Mesh, ct_RHO, ct_G, ct_Depth, ct_XEFF, ct_YEFF, ct_ZEFF, ct_Indiq_solver,&
    & ct_maxIterations, ct_restartParam, ct_TOL_GMRES)
!
    USE MIDENTIFICATION
    USE COM_VAR
    USE PREPARE_MESH
    USE MMesh

!
    IMPLICIT NONE

      ! --- Environment --------------------------------------------------
    ! Volumique mass of the fluid in KG/M**3
    REAL :: ct_RHO
    ! Gravity
    REAL :: ct_G
    ! depth of the domain
    REAL :: ct_Depth
    ! Coordinates where the waves are measured
    REAL :: ct_XEFF,ct_YEFF,ct_ZEFF

    ! Which solver: (0) direct solver, (1): GMRES
    INTEGER:: ct_Indiq_solver

    INTEGER :: ct_maxIterations, ct_restartParam
    REAL :: ct_TOL_GMRES

!   ID
    TYPE(TID) :: ID
!   Geometry
    INTEGER :: NF,NSYM
    REAL :: XF,YF
    TYPE(TMesh) :: Mesh

    RHO = ct_RHO
    G = ct_G
    DEPTH = ct_DEPTH
    XF = ct_XEFF
    YF = ct_YEFF

    Indiq_solver = ct_Indiq_solver
    restartParam = ct_restartParam
    TOL_GMRES = ct_TOL_GMRES
    maxIterations = ct_maxIterations


    XEFF=XF
    YEFF=YF
    NFA=Mesh%Npanels
    NP=Mesh%Npoints
    NSYMY=Mesh%Isym
    IF (NSYMY.NE.1) NSYMY=0
    IF (NSYMY.EQ.1) YEFF=0
    NF=NFA
    NSYM=NSYMY
!   Initialise Nemoh
    CALL ALLOCATE_DATA
    CALL PRE_PROC_MESH(Mesh)
!
    END SUBROUTINE INITIALIZE  

END MODULE INITIALIZATION
