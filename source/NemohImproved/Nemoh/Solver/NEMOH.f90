!--------------------------------------------------------------------------------------
!
!Copyright (C) 2014-2016 TopCoder Inc., All Rights Reserved.
!
!--------------------------------------------------------------------------------------

!--------------------------------------------------------------------------------------
!
!   NEMOH V1.0 - BVP solver - January 2014
!
!--------------------------------------------------------------------------------------
!
!   Copyright 2014 Ecole Centrale de Nantes, 1 rue de la No, 44300 Nantes, France
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
!   - A. Babarit  
!
!--------------------------------------------------------------------------------------
!


!   This is the main program for the Nemoh software.
!   Nemoh is a Boundary Element Methods (BEM) code dedicated to the computation of first order wave loads
!   on offshore structures (added mass, radiation damping, diffraction forces).
!   Unlike other BEM softwares, Nemoh's approach decouples the resolution of the linear free surface
!   Boundary Value Problem (BVP) and the definition of the boundary condition on the body (body condition).
!   This feature makes it easy to deal with flexible structure, hydroelasticity, generalised modes and
!   unconventional degrees of freedom with Nemoh.
!   Changes in version 1.2
!           Remove  IO and make this module contains a subroutine instead of a program
!   Changes in version 1.3 (Drift forces and QTF Implementation of Nemoh)
!           Added parameter to store drift forces and yaw moment
!
! Changes in version 1.4 (Code Acceleration of the Calculation of Influence Coefficients of Nemoh)
!       Allow fast code for influence coefficients by accepting a switch and passing it to SOLVE_BVP
!       subroutine
!
! Changes in version 1.5 (Implementation of Higher Order Panel Methods)
!       Added additional parameters and control to COMPUTE_NEMOH subroutine to handle higher order panel method
!       Added CHECK_IS_WAVE and GET_RAD_IDX function

! Changes in version 1.6 (Dipoles Implementation in NEMOH)
!       Added additional parameters and control to COMPUTE_NEMOH subroutine to handle dipoles implementation
!
! Changes in version 1.7 (Hydrodynamic Data Exporter Assembly v1.0)
!       Run the Mesh processor in order to fillmesh properties like center of buoyancy,
!       volume displacement, stiffness matrix
!
! Changes in version 1.8 (Irregular Frequencies Assembly)
!       Implement the removal of irregular frequencies using the potential formulation
!       with the extended integral equations method.  Added additional parameters to
!       control this implementation.
!
! Changes in version 1.9 (OpenWarp - Add Logging Functionality version 1.0)
!       Added support for logging.
!
!   @author yedtoss
!   @version 1.9
#include "logging.h"
MODULE NEMOH

    IMPLICIT NONE

    CONTAINS

    SUBROUTINE COMPUTE_NEMOH (ct_RHO, ct_G, ct_Depth, ct_XEFF, ct_YEFF, ct_ZEFF, ct_Indiq_solver,&
    &ct_maxIterations, ct_restartParam, ct_TOL_GMRES, mesh_P, mesh_X, mesh_cPanel, mesh_XM,&
    &mesh_N, mesh_A, mesh_Isym, Npoints, Npanels, Nbodies, Nproblems, NBCPanels,&
    & bc_NormalVelocity, bc_Omega, bc_Switch_Potential,&
    & bc_Switch_FreeSurface, bc_Switch_Kochin, bc_Switch_Type, NDS, Nintegration, Theta,  Ntheta,&
    & NFSPoints, NFSPanels,meshfs_X, meshfs_P, out_PHI,out_PRESSURE, out_HKochin, line, &
    & out_potential, n_potentials, n_tabulatedx, n_tabulatedz, n_points_simpson, drift_forces, &
    & yaw_moment, FastInfluence_Switch, higher_panel_method, n_beta, n_radiation, rad_case, &
    & beta, num_panel_higher_order, b_spline_order, is_thin_body, use_dipoles_implementation, &
    & Switch_YAW_MOMENT, Switch_DRIFT_FORCES, center_buoyancy, displacement, waterplane_area, stifness, &
    & is_interior_domain, remove_irregular_frequencies, log_level) BIND(C)

    USE MIdentification
    USE MMesh
    USE MBodyConditions
    !    USE iflport
    USE SOLVE_BEM
    USE OUTPUT
    USE INITIALIZATION
    USE COM_VAR
    USE OMP_LIB
    USE ISO_C_BINDING
    USE MESH_STATISTICS
    USE FLOGGING

    !

        ! --- Environment --------------------------------------------------
    ! Volumique mass of the fluid in KG/M**3
    REAL(C_FLOAT) :: ct_RHO
    ! Gravity
    REAL(C_FLOAT) :: ct_G
    ! depth of the domain
    REAL(C_FLOAT) :: ct_Depth
    ! Coordinates where the waves are measured
    REAL(C_FLOAT) :: ct_XEFF,ct_YEFF,ct_ZEFF

    ! Which solver: (0) direct solver, (1): GMRES
    INTEGER(C_INT):: ct_Indiq_solver

    INTEGER(C_INT) :: ct_maxIterations, ct_restartParam
    REAL(C_FLOAT) :: ct_TOL_GMRES

    INTEGER(C_INT) :: Npoints,Npanels,Nbodies, Nproblems, NBCPanels, Nintegration, NFSPoints, NFSPanels
     INTEGER(C_INT) :: n_beta, n_radiation
    INTEGER(C_INT) :: Ntheta, n_potentials, n_tabulatedx, n_tabulatedz, n_points_simpson

    INTEGER(C_INT) :: mesh_Isym, higher_panel_method, b_spline_order, num_panel_higher_order

    REAL(C_FLOAT) :: mesh_X(3, Npoints)                  ! Nodes coordinates
    REAL(C_FLOAT) :: mesh_N(3, Npanels)                  ! Normal vectors
    REAL(C_FLOAT) :: mesh_XM(3, Npanels)                 ! Centre of panels
    INTEGER(C_INT) :: mesh_P(4, Npanels)               ! Connectivities
    INTEGER(C_INT) :: mesh_cPanel(Npanels)            ! To which body belongs the panel
    REAL(C_FLOAT) :: mesh_A(Npanels)                    ! Area of panel
    INTEGER(C_INT) :: is_thin_body(Npanels)  ! whether a panel is a thin body or not
    INTEGER(C_INT) :: use_dipoles_implementation ! 1 to use thin dipoles implementation when available, 0 otherwise
    INTEGER(C_INT) :: remove_irregular_frequencies ! 1 to remove irregular frequencies, 0 otherwise
    INTEGER(C_INT) :: Switch_DRIFT_FORCES, Switch_YAW_MOMENT ! to compute drift forces and yaw moment respectively, 0 otherwise

    INTEGER(C_INT) :: is_interior_domain(Npanels) ! Array indicating whether or not a panel is in the body or the interior free surface domain.
                                                  ! 1 to indicate the interior, 0 for the body.

    INTEGER(C_INT) :: meshfs_P(4, NFSPanels)               ! Connectivities
    REAL(C_FLOAT) :: meshfs_X(3, NFSPoints)                  ! Nodes coordinates

    REAL(C_FLOAT) :: rad_case(n_radiation, 8)  ! The radiation case
    REAL(C_FLOAT) :: beta(n_beta)  ! The wave directions

    REAL(C_FLOAT) :: bc_Omega(Nproblems)
    COMPLEX(C_FLOAT_COMPLEX) :: bc_NormalVelocity(NBCPanels, Nproblems)
    INTEGER(C_INT) :: bc_Switch_Potential(Nproblems)
    INTEGER(C_INT) :: bc_Switch_FreeSurface(Nproblems)
    INTEGER(C_INT) :: bc_Switch_Kochin(Nproblems)
    INTEGER(C_INT) :: bc_Switch_Type(Nproblems)
    INTEGER(C_INT) :: FastInfluence_Switch(Nproblems)
    INTEGER(C_INT) :: log_level

    REAL(C_FLOAT) :: NDS(Nintegration, NBCPanels)
    REAL(C_FLOAT) :: Theta(Ntheta)

    REAL(C_FLOAT) :: line(Nintegration, Nproblems*2)
    REAL(C_FLOAT) :: out_potential(Nproblems, n_potentials)
    COMPLEX(C_FLOAT_COMPLEX) :: out_PHI(Nproblems, 1+NFSPoints)
    COMPLEX(C_FLOAT_COMPLEX) :: out_PRESSURE(Nproblems, NBCPanels)
    COMPLEX(C_FLOAT_COMPLEX) :: out_HKochin(Nproblems, Ntheta)
    REAL(C_FLOAT) :: drift_forces(Nproblems, Ntheta, 2)  ! The drift forces
    REAL(C_FLOAT) :: yaw_moment(Nproblems, Ntheta)  ! The yaw moment

    REAL(C_FLOAT) :: center_buoyancy(Nbodies , 3)
    REAL(C_FLOAT) :: displacement(Nbodies)
    REAL(C_FLOAT):: waterplane_area(Nbodies)
    REAL(C_FLOAT) :: stifness(Nbodies, 6, 6)
    REAL :: tmp_par(Nproblems) ! temporary variable per problem to be used in parallel loop
    INTEGER :: tmp_idx(Nproblems) ! temporary variable per problem to be used in parallel loop


    !   ID
    TYPE(TID)               :: ID
    !   Array sizes
    !INTEGER                 :: NFA,NSYMY    ! Number of panels and symmetry about xOz plane
    !   Body conditions
    TYPE(TBodyConditions)   :: BodyConditions
    !   Kochin function

    !REAL,DIMENSION(:),ALLOCATABLE :: Theta
    !   Meshes
    TYPE(TMesh) :: Mesh 
    TYPE(TMesh) :: MeshFS
    !   Nemoh
    !REAL                    :: T
    !   Results
    !REAL :: XEFF,YEFF
    !INTEGER :: Nintegration, NpanIsy
    INTEGER ::  NpanIsy
    !REAL,DIMENSION(:,:),ALLOCATABLE :: NDS
    COMPLEX,DIMENSION(:,:),ALLOCATABLE :: Force
    !REAL,DIMENSION(:),ALLOCATABLE :: line
    !   Locals
    REAL                    :: PI
    INTEGER                 :: c,i,j, WhereToSave,k
    INTEGER :: NumTargetProblem
    INTEGER :: ntemp ! temporary variable

    REAL:: stat_xg, stat_yg, stat_zg
    INTEGER:: stat_np, stat_nface
    REAL, DIMENSION(:), ALLOCATABLE :: stat_x, stat_y, stat_z
    INTEGER, DIMENSION(:, :), ALLOCATABLE :: stat_p
    INTEGER:: stat_tmp1, stat_tmp2

    ! Enable date output
    call log_set_default_output_date(.true.)

    ! Set log level
    call log_set_default_log_level(log_level)

    
    !log(info,*) "here, have some info"

    !WRITE(0,*) ' Start'
    if (logp(LOG_INFO)) then
        write(logu,*) trim(logl(LOG_INFO)), " Simulation started"
    endif

    !
    !   --- Initialisation -------------------------------------------------------------------------------------
    !
    PI=4.*ATAN(1.)
    !WRITE(*,*) ' '
    !WRITE(*,'(A)', advance='NO') '  -> Initialisation '
    if (logp(LOG_INFO)) then
        write(logu,*) trim(logl(LOG_INFO)), " -> Initialisation"
    endif

    ! Executing mesh statistics

    DO I = 1, Nbodies

        stat_zg = rad_case(6*I, 8)
        stat_yg = rad_case(6*I, 7)
        stat_xg = rad_case(6*I, 6)
        stat_nface = count(mesh_cPanel == I)
        stat_np = 4*stat_nface
        ALLOCATE(stat_x(stat_np), stat_y(stat_np), stat_z(stat_np), stat_p(stat_nface, 4))

        stat_tmp1 =   0
        stat_tmp2 =   1

        ! Extracting connectivities information for body I

        DO J = 1, Npanels
            IF (mesh_cPanel(J) == I .and. is_interior_domain(J) == 0) THEN
                stat_tmp1 = stat_tmp1 + 1
                DO K = 1, 4
                    stat_x(stat_tmp2 + K -1) = mesh_X(1, mesh_P(K, J))
                    stat_y(stat_tmp2 + K -1) = mesh_X(2, mesh_P(K, J))
                    stat_z(stat_tmp2 + K -1) = mesh_X(3, mesh_P(K, J))
                END DO
                stat_p(stat_tmp1, 1) = stat_tmp2
                stat_p(stat_tmp1, 2) = stat_tmp2 + 1
                stat_p(stat_tmp1, 1) = stat_tmp2 + 2
                stat_p(stat_tmp1, 2) = stat_tmp2 + 3
                stat_tmp2 = stat_tmp2 + 4
            END IF

        END DO

        ! Computing the volume displacement, center of buoyancy and  stiffness matrix

        CALL COMPUTE_MESH_STATISTICS(mesh_Isym, ct_XEFF,ct_YEFF, 0., 1., stat_xg, stat_yg, stat_zg, &
        & stat_np, stat_nface, stat_x, stat_y, stat_z, stat_p, ct_RHO, ct_G, &
        & center_buoyancy(I, :), displacement(I), waterplane_area(I), stifness(I, :, :))

    END DO


    TABULATION_JZ = n_tabulatedz
    TABULATION_IR = n_tabulatedx
    TABULATION_NPINTE = n_points_simpson

    IF (TABULATION_NPINTE < 100) THEN
        if (logp(LOG_ERROR)) then
            write(logu,*) trim(logl(LOG_ERROR)), " Please choose a number of sub intervals greater than 100"
        endif
        !WRITE(*,*)'Error: Please choose a number of sub intervals greater than 100'
        STOP
    END IF

    IF (TABULATION_NPINTE >= 1000000) THEN
        !WRITE(*,*)'Error: Please choose a number of sub intervals less than 1 000 000'
        !WRITE(*,*)'because all surfaces will be thin'
        if (logp(LOG_ERROR)) then
            write(logu,*) trim(logl(LOG_ERROR)), " Please choose a number of sub intervals greater than 100"
            write(logu,*) trim(logl(LOG_ERROR)), " because all surfaces will be thin"
        endif
        STOP
    END IF

    IF (TABULATION_IR < 328) THEN
        !WRITE(*,*)'Error: Please choose a number of points in x direction greater than or equal to 328'
        if (logp(LOG_ERROR)) then
            write(logu,*) trim(logl(LOG_ERROR)), " Please choose a number of points in x direction greater than or equal to 328"
        endif
        STOP
    END IF

    IF (TABULATION_JZ < 46) THEN
        !WRITE(*,*)'Error: Please choose a number of points in z direction greater than or equal to  46'
        if (logp(LOG_ERROR)) then
            write(logu,*) trim(logl(LOG_ERROR)), " Please choose a number of points in z direction greater than or equal to  46"
        endif
        STOP
    END IF


    CALL CreateTMesh(Mesh,Npoints,Npanels,Nbodies)
    CALL CreateTBodyConditions(BodyConditions,Nproblems,Npanels)

    Mesh%Isym = mesh_Isym


    Mesh%X = mesh_X
    Mesh%N = mesh_N
    Mesh%XM = mesh_XM
    Mesh%P = mesh_P
    Mesh%cPanel = mesh_cPanel
    Mesh%A = mesh_A

    BodyConditions%Nproblems=Nproblems
    BodyConditions%Npanels=NBCPanels

    BodyConditions%Omega = bc_Omega
    BodyConditions%Switch_Type = bc_Switch_Type
    BodyConditions%Switch_Potential = bc_Switch_Potential
    BodyConditions%Switch_FreeSurface = bc_Switch_FreeSurface
    BodyConditions%Switch_Kochin = bc_Switch_Kochin

    BodyConditions%NormalVelocity = bc_NormalVelocity

    !   Initialise Nemoh
    CALL INITIALIZE(ID,NFA,NSYMY,XEFF,YEFF,Mesh, ct_RHO, ct_G, ct_Depth, ct_XEFF, ct_YEFF, ct_ZEFF, ct_Indiq_solver,&
    & ct_maxIterations, ct_restartParam, ct_TOL_GMRES)

    is_rankine_dipoles_computed = 0
    is_rankine_sources_computed = 0

    IF(use_dipoles_implementation==1 .and. higher_panel_method == 1) THEN
        higher_panel_method = 0
        !write(*, *) 'Warning: Higher panel method does not support yet the thin dipoles implementation'
        !write(*, *) 'Warning: Falling back to the low order method.'
        if (logp(LOG_WARNING)) then
            write(logu,*) trim(logl(LOG_WARNING)), " Higher panel method does not support yet", &
            & "the thin dipoles implementation"
            write(logu,*) trim(logl(LOG_WARNING)), " Falling back to the low order method."
        endif
        
    END IF

    IF(remove_irregular_frequencies==1 .and. higher_panel_method == 1) THEN
        higher_panel_method = 0
        !write(*, *) 'Warning: Higher panel method does not support yet the removal of irregular frequencies'
        !write(*, *) 'Warning: Falling back to the low order method.'
        if (logp(LOG_WARNING)) then
            write(logu,*) trim(logl(LOG_WARNING)), " Higher panel method does not support yet", &
            & "the removal of irregular frequencies"
            write(logu,*) trim(logl(LOG_WARNING)), " Falling back to the low order method."
        endif
    END IF

    IF(use_dipoles_implementation==1 .and. remove_irregular_frequencies == 1) THEN
        use_dipoles_implementation = 0
        !write(*, *) 'Warning: Dipoles Implementation and Irregular Frequencies removal are not compatible'
        !write(*, *) 'Warning: Disabling Dipoles Implementation'
        if (logp(LOG_WARNING)) then
            write(logu,*) trim(logl(LOG_WARNING)), " Dipoles Implementation and Irregular", &
            & "Frequencies removal are not compatible"
            write(logu,*) trim(logl(LOG_WARNING)), " Disabling Dipoles Implementation"
        endif
    END IF
    ! Don't precompute green function tabulation if we want to use ode in the infinite depth case
    ntemp = 0
    DO I = 1, Nproblems
        IF(FastInfluence_Switch(I) == 0) THEN

            ntemp = ntemp  + 1
            ! EXIT do not exit anymore as we want to know if all problems use normal influence
            ! computation

        END IF
    END DO
    IF (Depth.NE.0. .or. ntemp /= 0 .or. higher_panel_method == 1 .or. use_dipoles_implementation == 1 &
    & .or. remove_irregular_frequencies == 1) THEN
        CALL CREK(TABULATION_NPINTE, TABULATION_IR, TABULATION_JZ)
    END IF

    IF ((higher_panel_method == 1 .or. use_dipoles_implementation == 1 &
    & .or. remove_irregular_frequencies == 1) .and. ntemp /= Nproblems) THEN
        !write(*, *) 'Warning: Higher panel method does not support yet the computation'
        !write(*, *) ' of influence coefficients using the ode technique. Falling back to default mode.'
        if (logp(LOG_WARNING)) then
            write(logu,*) trim(logl(LOG_WARNING)), " Higher panel method does not support", &
            & "yet the computation"
            write(logu,*) trim(logl(LOG_WARNING)), " of influence coefficients using the ode", &
            & "technique. Falling back to default mode."
        endif
    END IF

    !WRITE(*,'(A)', advance='NO') '.'
    !   Initialise Force matrix


    !   Initialise free surface calculation points

    MeshFS%Npoints = NFSPoints
    MeshFS%Npanels = NFSPanels
    IF (MeshFS%Npoints.GT.0) THEN
        CALL CreateTMesh(MeshFS,NFSPoints,NFSPanels,1)
        MeshFS%X = meshfs_X
        MeshFS%P = meshFS_P
    END IF
    !   Initialise results table
    ALLOCATE(Force(Nintegration,Bodyconditions%Nproblems))
    Force(:,:)=0.
    !WRITE(*,*) '. Done !'
    !WRITE(*,*) ' '
    if (logp(LOG_INFO)) then
        write(logu,*) trim(logl(LOG_INFO)), " -> Initialisation is done!"
    endif
    !


    NpanIsy = Mesh%Npanels*2**Mesh%Isym

    ALLOCATE(ProblemSavedAt(Bodyconditions%Nproblems), ProblemToSaveLocation(Bodyconditions%Nproblems))

    ProblemSavedAt = -3
    ProblemToSaveLocation = -1
    WhereToSave = 1
    NumTargetProblem = 0



    DO i=1,BodyConditions%Nproblems

        IF(ProblemSavedAt(i) == -3) THEN

            ProblemSavedAt(i) = -2

            DO j=i+1, BodyConditions%Nproblems

                IF(ProblemSavedAt(j) /= -3) THEN
                    CYCLE
                END IF

                IF(ABS(BodyConditions%Omega(j)-BodyConditions%Omega(i)) <= 1.E-20  ) THEN

                    ProblemSavedAt(i) = -1
                    ProblemSavedAt(j) = WhereToSave

                END IF


            END DO

            IF(ProblemSavedAt(i) == -1) THEN

                ProblemToSaveLocation(i) = WhereToSave

                WhereToSave = WhereToSave + 1
                NumTargetProblem = NumTargetProblem + 1

            END IF

            IF(ProblemSavedAt(i) == -2) THEN

                NumTargetProblem = NumTargetProblem + 1

            END IF

        END IF

    END DO

    IF(WhereToSave > 1) THEN

        ALLOCATE(InfluenceMatrixCache(WhereToSave-1,2,NFA,NFA))
        ALLOCATE(CQCache(WhereToSave-1, 101), QQCache(WhereToSave-1, 101))
        ALLOCATE( SM1Cache(WhereToSave-1, NFA, NFA), SM2Cache(WhereToSave-1,NFA, NFA))
        ALLOCATE( SP1Cache(WhereToSave-1, NFA, NFA), SP2Cache(WhereToSave-1,NFA, NFA))
        ALLOCATE(ARCache(WhereToSave-1, 31), AMBDACache(WhereToSave-1, 101))
        ALLOCATE(NEXPCache(WhereToSave-1))


    END IF

    !   --- Solve BVPs and calculate forces ---------------------------------------------------------
    !
    !WRITE(*,*) ' -> Solve BVPs and calculate forces '
    !WRITE(*,*) ' '
    if (logp(LOG_INFO)) then
        write(logu,*) trim(logl(LOG_INFO)), " -> Solving BVPs and calculate forces"
    endif



    !$OMP PARALLEL DO
    DO j=1,BodyConditions%Nproblems

        !       Solve BVP

        IF(ProblemSavedAt(j) == -2 .or. ProblemSavedAt(j) == -1) THEN

            tmp_idx(j) = GET_RAD_IDX(J, n_beta, n_radiation)
            tmp_par(j) = 0
            if(n_beta > 0 .and. tmp_idx(j) <= n_beta) tmp_par(j) = beta(tmp_idx(j))

            CALL SOLVE_BVP(j,ID,2.*PI/BodyConditions%Omega(j),&
                & BodyConditions%Switch_Kochin(j),NTheta,Theta,&
                & BodyConditions%Switch_Freesurface(j),MeshFS,BodyConditions%Switch_Potential(j), &
                & Nintegration, NpanIsy, NDS, Force, BodyConditions%NormalVelocity, omp_get_thread_num(), &
                & out_PHI(j, :), out_PRESSURE(j, :), out_HKochin(j, :), out_potential(j, :), drift_forces(j,:, :), &
                & yaw_moment(j,:), FastInfluence_Switch(j), higher_panel_method, tmp_par(j),  &
                & Mesh, CHECK_IS_WAVE(J, n_beta, n_radiation), GET_RAD_IDX(J, n_beta, n_radiation),  &
                & rad_case, num_panel_higher_order, b_spline_order, is_thin_body, use_dipoles_implementation, &
                & Switch_DRIFT_FORCES, Switch_YAW_MOMENT, is_interior_domain, remove_irregular_frequencies)

            !WRITE(0,'(A,I5,A,I5,A)', advance='YES') ' Problem ',j,' / ',BodyConditions%Nproblems,' ... Done !'
            if (logp(LOG_INFO)) then
                write(logu, '(A, A,I5,A,I5,A)', advance='YES') trim(logl(LOG_INFO)), ' Problem ',j,' / ', &
                & BodyConditions%Nproblems,' ... Done !'
            endif



        END IF


    END DO
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO
    DO j=1,BodyConditions%Nproblems 

        IF(ProblemSavedAt(j) >= 1) THEN

            tmp_idx(j) = GET_RAD_IDX(J, n_beta, n_radiation)
            tmp_par(j) = 0
            if(n_beta > 0 .and. tmp_idx(j) <= n_beta) tmp_par(j) = beta(tmp_idx(j))

            !       Solve BVP
            CALL SOLVE_BVP(j,ID,2.*PI/BodyConditions%Omega(j),&
                & BodyConditions%Switch_Kochin(j),NTheta,Theta,&
                & BodyConditions%Switch_Freesurface(j),MeshFS,BodyConditions%Switch_Potential(j), &
                & Nintegration, NpanIsy, NDS, Force, BodyConditions%NormalVelocity, omp_get_thread_num(), &
                & out_PHI(j, :), out_PRESSURE(j, :), out_HKochin(j, :), out_potential(j, :), drift_forces(j,:, :), &
                & yaw_moment(j,:), FastInfluence_Switch(j), higher_panel_method, tmp_par(j),  &
                & Mesh, CHECK_IS_WAVE(J, n_beta, n_radiation), GET_RAD_IDX(J, n_beta, n_radiation), &
                & rad_case, num_panel_higher_order, b_spline_order, is_thin_body, use_dipoles_implementation, &
                & Switch_DRIFT_FORCES, Switch_YAW_MOMENT, is_interior_domain, remove_irregular_frequencies)

            !WRITE(0,'(A,I5,A,I5,A)', advance='YES') ' Problem ',j,' / ',BodyConditions%Nproblems,' ... Done !'
            if (logp(LOG_INFO)) then
                write(logu, '(A, A,I5,A,I5,A)', advance='YES') trim(logl(LOG_INFO)), ' Problem ',j,' / ', &
                & BodyConditions%Nproblems,' ... Done !'
            endif

        END IF
    END DO
    !$OMP END PARALLEL DO


    !WRITE(*,*) ' ' 


    !
    !   --- Save results -------------------------------------------------------
    !
    DO c=1,Nintegration
        DO j=1,BodyConditions%Nproblems
            IF (Bodyconditions%Switch_type(j).NE.-1) THEN
                line(c, 2*j-1)=ABS(Force(c,j))
                line(c, 2*j)=ATAN2(AIMAG(Force(c,j)),REAL(Force(c,j)))
            ELSE
                line(c, 2*j-1)=-AIMAG(Force(c,j))/BodyConditions%Omega(j)
                line(c, 2*j)=REAL(Force(c,j))
            END IF
        END DO
    END DO
    !
    !   --- Finalize -----------------------------------------------
    !


    DEALLOCATE(Force)
    DEALLOCATE(ProblemSavedAt, ProblemToSaveLocation, InfluenceMatrixCache)
    DEALLOCATE(CQCache, QQCache, SM1Cache, SM2Cache, SP1Cache, SP2Cache)
    DEALLOCATE(NEXPCache, ARCache, AMBDACache)
    CALL DEALLOCATE_DATA
    !WRITE(0,*) 'NEMOH Solver completed.'
    if (logp(LOG_INFO)) then
        write(logu,*) trim(logl(LOG_INFO)), " NEMOH Solver completed."
    endif

    END SUBROUTINE COMPUTE_NEMOH

    INTEGER FUNCTION CHECK_IS_WAVE(J, n_beta, n_radiation)
    INTEGER:: is_wave, J, tmp, n_beta, n_radiation

    tmp = modulo(J-1, n_beta+n_radiation)

    IF (tmp >= n_beta) THEN
        CHECK_IS_WAVE = 0
    ELSE
        CHECK_IS_WAVE = 1
    END IF

    END FUNCTION

    INTEGER FUNCTION GET_RAD_IDX(J, n_beta, n_radiation)
    ! This return the radiation index or the beta index
    ! Radiation index is return when tmp>=n_beta
    ! and beta index otherwise
    INTEGER:: rad_case_idx, J, tmp, n_beta, n_radiation

    tmp = modulo(J-1, n_beta+n_radiation)

    IF (tmp >= n_beta) THEN
        GET_RAD_IDX = modulo(tmp - n_beta, n_radiation) + 1
    ELSE
        GET_RAD_IDX = modulo(tmp, n_beta) + 1
    END IF


    END FUNCTION
!
END MODULE NEMOH
!----------------------------------------------------------------

SUBROUTINE CREK(NPINTE, n_tabulatedx, n_tabulatedz)
    USE COM_VAR
    IMPLICIT NONE
    INTEGER:: I,J,JZ,IR,NPINTE, n_tabulatedx, n_tabulatedz
    REAL:: AKZ,AKR
    JZ= n_tabulatedz
    IR= n_tabulatedx
    DO J=1,JZ
        IF(J .LE. 46) THEN
            AKZ=MIN(10**(J/5.-6),10**(J/8.-4.5),16.)
        ELSE
            AKZ = (16./(JZ - 46)) * (J - 46) + 1.5E6
        ENDIF
        XZ(J)=-AKZ
    END DO
    XR(1)=0.
    DO I=2,IR
        IF(I.LT.40)THEN
            AKR=MIN(10**((I-1.)/5-6),4./3.+ABS((I-32.)/3.))
        ELSE IF(I .LE. 328) THEN
            AKR=4./3.+ABS((I-32.)/3.)
        ELSE
            AKR = (100./(IR- 328)) *(I-328)
        ENDIF
        XR(I)=AKR
    END DO

    DO J=1,JZ
        DO I=1,IR
            CALL VNS(NPINTE,XZ(J),XR(I),I,J)
        END DO
    END DO

    RETURN

END SUBROUTINE
!---------------------------------------------------------

SUBROUTINE VNS(NPINTE,AKZ,AKR,I,J)
    USE COM_VAR
    USE ELEMENTARY_FNS
    IMPLICIT NONE
    INTEGER:: I,J,NPINTE,IT
    REAL:: AKZ,AKR,PI,CT,ST,CQIT,TETA
    REAL:: QQT(NPINTE),CQT(NPINTE)
    REAL:: FD1JX,FD1JZ,FD2JX,FD2JZ
    COMPLEX:: IM,C1,C2,ZIK,GZ,CEX
    IM=(0.,1.)
    PI=4.*ATAN(1.)
    CALL COFINT(NPINTE,CQT,QQT)
    FD1JX=0.
    FD1JZ=0.
    FD2JX=0.
    FD2JZ=0.
    DO IT=1,NPINTE
        TETA=QQT(IT)
        CQIT=CQT(IT)
        CT=COS(TETA)
        ST=SIN(TETA)
        ZIK=AKZ+IM*AKR*CT
        IF(REAL(ZIK)+30. <= 0) THEN
            CEX=(0.,0.)
        ELSE
            CEX=CEXP(ZIK)
        END IF

        GZ=GG(ZIK,CEX)
        C1=CQIT*(GZ-1./ZIK)
        C2=CQIT*CEX
        FD1JX=FD1JX+CT*AIMAG(C1)
        FD1JZ=FD1JZ+REAL(C1)
        FD2JX=FD2JX+CT*AIMAG(C2)
        FD2JZ=FD2JZ+REAL(C2)
    END DO
    APD1X(I,J)=FD1JX
    APD1Z(I,J)=FD1JZ
    APD2X(I,J)=FD2JX
    APD2Z(I,J)=FD2JZ
    RETURN

END SUBROUTINE
!-------------------------------------------------------------------

SUBROUTINE COFINT(NPINTE,CQT,QQT)
    USE COM_VAR
    USE FLOGGING
    IMPLICIT NONE
    INTEGER :: J,NPINTE, num_thin
    REAL:: PI,QQT(NPINTE),CQT(NPINTE)

    PI=4.*ATAN(1.)

    num_thin = 0

    DO J=1,NPINTE
        QQT(J)=-PI/2.+(J-1.)/(NPINTE-1.)*PI
        IF(J-1 <= 0) THEN
            CQT(J)=PI/(3.*(NPINTE-1.))
        ELSE IF(J-NPINTE >= 0) THEN
            CQT(J)=PI/(3.*(NPINTE-1.))
        ELSE IF(MOD(J,2) /= 0) THEN
            CQT(J)=2./(3.*(NPINTE-1.))*PI
        ELSE
            CQT(J)=4./(3.*(NPINTE-1.))*PI
        END IF

        IF (ABS(CQT(J)) < 1e-6) THEN
            num_thin = num_thin +1

        END IF
    END DO

    IF (num_thin > 0) THEN
        !write(*, *) 'Error, there are ', num_thin, ' sub intervals too small. There will make calculations wrong'
        !write(*, *) 'Please choose a smaller value of the number of sub intervals'
        if (logp(LOG_ERROR)) then
            write(logu,*) trim(logl(LOG_ERROR)), " there are ', num_thin, ' sub intervals too small.", &
            & " There will make calculations wrong"
            write(logu,*) trim(logl(LOG_ERROR)), " Please choose a smaller value of", &
            & "the number of sub intervals"
        endif
        STOP
    END IF

    RETURN
END SUBROUTINE
!------------------------------------------------------------------


