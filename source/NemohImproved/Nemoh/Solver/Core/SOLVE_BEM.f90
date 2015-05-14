!--------------------------------------------------------------------------------------
!
!Copyright (C) 2014-2015 TopCoder Inc., All Rights Reserved.
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
!   - J. Singh
!   - A. Babarit 
!
!--------------------------------------------------------------------------------------


!   This module solves the linear BVP (Boundary Value Problem) for the potential and calculate pressure
!   field, hydrodynamic coefficients, far field coefficients and wave elevation.
!   Changes in version 1.2
!              Allow this function to return it's result to variables instead
!              of writing it to files
!
! Changes in version 1.3 (Drift forces and QTF Implementation of Nemoh)
!    Added call to subroutines computing the drift forces and the yaw moment.
!    The incident wave angles used are the same as the one given for the kochin function
!    computation
!
! Changes in version 1.4 (Code Acceleration of the Calculation of Influence Coefficients of Nemoh)
!       Allow fast code for influence coefficients by accepting a switch and assigning it to SolverVar
!       for use in other subroutines
!
! Changes in version 1.5 (Implementation of Higher Order Panel Methods)
!       Added additional parameters and control to SOLVE_BVP subroutine to handle higher order panel method
!
! Changes in version 1.6 (Dipoles Implementation in NEMOH)
!       Added additional parameters and control to SOLVE_BVP subroutine to handle dipoles implementation
!
! Changes in version 1.7 (Hydrodynamic Data Exporter Assembly v1.0)
!       Added a dedicated variable to toggle drift forces.
!
! Changes in version 1.8 (Irregular Frequencies Assembly)
!       Implement the removal of irregular frequencies using the potential formulation
!       with the extended integral equations method
!
!   @author yedtoss, TCSASSEMBLER
!   @version 1.8
MODULE SOLVE_BEM
    !
    IMPLICIT NONE
!
CONTAINS
    !
    SUBROUTINE SOLVE_BVP(ProblemNumber,ID,Period,Switch_Kochin,&
        &NTheta,Theta,Switch_FS,MeshFS,Switch_potential, &
        & Nintegration, NpanIsy, NDS, Force,NormalVelocity, ThreadId, PHI, PRESSURE, HKochin, &
        & out_potential, drift_forces, yaw_moment, Switch_FastInfluence, higher_panel_method, &
        & beta, Mesh, is_wave, rad_case_idx, rad_case, num_panel_higher_order, b_spline_order, &
        & is_thin_body, use_dipoles_implementation, Switch_DRIFT_FORCES, Switch_YAW_MOMENT, &
        & is_interior_domain, remove_irregular_frequencies)
        !
        USE MIDENTIFICATION
        USE COM_VAR
        USE MMesh
        USE OUTPUT
        USE SOLVE_BEM_INFD_DIRECT
        USE SOLVE_BEM_FD_DIRECT
        USE SOLVE_BEM_INFD_ITERATIVE
        USE SOLVE_BEM_FD_ITERATIVE
        USE KOCHIN
        USE POTENTIAL_DOMAIN
        USE COMPUTE_MEAN_DRIFT_FORCES_INF
        USE COMPUTE_MEAN_DRIFT_FORCES_FIN
        USE COMPUTE_YAW_MOMENT_INF
        USE COMPUTE_YAW_MOMENT_FIN
        USE SOLVE_BEM_HIGH
        USE SOLVE_BEM_DIPOLES_THIN
        USE SOLVE_BEM_IRREGULAR
        USE COMMON_TYPE
        !
        IMPLICIT NONE
        !       Inputs/outputs
        INTEGER :: ProblemNumber, Nintegration, NpanIsy
        REAL,DIMENSION(:,:) :: NDS
        COMPLEX,DIMENSION(:,:),ALLOCATABLE :: Force
        COMPLEX,DIMENSION(:,:),ALLOCATABLE :: NormalVelocity
        TYPE(TID) :: ID
        REAL :: Period
        COMPLEX,DIMENSION(:),ALLOCATABLE :: NVEL
        COMPLEX :: HKochin(:), PRESSURE(:)
        COMPLEX :: PHI(:)
        INTEGER :: Switch_Kochin, Switch_DRIFT_FORCES, Switch_YAW_MOMENT
        INTEGER :: NTheta
        REAL,DIMENSION(*) :: THeta
        INTEGER :: Switch_FS, higher_panel_method, num_panel_higher_order, b_spline_order ! description is clear
        INTEGER, VALUE :: ThreadId
        TYPE(TMesh) :: MeshFS, Mesh ! Mesh to get points and input Mesh respectively
        INTEGER :: Switch_potential, Switch_FastInfluence
        INTEGER:: is_wave, rad_case_idx ! is_wave problem, radiation case index for this problem respectively
        REAL:: rad_case(:, :) ! radiation case
        !       For solver (DIRECT or GMRES)
        INTEGER :: NEXP
        !       Locals
        INTEGER :: i,j
        REAL :: PI,W
        COMPLEX,PARAMETER :: II=CMPLX(0.,1.)
        REAL :: kwave
        REAL :: AKH,AMH, AM0
        REAL :: beta ! The wave incident angle
        COMPLEX,DIMENSION(1+MeshFS%Npoints) :: ETA
        REAL :: out_potential(:)

        REAL, DIMENSION(:, :) :: drift_forces
        REAL , DIMENSION(:):: yaw_moment
        TYPE(ParamsCommonInf) :: param
        TYPE(ParamsCommonInf) :: params
        TYPE(TEnvironment) :: Environment
        INTEGER, DIMENSION(:):: is_thin_body  ! whether a panel is a thin body or not
        INTEGER, DIMENSION(:):: is_interior_domain ! Array indicating whether or not a panel is in the body
                                                     ! or the interior free surface domain. 1 to indicate the interior, 0 for the body.
        INTEGER :: use_dipoles_implementation

        INTEGER :: remove_irregular_frequencies

        TYPE(TempVar) :: SolverVar
        CALL ALLOCATE_TEMP(SolverVar)

        ThreadId = ThreadId + 200



        !
        PI=4.*ATAN(1.)
        !       Define period
        SolverVar%T=Period
        W=2.*PI/SolverVar%T
        kwave=w*w/g
        SolverVar%ProblemNumber = ProblemNumber
        SolverVar%Switch_FastInfluence = Switch_FastInfluence

        ALLOCATE(NVEL(NFA*2**NSYMY))
        DO j=1,NpanIsy
            NVEL(j)=NormalVelocity(j,ProblemNumber)
        END DO


        params%SolverVar = SolverVar
        params%w = W
        Environment%Depth  = Depth
        Environment%Xeff = Xeff
        Environment%Yeff = Yeff

        params%Environment = Environment
        params%Mesh = Mesh
        params%beta = beta
        params%is_wave = is_wave
        params%rad_case_idx = rad_case_idx
        params%rad_case = rad_case
        params%NVEL = NVEL
        params%n_integration = Nintegration
        params%NpanIsy = NpanIsy

        params%NSPLIN = num_panel_higher_order
        params%KSPLIN = b_spline_order
        params%is_infinite = 0
        params%is_thin_body = is_thin_body
        params%is_interior_domain = is_interior_domain


        IF(higher_panel_method ==  1 .and. NSYMY.EQ.1) THEN
            write(*,*) 'Error: Symetric body computation with high panel method is not supported.'
            write(*,*)' Error: Set the symetric settings to 0 in order to run your problem by ignoring the symetry'
            STOP
        END IF

        IF(use_dipoles_implementation ==  1 .and. NSYMY.EQ.1) THEN
            write(*,*) 'Error: Symetric body computation with thin dipoles method is not supported.'
            write(*,*)' Error: Set the symetric settings to 0 in order to run your problem by ignoring the symetry'
            STOP
        END IF

        !       Select calculation in function of water depth
        IF ((Depth.EQ.0.).OR.(kwave*Depth.GE.20.)) THEN
            !           Calculate wave number
            kwave=w*w/g
            AMH=kwave*Depth
            AM0=W**2/G
            params%kwave = kwave
            params%AMH = AMH
            params%AKH = 0
            params%is_infinite = 1
            params%AM0=AM0



            IF (NSYMY.EQ.1 .and. Switch_FastInfluence == 1) THEN
                write(*,*) 'Error: Symetric body computation with ode influence coefficients is not supported.'
                write(*,*)' Error: Set the symetric settings to 0 in order to run your problem by ignoring the symetry'
                STOP
            END IF


            IF(use_dipoles_implementation ==  1) THEN

                CALL SOLVE_POTENTIAL_DIPOLES_THIN(params, PRESSURE, PHI, MeshFS, Switch_FS, Switch_Potential, NDS, &
                & Force(:, ProblemNumber),  out_potential)

            ELSEIF(higher_panel_method ==  1) THEN

                CALL SOLVE_POTENTIAL_HIGH(params, PRESSURE, PHI, MeshFS, Switch_FS, Switch_Potential, NDS, &
                & Force(:, ProblemNumber),  out_potential)

            ELSE IF(remove_irregular_frequencies == 1) THEN

                CALL SOLVE_POTENTIAL_IRREGULAR(params, PRESSURE, PHI, MeshFS, Switch_FS, Switch_Potential, NDS, &
                & Force(:, ProblemNumber),  out_potential)

            ELSE

                !           Solve with direct method ?
                IF (Indiq_solver == 0) CALL SOLVE_POTENTIAL_INFD_DIRECT(NVEL, SolverVar)
                !           Solve using GMRES ?
                IF (Indiq_solver == 1) THEN
                    CALL SOLVE_POTENTIAL_INFD_ITERATIVE(NVEL, SolverVar)
                END IF

            END IF



        ELSE
            !           Calculate wave number
            AKH=w*w*Depth/G                                                 
            AMH=X0(AKH)
            kwave=AMH/Depth
            AM0 = AMH/Depth
            AKH=AMH*TANH(AMH)

            params%kwave = kwave
            params%AMH = AMH
            params%AKH = AKH
            params%AM0=AM0


            IF(use_dipoles_implementation ==  1) THEN

                CALL SOLVE_POTENTIAL_DIPOLES_THIN(params, PRESSURE, PHI, MeshFS, Switch_FS, Switch_Potential, NDS, &
                & Force(:, ProblemNumber),  out_potential)

            ELSEIF(higher_panel_method ==  1) THEN
                CALL SOLVE_POTENTIAL_HIGH(params, PRESSURE, PHI, MeshFS, Switch_FS, Switch_Potential, NDS, &
                & Force(:, ProblemNumber), out_potential)

            ELSE IF(remove_irregular_frequencies == 1) THEN

                CALL SOLVE_POTENTIAL_IRREGULAR(params, PRESSURE, PHI, MeshFS, Switch_FS, Switch_Potential, NDS, &
                & Force(:, ProblemNumber),  out_potential)

            ELSE

                !           Solve with direct method ?
                IF (Indiq_solver == 0) CALL SOLVE_POTENTIAL_FD_DIRECT(NVEL,AMH,NEXP, SolverVar)
                !           Solve using GMRES ?
                IF (Indiq_solver == 1) THEN
                    CALL SOLVE_POTENTIAL_FD_ITERATIVE(NVEL,AMH,NEXP, SolverVar)
                    !STOP
                END IF

            END IF

        END IF

        IF(higher_panel_method ==  1 .or. use_dipoles_implementation ==1 .or.  remove_irregular_frequencies==1) THEN
            CALL REALEASE()
            RETURN
        END IF

        !       Assemble pressure
        DO i=1,NFA
            PRESSURE(i)=RHO*II*W*SolverVar%ZPB(i) !*AIRE(i)
            IF (NSYMY.EQ.1) THEN
                PRESSURE(i+NFA)=RHO*II*W*SolverVar%ZPS(i) !*AIRE(i)
            END IF
        END DO

        !       Save free surface elevation
        IF (Switch_FS.EQ.1) THEN
            DO j=1,MeshFS%Npoints
                CALL COMPUTE_POTENTIAL_DOMAIN(PHI(j),MeshFS%X(1,j),MeshFS%X(2,j),0.,kwave,AMH,NEXP, SolverVar)
                ETA(j)=II*W/G*PHI(j)
            END DO
        END IF

        !       Save output
        IF (Switch_Potential.EQ.1) THEN
            CALL SAVE_POTENTIAL(PRESSURE, out_potential)
        END IF

        !       Calculate force coefficients
        DO i=1,Nintegration
            DO j=1,NpanIsy
                Force(i,ProblemNumber)=Force(i,ProblemNumber)-PRESSURE(j)*NDS(i,j)
            END DO
        END DO



        !       Compute Kochin Functions
        IF (Switch_Kochin.EQ.1) THEN
            DO j=1,NTheta
                CALL COMPUTE_KOCHIN(kwave,Theta(j),HKochin(j), SolverVar)
            END DO
        END IF

        param%SolverVar = SolverVar
        param%w = W
        param%kwave = kwave
        param%AMH = AMH
        param%AKH = AKH
        param%AM0 = AM0



        ! Computing drift forces
        IF (Switch_DRIFT_FORCES .EQ. 1) THEN

            IF ((Depth.EQ.0.).OR.(kwave*Depth.GE.20.)) THEN

                DO j=1,NTheta

                param%beta = Theta(j)

                    CALL compute_mean_drift_forces_infinite(param, drift_forces(j, :))

               END DO

            ELSE
                DO j=1,NTheta

                param%beta = Theta(j)
                CALL compute_mean_drift_forces_finite(param, drift_forces(j, :))

                END DO

            ENDIF


        END IF

        ! Computing yaw_moment
        IF (Switch_YAW_MOMENT .EQ. 1) THEN

            IF ((Depth.EQ.0.).OR.(kwave*Depth.GE.20.)) THEN

            DO j=1,NTheta

                param%beta = Theta(j)

             CALL compute_yaw_moment_infinite(param, yaw_moment(j))

             END DO

            ELSE

               DO j=1,NTheta

                param%beta = Theta(j)
                CALL compute_yaw_moment_finite(param, yaw_moment(j))

                END DO

            ENDIF


        END IF




        CALL REALEASE()

        CONTAINS

        SUBROUTINE REALEASE()

            CALL DEALLOCATE_TEMP(SolverVar)
            DEALLOCATE(NVEL)

        END SUBROUTINE
    !
    END SUBROUTINE SOLVE_BVP
    !---------------------------------------------------------
!----------------------------------------------------------------
!    SUBROUTINE CREK(NPINTE)
!
!        USE COM_VAR
!        IMPLICIT NONE
!        INTEGER:: I,J,JZ,IR,NPINTE
!        REAL:: AKZ,AKR
!        JZ=46
!        IR=328
!        DO J=1,JZ
!            AKZ=AMIN1(10**(J/5.-6),10**(J/8.-4.5),16.)
!            XZ(J)=-AKZ     
!        END DO
!        XR(1)=0.
!        DO I=2,IR
!            IF(I.LT.40)THEN
!                AKR=AMIN1(10**((I-1.)/5-6),4./3.+ABS((I-32.)/3.))
!            ELSE
!                AKR=4./3.+ABS((I-32.)/3.)
!            ENDIF
!            XR(I)=AKR
!        END DO
!        DO J=1,JZ
!            DO I=1,IR
!                CALL VNS(NPINTE,XZ(J),XR(I),I,J)
!            END DO
!        END DO
!
!        RETURN
!
!    END SUBROUTINE
!---------------------------------------------------------                                                                    
!    SUBROUTINE VNS(NPINTE,AKZ,AKR,I,J)    
!
!        USE COM_VAR
!        USE ELEMENTARY_FNS       
!   
!        INTEGER:: I,J,NPINTE,IT
!        REAL:: AKZ,AKR,PI,CT,ST,CQIT,TETA
!        REAL:: QQT(NPINTE),CQT(NPINTE)
!        REAL:: FD1JX,FD1JZ,FD2JX,FD2JZ
!        COMPLEX:: IM,C1,C2,ZIK,GZ,CEX     
!        IM=(0.,1.)                                                                
!        PI=4.*ATAN(1.)
!        CALL COFINT(NPINTE,CQT,QQT)         
!        FD1JX=0.                                                              
!        FD1JZ=0.                                                              
!        FD2JX=0.                                                              
!        FD2JZ=0.                                                              
!        DO IT=1,NPINTE                                           
!            TETA=QQT(IT)                                                           !
!            CQIT=CQT(IT)
!            CT=COS(TETA)
!            ST=SIN(TETA)
!            ZIK=AKZ+IM*AKR*CT
!            IF(REAL(ZIK)+30.)2,2,1
!              2 CEX=(0.,0.)
!            GOTO 3
!              1 CEX=CEXP(ZIK)                                               
!              3 GZ=GG(ZIK,CEX)                      
!            C1=CQIT*(GZ-1./ZIK)
!            C2=CQIT*CEX
!            FD1JX=FD1JX+CT*AIMAG(C1)
!            FD1JZ=FD1JZ+REAL(C1)
!            FD2JX=FD2JX+CT*AIMAG(C2)
!            FD2JZ=FD2JZ+REAL(C2)
!        END DO                                                                
!        APD1X(I,J)=FD1JX                                               
!        APD1Z(I,J)=FD1JZ                                               
!        APD2X(I,J)=FD2JX                                               
!        APD2Z(I,J)=FD2JZ                                              
!        RETURN 
!!                                                                   
!    END SUBROUTINE  
!-------------------------------------------------------------------
!    SUBROUTINE COFINT(NPINTE,CQT,QQT)   
!
!    USE COM_VAR
!
!    INTEGER :: J,NPINTE
!    REAL:: PI,QQT(NPINTE),CQT(NPINTE)
!    PI=4.*ATAN(1.)
!      DO 160 J=1,NPINTE   
!      QQT(J)=-PI/2.+(J-1.)/(NPINTE-1.)*PI
!      IF(J-1)161,161,162
!  161 CQT(J)=PI/(3.*(NPINTE-1.))
!      GOTO 160
!  162 IF(J-NPINTE)163,161,161
!  163 IF(MOD(J,2))164,165,164
!  164 CQT(J)=2./(3.*(NPINTE-1.))*PI
!      GOTO 160
!  165 CQT(J)=4./(3.*(NPINTE-1.))*PI
!  160 CONTINUE
!    RETURN                                                                    
!
!    END SUBROUTINE   
!------------------------------------------------------------------
END MODULE
