!--------------------------------------------------------------------------------------
!
!Copyright (C) 2014 TopCoder Inc., All Rights Reserved.
!
!--------------------------------------------------------------------------------------

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
!   - J. Singh
!   - P. Guevel
!   - J.C. Daubisse 
!
!--------------------------------------------------------------------------------------


!   This module solves the infinite linear BVP (Boundary Value Problem) for the potential.
!   It used an direct solver for system of linear equations using an LU Factorization.
!   This direct solver is thus the Gauss method
!
! Changes in version 1.2 (Code Acceleration of the Calculation of Influence Coefficients of Nemoh)
!       Allow fast code for influence coefficients to be run according to a switch
!
! Changes in version 1.3 (Implementation of Higher Order Panel Methods)
!       Added COMMON_TYPE module as dependency

!   @author yedtoss
!   @version 1.3
MODULE SOLVE_BEM_INFD_DIRECT

    USE COMMON_TYPE
    USE COM_VAR
    USE COMPUTE_GREEN_INFD
    USE ELEMENTARY_FNS
    USE M_SOLVER
    USE COMPUTE_INFLUENCE_ODE
    IMPLICIT NONE

CONTAINS
    !--------------------------------------------------------------------------!
    SUBROUTINE SOLVE_POTENTIAL_INFD_DIRECT(NVEL,SolverVar)
        !In this subroutine the linear system Ax=b is constructed

        COMPLEX,DIMENSION(*) :: NVEL
        INTEGER :: BX,ISYM,NJJ
        INTEGER ::I,J,ISP,IFP
        INTEGER:: NP1,I1,JJ, N
        REAL:: GM, BJ, DIJ, AKK
        REAL:: W,tdepth
        REAL:: PCOS,PSIN
        REAL:: AM0,PI,DPI,ZERO
        COMPLEX:: ZOL(IMX,2)

        TYPE(TempVar), TARGET :: SolverVar
        REAL, POINTER :: T
        COMPLEX, DIMENSION(:), POINTER :: ZPB,ZPS
        COMPLEX, DIMENSION(:), POINTER :: ZIGB,ZIGS
        COMPLEX, DIMENSION(:, :), POINTER :: ZIJ
        REAL, POINTER :: FSP,FSM,VSXP,VSYP,VSZP,VSXM,VSYM,VSZM
        REAL, POINTER :: SP1,SM1,SP2,SM2
        REAL, POINTER :: VSXP1,VSXP2,VSYP1,VSYP2,VSZP1,VSZP2
        REAL, POINTER :: VSXM1,VSXM2,VSYM1,VSYM2,VSZM1,VSZM2
        INTEGER, POINTER:: NQ
        REAL, POINTER:: CQ(:),QQ(:),AMBDA(:),AR(:)
        T => SolverVar%T
        ZPB => SolverVar%ZPB
        ZPS => SolverVar%ZPS
        ZIGB => SolverVar%ZIGB
        ZIGS => SolverVar%ZIGS
        ZIJ => SolverVar%ZIJ
        FSP => SolverVar%FSP
        FSM => SolverVar%FSM
        VSXP => SolverVar%VSXP
        VSYP => SolverVar%VSYP
        VSZP => SolverVar%VSZP
        VSXM => SolverVar%VSXM
        VSYM => SolverVar%VSYM
        VSZM => SolverVar%VSZM
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
        NQ => SolverVar%NQ
        CQ => SolverVar%CQ
        QQ => SolverVar%QQ
        AMBDA => SolverVar%AMBDA
        AR => SolverVar%AR


        NJJ=NSYMY+1
        PI=4.*ATAN(1.)
        DPI=2.*PI
        W=DPI/T
        ZERO=0.0
        AM0=W**2/G
        tdepth=1.E+20
        !--------------Initilizations---------------
        CQ=0.0
        QQ=0.0
        VSXP=0.
        VSXM=0.
        VSYP=0.
        VSYM=0.
        VSZP=0.
        VSZM=0.
        ZOL=CMPLX(0.,0.)
        ZIJ=CMPLX(0.,0.)
        !--------------------------------------------

        IF(ProblemSavedAt(SolverVar%ProblemNumber) < 0 ) THEN
            GM=0.
            NP1=NP-1
            DO I=1,NP1
                I1=I+1
                DO JJ=1,NJJ
                    BJ=(-1.)**(JJ+1)
                    DO J=I1,NP
                        DIJ=SQRT((X(J)-X(I))**2+(Y(I)-Y(J)*BJ)**2)
                        GM=MAX(DIJ,GM)
                    END DO
                END DO
            END DO
            AKK=AM0*GM
            CALL CINT_INFD(AKK,N, SolverVar)
            !N=NQ Bug Fixes
            NQ = N
            IF(ProblemSavedAt(SolverVar%ProblemNumber) == -1 ) THEN
                DO I=1,101
                    CQCache(ProblemToSaveLocation(SolverVar%ProblemNumber), I)  = CQ(I)
                    QQCache(ProblemToSaveLocation(SolverVar%ProblemNumber), I) =  QQ(I)
                END DO

            END IF

        ELSE

            DO I=1,101
                CQ(I) = CQCache(ProblemSavedAt(SolverVar%ProblemNumber), I)
                QQ(I) = QQCache(ProblemSavedAt(SolverVar%ProblemNumber), I)
            END DO

        END IF

        !Construction of the influence matrix
        DO ISYM=1,NJJ
            BX=(-1)**(ISYM+1)

            IF(ProblemSavedAt(SolverVar%ProblemNumber) < 0 ) THEN

                IF (SolverVar%Switch_FastInfluence ==1) THEN
                    CALL compute_influence_infinite_dipoles(W, ZIJ) ! Call the subroutine to compute influence

                    is_rankine_dipoles_computed = 1

                ELSE
                    DO ISP=1,IMX



                        IF(ISYM.EQ.1)THEN
                            DO IFP=1,IMX

                                call VAVINFD(1, 1,XG(ISP),YG(ISP),ZG(ISP),ISP,IFP, SolverVar)  !1/r+1/r1

                                call VNSINFD(1, 1,ISP,IFP,XG(ISP),YG(ISP),ZG(ISP), SolverVar)

                                PCOS=VSXP1*XN(ISP)+VSYP1*YN(ISP)+VSZP1*ZN(ISP)
                                PSIN=VSXP2*XN(ISP)+VSYP2*YN(ISP)+VSZP2*ZN(ISP)
                                ZIJ(ISP,IFP)=CMPLX(PCOS,PSIN)
                            END DO
                        ELSE
                            DO IFP=1,IMX
                                call VAVINFD(1, 1,XG(ISP),YG(ISP),ZG(ISP),ISP,IFP, SolverVar)

                                call VNSINFD(1, 1,ISP,IFP,XG(ISP),YG(ISP),ZG(ISP), SolverVar)

                                PCOS=VSXM1*XN(ISP)+VSYM1*YN(ISP)+VSZM1*ZN(ISP)
                                PSIN=VSXM2*XN(ISP)+VSYM2*YN(ISP)+VSZM2*ZN(ISP)
                                ZIJ(ISP,IFP)=CMPLX(PCOS,PSIN)
                            END DO
                        ENDIF

                    END DO

               END IF

            END IF

            DO I=1,IMX
                IF (NSYMY.EQ.1) THEN
                    ZIJ(I,IMX+1)=(NVEL(I)+BX*NVEL(I+NFA))*0.5
                ELSE
                    ZIJ(I,IMX+1)=NVEL(I)
                END IF
            ENDDO

            IF(ProblemSavedAt(SolverVar%ProblemNumber) == -1 ) THEN

                DO I=1,IMX
                    DO J=1,IMX

                        InfluenceMatrixCache(ProblemToSaveLocation(SolverVar%ProblemNumber), ISYM, I, J) = ZIJ(I,J)

                    END DO

                END DO

            END IF


            IF(ProblemSavedAt(SolverVar%ProblemNumber) >= 0 ) THEN

                DO I=1,IMX
                    DO J=1,IMX

                        ZIJ(I,J) = InfluenceMatrixCache(ProblemSavedAt(SolverVar%ProblemNumber), ISYM, I, J)

                    END DO

                END DO

            END IF
 
            !------------------------------------------------!
            CALL GAUSSZ(ZIJ,NFA,IMX,IMX+1)
            !------------------------------------------------!

            DO I=1,IMX
                ZOL(I,(ISYM-1)+1)=ZIJ(I,IMX+1)
            END DO

        END DO

        ZIGB=(0.,0.)
        ZIGS=(0.,0.)

        DO I=1,IMX
            IF(NSYMY .EQ. 0)THEN
                ZIGB(I)=ZOL(I,1)    ! Factor 2 is removed in comparison with previous version of Nemoh because
                                    ! Normal velocity is halved in Nemoh (because of symmetry)
                ZIGS(I)=0.0
            ELSE
                ZIGB(I)=(ZOL(I,1)+ZOL(I,2))
                ZIGS(I)=(ZOL(I,1)-ZOL(I,2))
            ENDIF
        END DO

        !computation of potential phi=S*sigma on the boundary
        ZPB=(0.,0.)
        ZPS=(0.,0.)

        IF(ProblemSavedAt(SolverVar%ProblemNumber) < 0 ) THEN
            N = ProblemToSaveLocation(SolverVar%ProblemNumber)
            IF (SolverVar%Switch_FastInfluence ==1) THEN

                CALL compute_influence_infinite_sources(W, ZIJ) ! Call the subroutine to compute influence
                !!$OMP CRITICAL
                !We are safe without a critical or atomic operation. Although we are
                ! updating it's value, this will be it's final value. And If a thread is still seeing 0 instead
                ! of 1 after this update, it won't be a issue.  Any thread updating it is doing so with the same
                ! value
                is_rankine_sources_computed = 1
                !!$OMP END CRITICAL
                DO I=1,IMX
                    DO J=1,IMX
                        ZPS(I)= ZPS(I) -  ZIJ(I, J) * ZIGS(J)
                        ZPB(I)= ZPB(I) -  ZIJ(I, J) * ZIGB(J)
                        IF(ProblemSavedAt(SolverVar%ProblemNumber) == -1 ) THEN
                            SM1Cache(N, I, J) = REAL(ZIJ(I,J))
                            SM2Cache(N, I, J) = AIMAG(ZIJ(I,J))
                        END IF
                    END DO
                END DO

            ELSE
                DO I=1,IMX
                    DO J=1,IMX

                        call VAVINFD(1, 1,XG(I),YG(I),ZG(I),I,J, SolverVar)

                        call VNSINFD(1, 1,I,J,XG(I),YG(I),ZG(I), SolverVar)

                        IF(ProblemSavedAt(SolverVar%ProblemNumber) == -1 ) THEN
                            SM1Cache(N, I, J) = SM1
                            SM2Cache(N, I, J) = SM2
                            SP1Cache(N, I, J) = SP1
                            SP2Cache(N, I, J) = SP2
                        END IF



                        ZPB(I)=ZPB(I)+0.5*(ZIGB(J)*CMPLX(SP1+SM1,SP2+SM2)&
                            +ZIGS(J)*CMPLX(SP1-SM1,SP2-SM2))
                        ZPS(I)=ZPS(I)+0.5*(ZIGS(J)*CMPLX(SP1+SM1,SP2+SM2)&
                            +ZIGB(J)*CMPLX(SP1-SM1,SP2-SM2))
                    END DO

                END DO
           END IF

        ELSE

            N = ProblemSavedAt(SolverVar%ProblemNumber)


            IF (SolverVar%Switch_FastInfluence ==1) THEN
                DO I=1,IMX
                    DO J=1,IMX
                        ZPS(I)= ZPS(I) -  CMPLX(SM1Cache(N, I, J), SM2Cache(N, I, J)) * ZIGS(J)
                        ZPB(I)= ZPB(I) -  CMPLX(SM1Cache(N, I, J), SM2Cache(N, I, J)) * ZIGB(J)
                    END DO
                END DO
            ELSE

            DO I=1,IMX
                DO J=1,IMX

                    ZPB(I)=ZPB(I)+0.5*(ZIGB(J)*CMPLX(SP1Cache(N,I,J)+SM1Cache(N,I,J),SP2Cache(N,I,J)+SM2Cache(N,I,J))&
                        +ZIGS(J)*CMPLX(SP1Cache(N,I,J)-SM1Cache(N,I,J),SP2Cache(N,I,J)-SM2Cache(N,I,J)))
                    ZPS(I)=ZPS(I)+0.5*(ZIGS(J)*CMPLX(SP1Cache(N,I,J)+SM1Cache(N,I,J),SP2Cache(N,I,J)+SM2Cache(N,I,J))&
                        +ZIGB(J)*CMPLX(SP1Cache(N,I,J)-SM1Cache(N,I,J),SP2Cache(N,I,J)-SM2Cache(N,I,J)))
                END DO

            END DO

            ENDIF



        END IF
    !      do i=1,IMX
    !      print*,'ZPB,ZPS',I,ZPB(I),ZPS(I)
    !      enddo

    END SUBROUTINE

END MODULE
