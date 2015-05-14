!--------------------------------------------------------------------------------------
!
!Copyright (C) 2014 TopCoder Inc., All Rights Reserved.
!
!--------------------------------------------------------------------------------------

!   This module solves the finite linear BVP (Boundary Value Problem) for the potential.
!   It used an iterative solver for system of linear equations called GMRES
!
! Changes in version 1.2 (Implementation of Higher Order Panel Methods)
!       Added COMMON_TYPE module as dependency

!   @author yedtoss
!   @version 1.2

MODULE SOLVE_BEM_FD_ITERATIVE

    USE COMMON_TYPE
    USE COM_VAR
    USE ELEMENTARY_FNS
    USE COMPUTE_GREEN_FD
    USE M_SOLVER
    USE GMRES
    IMPLICIT NONE

CONTAINS
    !--------------------------------------------------------------------------!
    SUBROUTINE SOLVE_POTENTIAL_FD_ITERATIVE(NVEL,AMH,NEXP, SolverVar)

        !In this subroutine the linear system Ax=b is constructed

        !    LOGICAL:: RHS
        COMPLEX,DIMENSION(*) :: NVEL
        INTEGER:: BX,ISYM,NJJ
        INTEGER:: I,J,ISP,IFP
        INTEGER:: NP1,I1,JJ, N
        REAL:: GM, BJ, DIJ, AKK
        INTEGER:: NEXP
        REAL:: W
        REAL:: PCOS,PSIN
        REAL:: AM0,AKH,AMH,PI,DPI,ZERO
        COMPLEX:: ZOL(IMX,2)

        integer :: lwork(2)
        integer :: revcom, colx, coly, colz, nbscal
        integer :: irc(5,2), icntl(8), info(3,2)
        TYPE(GMRESVAR):: gmres_env(2)

        integer :: matvec, precondLeft, precondRight, dotProd
        parameter (matvec=1, precondLeft=2, precondRight=3, dotProd=4)

        complex, dimension(:), allocatable ::  work
        real ::  cntl(5), rinfo(2,2)
        complex :: ONE, ZZERO

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

        lwork = (restartParam+10)**2 + (restartParam+10)*(IMX+5) + 5*IMX + 1
        ALLOCATE(work(lwork(1)))
        ONE = CMPLX(1., 0.)
        ZZERO = CMPLX(0., 0.)



        NJJ=NSYMY+1
        PI=4.*ATAN(1.)
        DPI=2.*PI
        W=DPI/T
        ZERO=0.0
        AKH=W**2*Depth/G
        !       AMH=X0(AKH)
        AKH=AMH*TANH(AMH)
        AM0=AMH/Depth


        IF(ABS(W).LT.1.E-4)THEN
            WRITE(*,*)'ABS(WR)  = ',ABS(W),' < 1.E-4'
            STOP
        ENDIF

        !--------------Initilizations---------------
        CQ=0.0
        QQ=0.0
        AMBDA=0.0
        AR=0.0
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
            CALL CINT_FD(AKK,N)
            !N=NQ Bug Fixes
            NQ = N
            CALL LISC(AKH,AMH,NEXP, SolverVar)

            IF(ProblemSavedAt(SolverVar%ProblemNumber) == -1 ) THEN
                DO I=1,31
                    ARCache(ProblemToSaveLocation(SolverVar%ProblemNumber), I)  = AR(I)
                    AMBDACache(ProblemToSaveLocation(SolverVar%ProblemNumber), I) =  AMBDA(I)

                END DO
                NEXPCache(ProblemToSaveLocation(SolverVar%ProblemNumber)) = NEXP

            END IF

        ELSE

            DO I=1,31
                AR(I) = ARCache(ProblemSavedAt(SolverVar%ProblemNumber), I)
                AMBDA(I) = AMBDACache(ProblemSavedAt(SolverVar%ProblemNumber), I)
            END DO

            NEXP = NEXPCache(ProblemSavedAt(SolverVar%ProblemNumber))

        END IF

        !Construction of the influence matrix
        DO ISYM=1,NJJ
            BX=(-1)**(ISYM+1)
            IF(ProblemSavedAt(SolverVar%ProblemNumber) < 0 ) THEN
                DO ISP=1,IMX
                    IF(ISYM.EQ.1)THEN
                        DO IFP=1,IMX
                            call VAVFD(1, 2,XG(ISP),YG(ISP),ZG(ISP),ISP,IFP, SolverVar)  !1/r+1/r1
                            call VNSFD(1, AM0,AMH,NEXP,ISP,IFP,XG(ISP),YG(ISP),ZG(ISP),  SolverVar)

                            PCOS=VSXP1*XN(ISP)+VSYP1*YN(ISP)+VSZP1*ZN(ISP)
                            PSIN=VSXP2*XN(ISP)+VSYP2*YN(ISP)+VSZP2*ZN(ISP)
                            ZIJ(ISP,IFP)=CMPLX(PCOS,PSIN)

                        END DO

                    ELSE
                        DO IFP=1,IMX
                            call VAVFD(1, 2,XG(ISP),YG(ISP),ZG(ISP),ISP,IFP, SolverVar)

                            call VNSFD(1, AM0,AMH,NEXP,ISP,IFP,XG(ISP),YG(ISP),ZG(ISP), SolverVar)

                            PCOS=VSXM1*XN(ISP)+VSYM1*YN(ISP)+VSZM1*ZN(ISP)
                            PSIN=VSXM2*XN(ISP)+VSYM2*YN(ISP)+VSZM2*ZN(ISP)
                            ZIJ(ISP,IFP)=CMPLX(PCOS,PSIN)
                        END DO
                    ENDIF
                END DO

            END IF


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

            !******************************************************
            !* Initialize the control parameters to default value
            !******************************************************

            call init_cgmres(icntl,cntl, gmres_env(ISYM))

            !************************
            !c Tune some parameters
            !************************

            ! Save the convergence history on standard output
            !icntl(3) = 30
            icntl(2) = 0
            ! Maximum number of iterations
            icntl(7) = maxIterations

            ! Tolerance
            cntl(1) = TOL_GMRES

            ! preconditioner location
            icntl(4) = 0
            ! orthogonalization scheme
            icntl(5)=0
            ! initial guess
            icntl(6) = 0
            ! residual calculation strategy at restart
            icntl(8) = 1
            !* Initialise the right hand side

            work = CMPLX(0., 0.)

            work(1:IMX) = ONE

            DO I=1,IMX
                IF (NSYMY.EQ.1) THEN
                    work(IMX+I)=(NVEL(I)+BX*NVEL(I+NFA))*0.5
                ELSE
                    work(IMX+I)=NVEL(I)
                END IF
            ENDDO

            DO

                call drive_cgmres(IMX,IMX,restartParam,lwork(ISYM),work, &
                    irc(:, ISYM),icntl,cntl,info(:, ISYM),rinfo(:, ISYM), gmres_env(ISYM))
                revcom = irc(1, ISYM)
                colx   = irc(2, ISYM)
                coly   = irc(3, ISYM)
                colz   = irc(4, ISYM)
                nbscal = irc(5, ISYM)

                if (revcom == matvec) then
                    ! perform the matrix vector product
                    !        work(colz) <-- A * work(colx)
                    call cgemv('N',IMX,IMX,ONE,ZIJ,IMX,work(colx),1, &
                        ZZERO,work(colz),1)

                else if (revcom == precondLeft) then
                    ! perform the left preconditioning
                    !         work(colz) <-- M^{-1} * work(colx)
                    call ccopy(IMX,work(colx),1,work(colz),1)
                    call ctrsm('L','L','N','N',IMX,1,ONE,ZIJ,IMX,work(colz),IMX)

                else if (revcom == precondRight) then
                    ! perform the right preconditioning
                    call ccopy(IMX,work(colx),1,work(colz),1)
                    call ctrsm('L','U','N','N',IMX,1,ONE,ZIJ,IMX,work(colz),IMX)

                else if (revcom == dotProd) then
                    !      perform the scalar product
                    !      work(colz) <-- work(colx) work(coly)

                    call cgemv('C',IMX,nbscal,ONE,work(colx),IMX, &
                        work(coly),1,ZZERO,work(colz),1)

                ELSE
                    exit
                endif

            END DO


            DO I=1,IMX
                ZOL(I,(ISYM-1)+1)=work(I)
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
            DO I=1,IMX
                DO J=1,IMX
                    call VAVFD(1, 2,XG(I),YG(I),ZG(I),I,J, SolverVar)

                    call VNSFD(1, AM0,AMH,NEXP,I,J,XG(I),YG(I),ZG(I), SolverVar)

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

        ELSE

            N = ProblemSavedAt(SolverVar%ProblemNumber)



            DO I=1,IMX
                DO J=1,IMX

                    ZPB(I)=ZPB(I)+0.5*(ZIGB(J)*CMPLX(SP1Cache(N,I,J)+SM1Cache(N,I,J),SP2Cache(N,I,J)+SM2Cache(N,I,J))&
                        +ZIGS(J)*CMPLX(SP1Cache(N,I,J)-SM1Cache(N,I,J),SP2Cache(N,I,J)-SM2Cache(N,I,J)))
                    ZPS(I)=ZPS(I)+0.5*(ZIGS(J)*CMPLX(SP1Cache(N,I,J)+SM1Cache(N,I,J),SP2Cache(N,I,J)+SM2Cache(N,I,J))&
                        +ZIGB(J)*CMPLX(SP1Cache(N,I,J)-SM1Cache(N,I,J),SP2Cache(N,I,J)-SM2Cache(N,I,J)))
                END DO

            END DO

        END IF

    
    END SUBROUTINE
!----------------------------------------------------------------------

END MODULE

!*  some control parametrs definition for GMRES
!* icntl    (input) INTEGER array. length 6
!*            icntl(1) : stdout for error messages
!*            icntl(2) : stdout for warnings
!*            icntl(3) : stdout for convergence history
!*            icntl(4) : 0 - no preconditioning
!*                       1 - left preconditioning
!*                       2 - right preconditioning
!*                       3 - double side preconditioning
!*                       4 - error, default set in Init
!*            icntl(5) : 0 - modified Gram-Schmidt
!*                       1 - iterative modified Gram-Schmidt
!*                       2 - classical Gram-Schmidt
!*                       3 - iterative classical Gram-Schmidt
!*            icntl(6) : 0 - default initial guess x_0 = 0 (to be set)
!*                       1 - user supplied initial guess
!*            icntl(7) : maximum number of iterations
!*            icntl(8) : 1 - default compute the true residual at each restart
!*                       0 - use recurence formaula at restart
!*
!* cntl     (input) real array, length 5
!*            cntl(1) : tolerance for convergence
!*            cntl(2) : scaling factor for normwise perturbation on A
!*            cntl(3) : scaling factor for normwise perturbation on b
!*            cntl(4) : scaling factor for normwise perturbation on the
!*                      preconditioned matrix
!*            cntl(5) : scaling factor for normwise perturbation on 
!*                      preconditioned right hand side
