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


!   This module computes the potential for problem already solved.
!
! Changes in version 1.2 (Implementation of Higher Order Panel Methods)
!       Added COMMON_TYPE module as dependency

!   @author yedtoss
!   @version 1.2
MODULE POTENTIAL_DOMAIN

    USE COMMON_TYPE
    IMPLICIT NONE

CONTAINS


    SUBROUTINE COMPUTE_POTENTIAL_DOMAIN(PHI,XC,YC,ZC,AM0,AMH,NEXP, SolverVar)
        !
        USE MIDENTIFICATION
        USE COM_VAR
        USE ELEMENTARY_FNS
        USE COMPUTE_GREEN_FREESURFACE
        !

        !
        INTEGER::J,NEXP
        REAL::XC,YC,ZC,AM0,AMH
        REAL:: S1B,S2B,S1S,S2S,C,D,POT1,POT2
        REAL:: W,PI,DPI
        COMPLEX :: PHI

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
        !
        PI=4.*ATAN(1.)
        DPI=2.*PI
        W=DPI/T
        POT1=0.0
        POT2=0.0

        DO J=1,IMX
            SP1=0.
            SP2=0.
            SM1=0.
            SM2=0.
            IF((Depth .EQ. 0.) .OR. (AMH .GE. 20))  THEN
                CALL VVV(1, 1,J,XC,YC,ZC, SolverVar) ! finite part
                CALL VNV(1, J,XC,YC,ZC, SolverVar)  !infinite part
            ELSE
                CALL VVV(1, 2,J,XC,YC,ZC, SolverVar)
                CALL VNVF(1, AM0,AMH,NEXP,J,XC,YC,ZC, SolverVar)
            ENDIF
            S1B=REAL(ZIGB(J))
            S1S=REAL(ZIGS(J))
            S2B=AIMAG(ZIGB(J))
            S2S=AIMAG(ZIGS(J))
            IF(NSYMY.EQ.0)THEN
                C=S1B*SP1
                D=S2B*SP2
                POT1=POT1+C-D
                C=S1B*SP2
                D=S2B*SP1
                POT2=POT2+C+D
            ELSE
                C=S1B*(SP1+SM1)+S1S*(SP1-SM1)
                D=S2B*(SP2+SM2)+S2S*(SP2-SM2)
                POT1=POT1+(C-D)*0.5
                C=S1B*(SP2+SM2)+S1S*(SP2-SM2)
                D=S2B*(SP1+SM1)+S2S*(SP1-SM1)
                POT2=POT2+(C+D)*0.5
            ENDIF
        END DO

        PHI=CMPLX(POT1,POT2)
        RETURN
    !      PRE1=-POT2*W/G
    !      PRE2=POT1*W/G
    !      POMOD=SQRT(PRE1**2+PRE2**2)
    !
    END SUBROUTINE




END MODULE POTENTIAL_DOMAIN
