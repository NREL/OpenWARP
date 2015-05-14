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
!   - A. Babarit 
!
!--------------------------------------------------------------------------------------


!   This module handles the output of Nemoh computation.
!   It saves the potential on body surfaces for problem XX.
!
!    Changes in version 1.2:
!           Remove all subroutine for writing. Added subroutine SAVE_POTENTIAL
!
!   @author TCSASSEMBLER
!   @version 1.2
MODULE OUTPUT

CONTAINS


    SUBROUTINE SAVE_POTENTIAL(PRESSURE, out_potential)
    USE MIDENTIFICATION
    USE MMesh
    USE COM_VAR
    IMPLICIT NONE
    INTEGER:: i, idx
    COMPLEX :: PRESSURE(:)
    REAL :: out_potential(:)

    idx = 1

    DO i=1,NP
        out_potential(idx) =  X(i)
        out_potential(idx + 1) =  Y(i)
        out_potential(idx + 2) =  Z(i)
        out_potential(idx + 3) =  0.
        out_potential(idx + 4) =  0.
        idx = idx + 5
    END DO
    IF (NSYMY.EQ.1) THEN
        DO i=1,NP
            out_potential(idx) =  X(i)
            out_potential(idx + 1) =  -Y(i)
            out_potential(idx + 2) =  Z(i)
            out_potential(idx + 3) =  0.
            out_potential(idx + 4) =  0.
            idx = idx + 5
        END DO
    END IF
    DO i=1,NFA
        out_potential(idx) =  M1(i)
        out_potential(idx + 1) =  M2(i)
        out_potential(idx + 2) =  M3(i)
        out_potential(idx + 3) =  M4(i)
        idx = idx + 4
    END DO
    IF (NSYMY.EQ.1) THEN
        DO i=1,NFA
            out_potential(idx) =  M1(i)+NFA
            out_potential(idx + 1) =  M2(i)+NFA
            out_potential(idx + 2) =  M3(i)+NFA
            out_potential(idx + 3) =  M4(i)+NFA
            idx = idx + 4
        END DO
    END IF
    DO i=1,NFA
        out_potential(idx) =  XG(i)
        out_potential(idx + 1) =  YG(i)
        out_potential(idx + 2) =  ZG(i)
        out_potential(idx + 3) =  ABS(PRESSURE(i))
        out_potential(idx + 4) =  ATAN2(AIMAG(PRESSURE(i)),REAL(PRESSURE(i)))
        idx = idx + 5
    END DO
    IF (NSYMY.EQ.1) THEN
        DO i=1,NFA
            out_potential(idx) =  XG(i)
            out_potential(idx + 1) =  -YG(i)
            out_potential(idx + 2) =  ZG(i)
            out_potential(idx + 3) =  ABS(PRESSURE(i+NFA))
            out_potential(idx + 4) =  ATAN2(AIMAG(PRESSURE(i+NFA)),REAL(PRESSURE(i+NFA)))
            idx = idx + 5
        END DO
    END IF
    END SUBROUTINE
  
 END MODULE
