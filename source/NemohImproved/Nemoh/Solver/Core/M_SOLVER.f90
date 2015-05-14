!--------------------------------------------------------------------------------------
!
!Copyright (C) 2014 TopCoder Inc., All Rights Reserved.
!
!--------------------------------------------------------------------------------------

!--------------------------------------------------------------------------------------
!
!   Copyright 2014 Ecole Centrale de Nantes, 1 rue de la Noë, 44300 Nantes, France
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
!   - P. Guével
!   - J.C. Daubisse
!   - J. Singh  
!
!--------------------------------------------------------------------------------------


!   This module contains direct solver to find a solution for a system of complex linear equations
!   Currently, only GAUSS method is implemented.
!
!   @author TCSASSEMBLER
!   @version 1.1
MODULE M_SOLVER

    IMPLICIT NONE
 
CONTAINS
    !---------------------------------------------------------------------------!
    SUBROUTINE GAUSSZ(A,NMAX,N,M)
        INTEGER:: N,M,NMAX,I,J,K,L,IL
        !COMPLEX:: A(NMAX,1),C,P Bug Fixes
        COMPLEX:: A(N,M),C,P
        REAL:: EPS
        CHARACTER(LEN=*), PARAMETER  :: FMT = "(5X,'PIVOT INFERIEUR A ',1P,E16.6)"
        
        EPS=1.E-20
        DO J=1,N-1
            K=J
            DO I=J+1,N
                IF(ABS(A(K,J))-ABS(A(I,J)) < 0) THEN
                    K=I
                END IF

            END DO

            IF(K-J /= 0) THEN

                DO L=J,M
                    C=A(J,L)
                    A(J,L)=A(K,L)
                    A(K,L)=C

                END DO

            END IF

            IF(ABS(A(J,J))-EPS <= 0) THEN
                WRITE(*,FMT)EPS
                STOP

            END IF

            DO K=J+1,M
                P=A(J,K)/A(J,J)
                DO I=J+1,N
                    A(I,K)=A(I,K)-A(I,J)*P
                END DO
            END DO
        END DO

        IF(ABS(A(N,N))-EPS <= 0) THEN
            WRITE(*,FMT)EPS
            STOP
        END IF


        DO IL=N+1,M
            DO J=N,1,-1
                A(J,IL)=A(J,IL)/A(J,J)
                DO I=1,J-1
                    A(I,IL)=A(I,IL)-A(I,J)*A(J,IL)
                END DO
            END DO
        END DO
    END SUBROUTINE

END MODULE 
