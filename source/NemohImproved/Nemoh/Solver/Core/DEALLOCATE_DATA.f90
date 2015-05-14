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

!   This subroutine deallocates variables used to solve all problems
!
!   Changes in version 1.2 (Code Acceleration of the Calculation of Influence Coefficients of Nemoh)
!       Deallocate variables for caching rankine values.
!
!   @author yedtoss
!   @version 1.2
    SUBROUTINE DEALLOCATE_DATA

      USE COM_VAR

      DEALLOCATE(X,Y,Z,XG,YG,ZG)
      DEALLOCATE(AIRE,M1,M2,M3,M4)
      DEALLOCATE(XN,YN,ZN,DIST,TDIS)
      DEALLOCATE(XR, XZ, APD1X, APD1Z, APD2X, APD2Z)
      DEALLOCATE(rankine_sources_cache, rankine_dipoles_cache)
      !if(Indiq_solver .eq. 0)DEALLOCATE(ZIJ)
      !DEALLOCATE(ZPB,ZPS)
      !DEALLOCATE(ZIGB,ZIGS)
      
    END SUBROUTINE
