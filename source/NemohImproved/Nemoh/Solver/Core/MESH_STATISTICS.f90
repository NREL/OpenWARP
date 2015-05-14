!--------------------------------------------------------------------------------------
!
!Copyright (C) 2015 TopCoder Inc., All Rights Reserved.
!
!--------------------------------------------------------------------------------------

! Use to compute mesh properties such as displacement, buoyancy center, hydrostatic stiffness.
! It also makes estimate of masses and inertia matrix
! It only compute those statistics for a single body
!
!   @author TCSASSEMBLER
!   @version 1.0

MODULE MESH_STATISTICS

    IMPLICIT NONE

CONTAINS

    SUBROUTINE COMPUTE_MESH_STATISTICS(Nsym, tX, tY, Tcol, lambda, xG, yG, zG, Np, Nface, X, Y, Z, P, ct_RHO, ct_G, &
        & center_buoyancy, displacement, waterplane_area, stifness)
        USE HYDRO_NETTING


        REAL :: center_buoyancy(3)
        REAL :: displacement
        REAL:: waterplane_area
        REAL :: stifness(6, 6)

            !   Nombre de symétries du flotteur
        INTEGER :: Nsym

        !   Position du maillage dans le plan (x0y)
        REAL :: tX,tY
        !   Description du flotteur en faces
        INTEGER,PARAMETER :: nFacemx=3000
        !   Description de tout le flotteur
        INTEGER :: nFace                        ! Nombre reel de faces
        REAL,DIMENSION(4,3,nFace) :: Coin        ! Coins des faces
        !   Description de la partie immergee
        INTEGER :: nFacem                        ! Nombre reel de faces
        REAL,DIMENSION(4,3,nFace) :: Coinm    ! Coins des faces
        !   Parametres de reglage du maillage
        INTEGER  nmaillage
        REAL maille
        REAL Tcol
        REAL :: lambda                                  ! Facteur d echelle
        !   Maillage proprement dit
        INTEGER,PARAMETER :: NFMX=20000                    ! Nombre de facettes max
        INTEGER,PARAMETER :: NPMX=20000                    ! Nombre de points max
        INTEGER :: Nmailmx        ! Nombre de facettes du maillage std max
        !   Maillage du corps
        INTEGER :: NF,NP
        INTEGER,DIMENSION(4,Nface) :: Facette
        REAL,DIMENSION(Np) :: X,Y,Z
        !   Partie immergee du maillage
        INTEGER :: NFm,NPm
        INTEGER,DIMENSION(4,Nface) :: Facettem
        REAL,DIMENSION(Np) :: Xm,Ym,Zm
        !   Calcul hydrostatique
        REAL DEPLACEMENT,XF,YF,ZF,SF
        REAL,DIMENSION(6,6) :: KH
        REAL :: xG,yG,zG
        !   Calcul coque
        REAL,DIMENSION(3,3) :: Icoque
        REAL,DIMENSION(3) :: Gcoque
        !
        INTEGER i,j
        INTEGER,DIMENSION(Nface, 4) :: P
         ! Volumique mass of the fluid in KG/M**3
        REAL :: ct_RHO
        ! Gravity
        REAL :: ct_G


        DO j=1,nFace
            DO i=1,4
                Coin(i,1,j)=1.0*X(p(j, i))
                Coin(i,2,j)=1.0*Y(p(j, i))
                Coin(i,3,j)=1.0*Z(p(j, i))
            END DO
        END DO

        !   Mise a l echelle
        IF ((lambda).NE.(1.)) THEN
            WRITE(*,'(A,F7.3)') '   - Mesh is scaled by ', lambda
            xG=xG*lambda
            yG=yG*lambda
            zG=zG*lambda
            tX=tX*lambda
            tY=tY*lambda
            DO i=1,nFace
                DO j=1,4
                    Coin(j,1,i)=Coin(j,1,i)*lambda
                    Coin(j,2,i)=Coin(j,2,i)*lambda
                    Coin(j,3,i)=Coin(j,3,i)*lambda
                END DO
            END DO
        END IF
        !   Translate
        xG=xG+tX
        yG=yG+tY
        DO i=1,nFace
            DO j=1,4
                Coin(j,1,i)=Coin(j,1,i)+tX
                Coin(j,2,i)=Coin(j,2,i)+tY
            END DO
        END DO


        WRITE(*,*) ' -> Calculate hydrostatics '
        WRITE(*,*) ' '
        Np=0
        Nf=nFace
        DO i=1,nFace
            DO j=1,4
                X(Np+j)=Coin(j,1,i)-xG
                Y(Np+j)=Coin(j,2,i)-yG
                Z(Np+j)=Coin(j,3,i)+Tcol
                Facette(j,i)=Np+j
            END DO
            Np=Np+4
        END DO
        !    Calcul et sauvegarde de la matrice de raideur hydrostatique
        CALL HYDRO(X,Y,Z,FACETTE,NF,DEPLACEMENT,XF,YF,ZF,SF,KH,Xm,Ym,Zm,NPm,FACETTEm,NFm)

        !    Prise en compte de la symétrie
        IF (Nsym.EQ.1) THEN
            DEPLACEMENT=2.0*DEPLACEMENT
            YF=0.
            SF=2.0*SF
            KH(3,3)=2.*KH(3,3)
            KH(3,4)=0.
            KH(4,3)=0.
            KH(3,5)=2.*KH(3,5)
            KH(5,3)=KH(3,5)
            KH(4,4)=2.*KH(4,4)
            KH(4,5)=0.
            KH(5,4)=0.
            KH(5,5)=2.*KH(5,5)
        END IF
        WRITE(*,'(A,I3)') '   - Coordinates of buoyancy centre '
        WRITE(*,'(A,F7.3,A)') '     XB = ',XF+xG,'m'
        WRITE(*,'(A,F7.3,A)') '     YB = ',YF+yG,'m'
        WRITE(*,'(A,F7.3,A)') '     ZB = ',ZF,'m'
        WRITE(*,'(A,E14.7,A)') '    - Displacement  = ',DEPLACEMENT,' m^3'
        WRITE(*,'(A,E14.7,A)') '    - Waterplane area  = ',SF, ' m^2'
        WRITE(*,*) ' '
        IF ((ABS(XF).GT.1E-02).OR.(ABS(YF).GT.1E-02)) THEN
            WRITE(*,*) ' '
            WRITE(*,*) ' !!! WARNING !!! '
            WRITE(*,*) ' '
            WRITE(*,'(A,I3)') ' Buoyancy center and gravity center are not vertically aligned. '
            WRITE(*,*) ' This is not an equilibrium position.'
            WRITE(*,'(A,F7.3,1X,A,F7.3)') ' XF = ',XF+xG,' XG = ',xG
            WRITE(*,'(A,F7.3,1X,A,F7.3)') ' YF = ',YF+yG,' YG = ',yG
        END IF

        KH(4,4)=KH(4,4)+deplacement*ct_RHO*ct_G*(ZF-ZG)
        KH(5,5)=KH(5,5)+deplacement*ct_RHO*ct_G*(ZF-ZG)

        center_buoyancy(1) = XF
        center_buoyancy(2) = YF
        center_buoyancy(3) = ZF
        displacement = DEPLACEMENT
        waterplane_area = SF
        stifness = KH

    END SUBROUTINE COMPUTE_MESH_STATISTICS


END MODULE MESH_STATISTICS
