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
!   - A. Babarit
!
!--------------------------------------------------------------------------------------


!   Module containing utilities used to create Mesh
!
!   @author TCSASSEMBLER
!   @version 1.1
MODULE NETTING

USE HYDRO_NETTING

IMPLICIT NONE

CONTAINS

    SUBROUTINE calCol(NFm,Xm,Ym,Zm,Facettem,Tcol,nFacemx)
        !
        IMPLICIT NONE
        !   Maillage de la partie mouille
        INTEGER NFm,nFacemx
        REAL,DIMENSION(*) :: Xm,Ym,Zm
        INTEGER,DIMENSION(4,*) :: Facettem
        !   Hauteur de col
        REAL Tcol
        !   Locales
        INTEGER :: nFacem2,Np
        REAL,DIMENSION(4,3,nFacemx) :: Coinm2
        INTEGER :: i,j


        nFacem2=NFm
        DO i=1,NFm
            DO j=1,4
                Coinm2(j,1,i)=Xm(Facettem(j,i))
                Coinm2(j,2,i)=Ym(Facettem(j,i))
                Coinm2(j,3,i)=Zm(Facettem(j,i))-Tcol
            END DO
        END DO
        DO i=1,Nfm
            IF (((Coinm2(1,3,i)+Tcol).GT.-1.0E-03).AND.((Coinm2(2,3,i)+Tcol).GT.-1.0E-03)) THEN
                Nfacem2=Nfacem2+1
                Coinm2(1,1,Nfacem2)=Coinm2(2,1,i)
                Coinm2(1,2,Nfacem2)=Coinm2(2,2,i)
                Coinm2(1,3,Nfacem2)=Coinm2(2,3,i)
                Coinm2(2,1,Nfacem2)=Coinm2(1,1,i)
                Coinm2(2,2,Nfacem2)=Coinm2(1,2,i)
                Coinm2(2,3,Nfacem2)=Coinm2(1,3,i)
                Coinm2(3,1,Nfacem2)=Coinm2(1,1,i)
                Coinm2(3,2,Nfacem2)=Coinm2(1,2,i)
                Coinm2(3,3,Nfacem2)=0.
                Coinm2(4,1,Nfacem2)=Coinm2(2,1,i)
                Coinm2(4,2,Nfacem2)=Coinm2(2,2,i)
                Coinm2(4,3,Nfacem2)=0.
            ELSE
                IF (((Coinm2(2,3,i)+Tcol).GT.-1.0E-03).AND.((Coinm2(3,3,i)+Tcol).GT.-1.0E-03)) THEN
                    Nfacem2=Nfacem2+1
                    Coinm2(1,1,Nfacem2)=Coinm2(3,1,i)
                    Coinm2(1,2,Nfacem2)=Coinm2(3,2,i)
                    Coinm2(1,3,Nfacem2)=Coinm2(3,3,i)
                    Coinm2(2,1,Nfacem2)=Coinm2(2,1,i)
                    Coinm2(2,2,Nfacem2)=Coinm2(2,2,i)
                    Coinm2(2,3,Nfacem2)=Coinm2(2,3,i)
                    Coinm2(3,1,Nfacem2)=Coinm2(2,1,i)
                    Coinm2(3,2,Nfacem2)=Coinm2(2,2,i)
                    Coinm2(3,3,Nfacem2)=0.
                    Coinm2(4,1,Nfacem2)=Coinm2(3,1,i)
                    Coinm2(4,2,Nfacem2)=Coinm2(3,2,i)
                    Coinm2(4,3,Nfacem2)=0.
                ELSE
                    IF (((Coinm2(3,3,i)+Tcol).GT.-1.0E-03).AND.((Coinm2(4,3,i)+Tcol).GT.-1.0E-03)) THEN
                        Nfacem2=Nfacem2+1
                        Coinm2(1,1,Nfacem2)=Coinm2(4,1,i)
                        Coinm2(1,2,Nfacem2)=Coinm2(4,2,i)
                        Coinm2(1,3,Nfacem2)=Coinm2(4,3,i)
                        Coinm2(2,1,Nfacem2)=Coinm2(3,1,i)
                        Coinm2(2,2,Nfacem2)=Coinm2(3,2,i)
                        Coinm2(2,3,Nfacem2)=Coinm2(3,3,i)
                        Coinm2(3,1,Nfacem2)=Coinm2(3,1,i)
                        Coinm2(3,2,Nfacem2)=Coinm2(3,2,i)
                        Coinm2(3,3,Nfacem2)=0.
                        Coinm2(4,1,Nfacem2)=Coinm2(4,1,i)
                        Coinm2(4,2,Nfacem2)=Coinm2(4,2,i)
                        Coinm2(4,3,Nfacem2)=0.
                    ELSE
                        IF (((Coinm2(4,3,i)+Tcol).GT.-1.0E-03).AND.((Coinm2(1,3,i)+Tcol).GT.-1.0E-03)) THEN
                            Nfacem2=Nfacem2+1
                            Coinm2(1,1,Nfacem2)=Coinm2(1,1,i)
                            Coinm2(1,2,Nfacem2)=Coinm2(1,2,i)
                            Coinm2(1,3,Nfacem2)=Coinm2(1,3,i)
                            Coinm2(2,1,Nfacem2)=Coinm2(4,1,i)
                            Coinm2(2,2,Nfacem2)=Coinm2(4,2,i)
                            Coinm2(2,3,Nfacem2)=Coinm2(4,3,i)
                            Coinm2(3,1,Nfacem2)=Coinm2(4,1,i)
                            Coinm2(3,2,Nfacem2)=Coinm2(4,2,i)
                            Coinm2(3,3,Nfacem2)=0.
                            Coinm2(4,1,Nfacem2)=Coinm2(1,1,i)
                            Coinm2(4,2,Nfacem2)=Coinm2(1,2,i)
                            Coinm2(4,3,Nfacem2)=0.
                        END IF
                    END IF
                END IF
            END IF
        END DO
        Np=0
        NFm=Nfacem2
        DO i=1,nFacem2
            DO j=1,4
                Xm(Np+j)=Coinm2(j,1,i)
                Ym(Np+j)=Coinm2(j,2,i)
                Zm(Np+j)=Coinm2(j,3,i)
                Facettem(j,i)=Np+j
            END DO
            Np=Np+4
        END DO

    END SUBROUTINE


    SUBROUTINE coque(X,Y,Z,facettes,NF,Deplacement,Icoque,Gcoque)

        IMPLICIT NONE

        !    Globale
        REAL,PARAMETER :: PI=3.141592653589
        !    Description du maillage
        INTEGER Nf
        REAL,DIMENSION(*) :: X,Y,Z
        INTEGER,DIMENSION(4,*) :: Facettes
        REAL Deplacement
        !    Caracteristiques de la coque
        REAL,DIMENSION(3) :: Gcoque
        REAL mcoque
        REAL,DIMENSION(3,3) :: Icoque
        REAL decoque
        !    Locales
        REAL,DIMENSION(3) :: U,V,W
        REAL,DIMENSION(3,Nf) :: CdG
        REAL :: N1,N2
        REAL,DIMENSION(Nf) :: Aire
        !    Indices
        INTEGER i,j

        mcoque=0.
        Gcoque(1)=0.
        Gcoque(2)=0.
        Gcoque(3)=0.
        DO i=1,nF
            CdG(1,i)=0.25*(X(Facettes(1,i))+X(Facettes(2,i))+X(Facettes(3,i))+X(Facettes(4,i)))
            CdG(2,i)=0.25*(Y(Facettes(1,i))+Y(Facettes(2,i))+Y(Facettes(3,i))+Y(Facettes(4,i)))
            CdG(3,i)=0.25*(Z(Facettes(1,i))+Z(Facettes(2,i))+Z(Facettes(3,i))+Z(Facettes(4,i)))
            U(1)=X(Facettes(2,i))-X(Facettes(1,i))
            U(2)=Y(Facettes(2,i))-Y(Facettes(1,i))
            U(3)=Z(Facettes(2,i))-Z(Facettes(1,i))
            V(1)=X(Facettes(3,i))-X(Facettes(1,i))
            V(2)=Y(Facettes(3,i))-Y(Facettes(1,i))
            V(3)=Z(Facettes(3,i))-Z(Facettes(1,i))
            CALL prdvct(U,V,W)
            N1=0.5*SQRT(W(1)*W(1)+W(2)*W(2)+W(3)*W(3))
            U(1)=X(Facettes(3,i))-X(Facettes(1,i))
            U(2)=Y(Facettes(3,i))-Y(Facettes(1,i))
            U(3)=Z(Facettes(3,i))-Z(Facettes(1,i))
            V(1)=X(Facettes(4,i))-X(Facettes(1,i))
            V(2)=Y(Facettes(4,i))-Y(Facettes(1,i))
            V(3)=Z(Facettes(4,i))-Z(Facettes(1,i))
            CALL prdvct(U,V,W)
            N2=0.5*SQRT(W(1)*W(1)+W(2)*W(2)+W(3)*W(3))
            Aire(i)=N1+N2
            mcoque=mcoque+Aire(i)
            Gcoque(1)=Gcoque(1)+Aire(i)*CdG(1,i)
            Gcoque(2)=Gcoque(2)+Aire(i)*CdG(2,i)
            Gcoque(3)=Gcoque(3)+Aire(i)*CdG(3,i)
        END DO
        decoque=Deplacement*1000.0/mcoque
        Gcoque(1)=Gcoque(1)/mcoque
        Gcoque(2)=Gcoque(2)/mcoque
        Gcoque(3)=Gcoque(3)/mcoque
        DO i=1,3
            DO j=1,3
                Icoque(i,j)=0.
            END DO
        END DO
        DO i=1,nF
            Icoque(1,1)=Icoque(1,1)+Aire(i)*decoque*((CdG(2,i)-Gcoque(2))**2+(CdG(3,i)-Gcoque(3))**2)
            Icoque(1,2)=Icoque(1,2)+Aire(i)*decoque*((CdG(1,i)-Gcoque(1))*(CdG(2,i)-Gcoque(2)))
            Icoque(1,3)=Icoque(1,3)+Aire(i)*decoque*((CdG(1,i)-Gcoque(1))*(CdG(3,i)-Gcoque(3)))
            Icoque(2,2)=Icoque(2,2)+Aire(i)*decoque*((CdG(1,i)-Gcoque(1))**2+(CdG(3,i)-Gcoque(3))**2)
            Icoque(2,3)=Icoque(2,3)+Aire(i)*decoque*((CdG(2,i)-Gcoque(2))*(CdG(3,i)-Gcoque(3)))
            Icoque(3,3)=Icoque(3,3)+Aire(i)*decoque*((CdG(1,i)-Gcoque(1))**2+(CdG(2,i)-Gcoque(2))**2)
        END DO
        !
        Icoque(2,1)=Icoque(1,2)
        Icoque(3,1)=Icoque(1,3)
        Icoque(3,2)=Icoque(2,3)

        RETURN

    END SUBROUTINE


    SUBROUTINE ExMaillage(ID,DSCRPT,X,Y,Z,NP,facette,NF,NFMX)

        USE MIdentification

        implicit none

        ! Entrees
        TYPE(TID) :: ID,DSCRPT
        INTEGER NP,NF,NFMX
        REAL,DIMENSION(*) :: X,Y,Z
        INTEGER,DIMENSION(4,*) :: facette
        ! Locales
        REAL,DIMENSION(NFMX) :: Xn,Yn,Zn    ! Normale a la facette
        REAL,DIMENSION(NFMX) :: Xg,Yg,Zg    ! CdG des facettes
        REAL,DIMENSION(NFMX) :: Aire        ! Aire des facettes
        REAL,DIMENSION(NFMX) :: Dist        ! Plus grande longueur des facettes
        INTEGER i
        REAL T                              ! Tirant d eau (valeur adimensionnalisation ACHIL3D)
        REAL,DIMENSION(3) :: u,v,w
        REAL norme,norme1,norme2,norme3

        ! Recherche du tirant d eau
        T=0.0
        do i=1,NP
            IF (ABS(Z(i)).GT.T) T=ABS(Z(i))
        end do
        ! Calcul des normales
        do i=1,NF
            u(1)=X(facette(2,i))-X(facette(1,i))
            u(2)=Y(facette(2,i))-Y(facette(1,i))
            u(3)=Z(facette(2,i))-Z(facette(1,i))
            norme1=u(1)**2+u(2)**2+u(3)**2
            v(1)=X(facette(3,i))-X(facette(1,i))
            v(2)=Y(facette(3,i))-Y(facette(1,i))
            v(3)=Z(facette(3,i))-Z(facette(1,i))
            norme2=v(1)**2+v(2)**2+v(3)**2
            call prodvect(u,v,w)
            Xn(i)=w(1)
            Yn(i)=w(2)
            Zn(i)=w(3)
            u(1)=X(facette(3,i))-X(facette(1,i))
            u(2)=Y(facette(3,i))-Y(facette(1,i))
            u(3)=Z(facette(3,i))-Z(facette(1,i))
            v(1)=X(facette(4,i))-X(facette(1,i))
            v(2)=Y(facette(4,i))-Y(facette(1,i))
            v(3)=Z(facette(4,i))-Z(facette(1,i))
            norme3=v(1)**2+v(2)**2+v(3)**2
            call prodvect(u,v,w)
            Xn(i)=0.5*(w(1)+Xn(i))
            Yn(i)=0.5*(w(2)+Yn(i))
            Zn(i)=0.5*(w(3)+Zn(i))
            Dist(i)=max(norme1,norme2,norme3)
            Aire(i)=sqrt(Xn(i)**2+Yn(i)**2+Zn(i)**2)
            norme=Aire(i)
            Xn(i)=Xn(i)/norme
            Yn(i)=Yn(i)/norme
            Zn(i)=Zn(i)/norme
            Xg(i)=0.25*(X(facette(1,i))+X(facette(2,i))+X(facette(3,i))+X(facette(4,i)))
            Yg(i)=0.25*(Y(facette(1,i))+Y(facette(2,i))+Y(facette(3,i))+Y(facette(4,i)))
            Zg(i)=0.25*(Z(facette(1,i))+Z(facette(2,i))+Z(facette(3,i))+Z(facette(4,i)))
        end do
        ! Maillage AQUAPLUS
        ! Fichier Tecplot
        open(10,file=ID%ID(1:ID%lID)//'/Mesh/'//DSCRPT%ID(1:DSCRPT%lID)//'.tec')
        write(10,*) 'ZONE N=',np,', E=',nf,' , F=FEPOINT,ET=QUADRILATERAL'
        do i=1,np
            write(10,'(6(2X,E14.7))') X(i),Y(i),Z(i),0.0,0.0,0.0
        end do
        do i=1,nf
            write(10,'(I6,3(2X,I6))') facette(1,i),facette(2,i),facette(3,i),facette(4,i)
        end do
        write(10,*) 'ZONE t="normales", F=POINT, I=',nf
        do i=1,nf
            write(10,'(6(2X,E14.7))') Xg(i),Yg(i),Zg(i),Xn(i),Yn(i),Zn(i)
        end do
        close(10)
        ! Fichier de maillage
        open(10,file=ID%ID(1:ID%lID)//'/Mesh/'//DSCRPT%ID(1:DSCRPT%lID)//'.dat')
        write(10,'(20X,I1,10X,I1)') 2,1
        do i=1,np
            write(10,'(10X,I4,3(10X,F14.7))') i,X(i),Y(i),Z(i)
        end do
        write(10,'(10X,I4,3(10X,F5.2))') 0,0.,0.,0.
        do i=1,nf
            write(10,'(4(10X,I4))') facette(1,i),facette(2,i),facette(3,i),facette(4,i)
        end do
        write(10,'(4(10X,I1))') 0,0,0,0
        close(10)
        ! Fichiers de maillage
        open(10,file=ID%ID(1:ID%lID)//'/Mesh/'//DSCRPT%ID(1:DSCRPT%lID)//'_info.dat')
        write(10,'(I7,1X,I7,A)') np,nf, ' Number of points and number of panels'
        close(10)

        RETURN

    end SUBROUTINE

    !*******************************************************************
    !
    !    Calcul du produit vectoriel w = u ^ v
    !
    !******************************************************************

    SUBROUTINE prodvect(u,v,w)

        IMPLICIT NONE

        REAL,DIMENSION(3) :: u,v,w

        w(1)=u(2)*v(3)-u(3)*v(2)
        w(2)=u(3)*v(1)-u(1)*v(3)
        w(3)=u(1)*v(2)-u(2)*v(1)

        RETURN

    END SUBROUTINE




    SUBROUTINE MAILLAGE(nFace,Coin,maille,X,Y,Z,NP,facette,NF)

        implicit none

        ! Entrees
        REAL maille                    !   Longueur maximale de maille
        INTEGER nFace
        REAL,DIMENSION(4,3,*) :: Coin
        ! Sorties
        INTEGER NP,NF
        REAL,DIMENSION(*) :: X,Y,Z
        INTEGER,DIMENSION(4,*) :: facette
        ! Locales
        INTEGER i,j,k
        INTEGER,PARAMETER :: ndpmx=10000,ndfmx=10000
        INTEGER :: ndp,ndf            !    ndp : nombre de points
                                    !    nf : nombre de facettes
        REAL,DIMENSION(3,ndpmx) :: noeud    !    noeuds : du maillage
        INTEGER,DIMENSION(4,ndfmx) :: face    !    faces : du maillage
        INTEGER,DIMENSION(ndfmx) :: indx    !    Index de maillage
        REAL,DIMENSION(3) :: X1,X2,X3,X4
        REAL norme1,norme2

        ! Analyse des faces. Si une face est trop petite par rapport a la taille de la maille
        ! elle ne sera pas mailleeï¿½et constituera une facette aï¿½part entiere.
        do i=1,Nface
            do j=1,3
                X1(j)=Coin(1,j,i)
                X2(j)=Coin(2,j,i)
                X3(j)=Coin(3,j,i)
                X4(j)=Coin(4,j,i)
            end do
            call nprodvect(X1,X2,X3,norme1)
            call nprodvect(X1,X3,X4,norme2)
            !        if (0.5*(norme1+norme2).LT.maille) then
            !            indx(i)=0
            !        else
            indx(i)=1
        !        end if
        end do
        ! Calcul du maillage
        np=0
        nf=0
        ! Maillage des faces
        do i=1,Nface
            if (indx(i).EQ.1) then
                do j=1,3
                    X1(j)=Coin(1,j,i)
                    X2(j)=Coin(2,j,i)
                    X3(j)=Coin(3,j,i)
                    X4(j)=Coin(4,j,i)
                end do
                call plan(X2,X3,X4,X1,maille,noeud,ndp,face,ndf)
                do j=1,ndp
                    X(j+np)=noeud(1,j)
                    Y(j+np)=noeud(2,j)
                    Z(j+np)=noeud(3,j)
                end do
                do j=1,ndf
                    do k=1,4
                        facette(k,j+nf)=face(k,j)+np
                    end do
                end do
                np=np+ndp
                nf=nf+ndf
            else
                do j=1,4
                    X(j+np)=Coin(j,1,i)
                    Y(j+np)=Coin(j,2,i)
                    Z(j+np)=Coin(j,3,i)
                end do
                do k=1,4
                    facette(k,1+nf)=np+k
                end do
                np=np+4
                nf=nf+1
            end if
        end do

        RETURN

    end SUBROUTINE
    !*************************************************************
    !
    !    Cette subroutine realise le maillage d'une face
    !    sommmets X1,X2,X3,X4 avec des mailles de dimension max maille.
    !        X3 peut etre confondu avec X4
    !        Les noeuds sont stockes dans noeud.
    !        Les facettes dans facette.
    !
    !**************************************************************

    SUBROUTINE plan(X1,X2,X3,X4,maille,noeud,n,facette,nf)

        INTEGER i,j,k,n,nf,nj,ni
        REAL X1,X2,X3,X4,Pi,Pf
        DIMENSION X1(3),X2(3),X3(3),X4(3),Pi(3),Pf(3)
        REAL noeud,maille
        INTEGER,DIMENSION(4,*) :: facette
        DIMENSION noeud(3,*)
        REAL norme1,norme2,norme3,norme
        REAL angle1,angle2,angle3
        REAL A1(3),A2(3),A3(3)
        !    REAL,PARAMETER :: njmn=3        ! Nombre minimum de points de discretisation
        INTEGER :: njmn

        OPEN(10,FILE='mesh.cal')
        DO i=1,5
            READ(10,*)
        END DO
        READ(10,*) njmn
        CLOSE(10)
        call nprodvect(X1,X2,X3,norme1)
        call nprodvect(X1,X2,X4,norme2)
        call nprodvect(X1,X3,X4,norme3)
        if ((norme1.EQ.0.0).AND.(norme2.EQ.0.0).AND.(norme3.EQ.0.0)) then
            nf=0
            n=0
        else
            if ((X3(1).EQ.X4(1)).AND.(X3(2).EQ.X4(2)).AND.(X3(3).EQ.X4(3))) then
                norme1=((X2(1)-X1(1))**2+(X2(2)-X1(2))**2+(X2(3)-X1(3))**2)**0.5
                norme2=((X3(1)-X2(1))**2+(X3(2)-X2(2))**2+(X3(3)-X2(3))**2)**0.5
                norme3=((X1(1)-X3(1))**2+(X1(2)-X3(2))**2+(X1(3)-X3(3))**2)**0.5
                Angle1=((X2(1)-X1(1))*(X3(3)-X1(3))+(X2(3)-X1(3))*(X3(3)-X1(3))+(X2(3)-X1(3))*(X3(3)-X1(3)))/norme1/norme3
                Angle2=((X3(1)-X2(1))*(X1(3)-X2(3))+(X3(3)-X2(3))*(X1(3)-X2(3))+(X3(3)-X2(3))*(X1(3)-X2(3)))/norme2/norme1
                Angle3=((X1(1)-X3(1))*(X2(3)-X3(3))+(X1(3)-X3(3))*(X2(3)-X3(3))+(X1(3)-X3(3))*(X2(3)-X3(3)))/norme3/norme2
                if ((angle1.GT.angle2).AND.(angle1.GT.angle3)) then
                    do i=1,3
                        A1(i)=X1(i)
                        A2(i)=X2(i)
                        A3(i)=X3(i)
                    end do
                else
                    if (angle2.GT.angle3) then
                        do i=1,3
                            A1(i)=X2(i)
                            A2(i)=X3(i)
                            A3(i)=X1(i)
                        end do
                    else
                        do i=1,3
                            A1(i)=X3(i)
                            A2(i)=X1(i)
                            A3(i)=X2(i)
                        end do
                    end if
                end if
                norme1=((A2(1)-A1(1))**2+(A2(2)-A1(2))**2+(A2(3)-A1(3))**2)**0.5
                norme2=((A3(1)-A2(1))**2+(A3(2)-A2(2))**2+(A3(3)-A2(3))**2)**0.5
                norme3=((A1(1)-A3(1))**2+(A1(2)-A3(2))**2+(A1(3)-A3(3))**2)**0.5
                if (norme1.LT.norme3) then
                    nj=int(norme1/maille)+njmn
                    n=0
                    do j=1,nj
                        do k=1,3
                            Pi(k)=A1(k)+(A2(k)-A1(k))*(j-1)/(nj-1)
                        end do
                        norme=((Pi(1)-A1(1))**2+(Pi(2)-A1(2))**2+(Pi(3)-A1(3))**2)**0.5
                        do k=1,3
                            Pf(k)=A1(k)+(A3(k)-A1(k))*norme/norme1
                        end do
                        do i=1,j
                            n=n+1
                            do k=1,3
                                if (j.GT.1) then
                                    noeud(k,n)=Pi(k)+(Pf(k)-Pi(k))*(i-1)/(j-1)
                                else
                                    noeud(k,n)=Pi(k)
                                end if
                            end do
                        end do
                    end do
                    nf=0
                    do j=1,nj-1
                        do i=1,j-1
                            nf=nf+1
                            facette(1,nf)=nf
                            facette(2,nf)=nf+j
                            facette(3,nf)=nf+j+1
                            facette(4,nf)=nf+1
                        end do
                        nf=nf+1
                        facette(1,nf)=nf
                        facette(2,nf)=nf+j
                        facette(3,nf)=nf+j+1
                        facette(4,nf)=nf
                    end do
                else
                    nj=int(norme3/maille)+njmn
                    n=0
                    do j=1,nj
                        do k=1,3
                            Pi(k)=A1(k)+(A3(k)-A1(k))*(j-1)/(nj-1)
                        end do
                        norme=((Pi(1)-A1(1))**2+(Pi(2)-A1(2))**2+(Pi(3)-A1(3))**2)**0.5
                        do k=1,3
                            Pf(k)=A1(k)+(A2(k)-A1(k))*norme/norme3
                        end do
                        do i=1,j
                            n=n+1
                            do k=1,3
                                if (j.GT.1) then
                                    noeud(k,n)=Pi(k)+(Pf(k)-Pi(k))*(i-1)/(j-1)
                                else
                                    noeud(k,n)=Pi(k)
                                end if
                            end do
                        end do
                    end do
                    nf=0
                    do j=1,nj-1
                        do i=1,j-1
                            nf=nf+1
                            facette(1,nf)=nf
                            facette(4,nf)=nf+j
                            facette(3,nf)=nf+j+1
                            facette(2,nf)=nf+1
                        end do
                        nf=nf+1
                        facette(1,nf)=nf
                        facette(3,nf)=nf+j
                        facette(2,nf)=nf+j+1
                        facette(4,nf)=nf
                    end do
                end if
            else
                ni=int(((X2(1)-X1(1))**2+(X2(2)-X1(2))**2+(X2(3)-X1(3))**2)**0.5/maille)+njmn
                nj=int(((X4(1)-X1(1))**2+(X4(2)-X1(2))**2+(X4(3)-X1(3))**2)**0.5/maille)+njmn
                n=0
                do j=1,nj
                    do k=1,3
                        Pi(k)=X1(k)+(X4(k)-X1(k))*(j-1)/(nj-1)
                        Pf(k)=X2(k)+(X3(k)-X2(k))*(j-1)/(nj-1)
                    end do
                    do i=1,ni
                        n=n+1
                        do k=1,3
                            noeud(k,n)=Pi(k)+(Pf(k)-Pi(k))*(i-1)/(ni-1)
                        end do
                    end do
                end do
                nf=0
                do j=1,nj-1
                    do i=1,ni-1
                        nf=nf+1
                        facette(1,nf)=i+(j-1)*ni
                        facette(2,nf)=i+1+(j-1)*ni
                        facette(4,nf)=i+j*ni
                        facette(3,nf)=i+1+j*ni
                    end do
                end do
            end if
        end if

    END SUBROUTINE

    !*******************************************************************
    !
    !    Calcule la norme du  produit vectoriel des vecteurs X1X2 et X1X3
    !
    !******************************************************************

    SUBROUTINE nprodvect(X1,X2,X3,norme)

        IMPLICIT NONE

        REAL X1(3),X2(3),X3(3)
        REAL N(3)
        REAL norme

        N(1)=(X2(2)-X1(2))*(X3(3)-X1(3))-(X2(3)-X1(3))*(X3(2)-X1(2))
        N(2)=(X2(3)-X1(3))*(X3(1)-X1(1))-(X2(1)-X1(1))*(X3(3)-X1(3))
        N(3)=(X2(1)-X1(1))*(X3(2)-X1(2))-(X2(2)-X1(2))*(X3(1)-X1(1))
        norme=sqrt(N(1)**2+N(2)**2+N(3)**2)

        RETURN

    END SUBROUTINE



END MODULE
