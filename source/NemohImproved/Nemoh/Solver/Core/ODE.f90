!--------------------------------------------------------------------------------------
!
!Copyright (C) 2014 TopCoder Inc., All Rights Reserved.
!
!--------------------------------------------------------------------------------------

!  Contest Code Acceleration of the Calculation of Influence Coefficients of Nemoh
!
!   @author yedtoss
!   @version 1.0

!********************************************************************
!*             Differential equations of order N                    *
!*             by Runge-Kutta method of order 4                     *
!* ---------------------------------------------------------------- *
!* Reference: "Analyse en Turbo Pascal versions 5.5 et 6.0  By Marc *
!*             DUCAMP et Alain REVERCHON - Eyrolles, Paris 1991"    *
!*             [BIBLI 03].                                          *
!*                                                                  *
!*                             F90 Version By J-P Moreau, Paris     *
!*                                    (www.jpmoreau.fr)             *
!* ---------------------------------------------------------------- *
!* SAMPLE RUN:                                                      *
!*                                                                  *
!* Example: integrate y"=(2y'y'+yy)/y from x=4 to x=6               *
!*          with initial conditions: y(4)=2 and y'(4)=-2tan(1)      *
!*                                                                  *
!* Exact solution is:   y = 2cos(1)/cos(x-5)                        *
!*                                                                  *
!*        DIFFERENTIAL EQUATION WITH 1 VARIABLE OF ORDER N          *
!*              of type y(n) = f(y(n-1),...,y',y,x)                 *
!*                                                                  *
!*   order of equation: 2                                           *
!*   begin value x    : 4                                           *
!*   end value x      : 6                                           *
!*   y0 value at x0   : 2                                           *
!*   y1 value at x0   : -3.114815                                   *
!*   number of points : 11                                          *
!*   finesse          : 20                                          *
!*                                                                  *
!*       X            Y                                             *
!* ---------------------------                                      *
!*   4.000000     2.000000                                          *
!*   4.200000     1.551018                                          *
!*   4.400000     1.309291                                          *
!*   4.600000     1.173217                                          *
!*   4.800000     1.102583                                          *
!*   5.000000     1.080605                                          *
!*   5.200000     1.102583                                          *
!*   5.400000     1.173217                                          *
!*   5.600000     1.309291                                          *
!*   5.800000     1.551018                                          *
!*   6.000000     2.000000                                          *
!*                                                                  *
!********************************************************************

module ODE

    USE COMMON_TYPE
    implicit none
    contains



function Equadifnc(fp, xi,xf,yi,m,n,fi, param)
!***************************************************************************
!*        SOLVING DIFFERENTIAL EQUATIONS WITH 1 VARIABLE OF ORDER N        *
!*                of type y(n) = f(y(n-1),y(n-2),...,y',y,x)               *
!* ----------------------------------------------------------------------- *
!*  INPUTS:                                                                *
!*    xi, xf    Begin, end values of variable x                            *
!*    yi        Begin values at xi (of f and derivatives)                  *
!*    m         Number of points to calculate                              *
!*    n         Order of differential equation                             *
!*    fi        finesse (number of intermediate points)                    *
!*                                                                         *
!*  OUTPUTS:                                                               *
!*    t         real vector storing m results for function y               *
!* ----------------------------------------------------------------------- *
!*  EXAMPLE:    y' = (2 y'y' + yy) / y with y(4) = 2, y'(4) = -2*tan(1)    *
!*              Exact solution:  y = (2 cos(1))/cos(x-5)                   *
!***************************************************************************
integer n,m,fi
real  xi,xf
complex yi(0:10),t(50)
complex Equadifnc(50)
real  h,x,a,b,c,d
integer i,j,k,ni,n1,n2
complex  ta(0:10),tb(0:10),tc(0:10),td(0:10),y(0:10),z(0:10)
TYPE(ParamsCommonInf) :: param

interface
complex function fp(x,y, param1)
USE COMMON_TYPE
real x
complex y(0:10)
TYPE(ParamsCommonInf) :: param1

end function fp
end interface
   if (fi<1) then
    return
   end if

   if ( m == 1 ) then
    h = (xf - xi) / fi
   else
    h = (xf - xi) / fi / (m-1)
   endif
   n1=n-1; n2=n-2
   t(1)=yi(0)
   do k=0, n1
     y(k)=yi(k); z(k)=yi(k)
   end do
   do i=1, m
     ni=(i-1)*fi-1
     do j=1, fi
       x=xi+h*(ni+j)
       do k=0, n1
         y(k)=z(k)
       end do
       a=h*fp(x,y, param)
       do k=0, n2
         ta(k)=h*y(k+1); y(k)=z(k)+ta(k)/2.d0
       end do
       y(n1)=z(n1)+a/2.d0
       x=x+h/2
       b=h*fp(x,y, param)
       do k=0, n2
         tb(k)=h*y(k+1); y(k)=z(k)+tb(k)/2.d0
       end do
       y(n1)=z(n1)+b/2.d0
       c=h*fp(x,y, param)
       do k=0, n2
         tc(k)=h*y(k+1); y(k)=z(k)+tc(k)
       end do
       y(n1)=z(n1)+c
       x=x+h/2
       d=h*fp(x,y, param)
       do k=0, n2
         z(k)=z(k)+(ta(k)+2.d0*tb(k)+2.d0*tc(k)+h*y(k+1))/6.d0
       end do
       z(n1)=z(n1)+(a+b+b+c+c+d)/6.d0
     end do
     t(i+1)=z(0)
   end do
   !call subroutine Affiche
   !call Affiche(m,xi,xf,t)




   Equadifnc = t
   return
end function equadifnc
end module ODE


!Example: y'=(2y'y'+yy)/y
!real*8 FUNCTION fp(x,y)
!real*8 x,y(10)
  !if (dabs(y(0))<1.d-12)  y(0)=1.d-12
  !fp=(2.d0*y(1)*y(1)+y(0)*y(0))/y(0)
!end




!print solutions
!Subroutine Affiche(m,xi,xf,t)
 ! integer i,m
  !real*8 h,x,xi,xf,t(50)
  !h=(xf-xi)/dfloat(m-1)
  !x=xi-h
  !print *,' '
  !print *, '      X         Y     '
  !print *, '----------------------'
  !do i=1, m
   ! x=x+h
    !write(*,100)  x, t(i)
  !end do
  !print *,' '
  !return
!100 format(' ',f9.6,'  ',f9.6)
!end

!End of file teqdifn.f90
