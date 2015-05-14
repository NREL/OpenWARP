!--------------------------------------------------------------------------------------
!
!Copyright (C) 2014 TopCoder Inc., All Rights Reserved.
!
!--------------------------------------------------------------------------------------

!   This module computes utility for numerical integration and derivation of functions
!   Currently the romberg integration algorithm is implemented
!   and the Derivative by center differences/Richardson extrapolation is implemented
!   They work for real function or complex function of real parameter
!
! Contest Drift forces and QTF Implementation of Nemoh
!
! Changes in version 1.1 (Implementation of Higher Order Panel Methods)
!       Added COMMON_TYPE module as dependency

! Changes in version 1.2 (Dipoles Implementation in NEMOH)
!       Switch from ParamsCommon to ParamsCommonInf and added COMPUTE_NEXP subroutine
!
!   @author yedtoss
!   @version 1.2


module UTILITY
    USE COMMON_TYPE
    implicit none



contains

    SUBROUTINE trapzd(func,a,b,s,n,param)
        INTEGER n
        REAL a,b,s,func
        EXTERNAL func
        INTEGER it,j
        REAL del,sum,tnm,x
        TYPE(ParamsCommonInf) :: param
        if (n.eq.1) then
            s=0.5*(b-a)*(func(a, param)+func(b, param))
        else
            it=2**(n-2)
            tnm=it
            del=(b-a)/tnm
            x=a+0.5*del
            sum=0.
            do j=1,it
                sum=sum+func(x, param)
                x=x+del
            end do
            s=0.5*(s+(b-a)*sum/tnm)
        endif
        return
    END SUBROUTINE trapzd

    FUNCTION qtrap(func,a,b, param)
        INTEGER JMAX
        REAL a,b,func,qtrap,EPS
        EXTERNAL func
        PARAMETER (EPS=1.e-6, JMAX=20)
        TYPE(ParamsCommonInf) :: param
        INTEGER j
        REAL olds
        olds=-1.e30
        do j=1,JMAX
            call trapzd(func,a,b,qtrap,j, param)
            if (abs(qtrap-olds).lt.EPS*abs(olds)) return
            if (qtrap.eq.0..and.olds.eq.0..and.j.gt.6) return
            olds=qtrap
        end do
    !pause 'too many steps in qtrap'
    END FUNCTION qtrap



    SUBROUTINE trapzdc(func,a,b,s,n, param)
        INTEGER n
        REAL a,b
        COMPLEX func, sum, s
        EXTERNAL func
        INTEGER it,j
        REAL del,tnm,x
        TYPE(ParamsCommonInf) :: param
        COMPLEX,PARAMETER :: II=CMPLX(0.,1.)
        COMPLEX,PARAMETER :: ZERO=CMPLX(0.,0.)
        if (n.eq.1) then
            s=0.5*(b-a)*(func(a, param)+func(b, param))
        else
            it=2**(n-2)
            tnm=it
            del=(b-a)/tnm
            x=a+0.5*del
            sum=zero
            do j=1,it
                sum=sum+func(x, param)
                x=x+del
            end do
            s=0.5*(s+(b-a)*sum/tnm)
        endif
        return
    END SUBROUTINE trapzdc

    FUNCTION qtrapc(func,a,b, param)
        INTEGER JMAX
        REAL a,b,EPS
        COMPLEX qtrapc, func
        EXTERNAL func
        PARAMETER (EPS=1.e-6, JMAX=20)
        INTEGER j
        COMPLEX olds
        TYPE(ParamsCommonInf) :: param
        COMPLEX,PARAMETER :: II=CMPLX(0.,1.)
        COMPLEX,PARAMETER :: ZERO=CMPLX(0.,0.)
        olds=zero
        do j=1,JMAX
            call trapzdc(func,a,b,qtrapc,j, param)
            if (abs(qtrapc-olds).lt.EPS*abs(olds)) return
            if (qtrapc.eq.0..and.olds.eq.0..and.j.gt.6) return
            olds=qtrapc
        end do
    !pause 'too many steps in qtrap'
    END FUNCTION qtrapc






    real function deriv_richardson(f,x,n,h_ini, param)
        integer, intent(in) :: n
        real, dimension (0:n,0:n) :: d
        real, intent(in) :: x
        real :: q, h_ini, h
        integer :: i, j
        real, external :: f
        TYPE(ParamsCommonInf) :: param

        h = h_ini

        do  i=0,n
            d(i,0)=(f(x+h, param)-f(x-h, param))/(2.0*h)
            q=4.0
            do  j=0,i-1
                d(i,j+1)=d(i,j)+(d(i,j)-d(i-1,j))/(q-1.0)
                deriv_richardson = d(i,j+1)
                q=4.0*q
            end do
            h=h/2.0
        end do
    end function deriv_richardson


    complex function deriv_richardson_complex(f,x,n,h_ini, param)
        integer, intent(in) :: n
        complex, dimension (0:n,0:n) :: d
        real, intent(in) :: x
        real :: q, h, h_ini
        integer :: i, j
        complex, external :: f
        TYPE(ParamsCommonInf) :: param

        h = h_ini

        do  i=0,n
            d(i,0)=(f(x+h, param)-f(x-h, param))/(2.0*h)
            q=4.0
            do  j=0,i-1
                d(i,j+1)=d(i,j)+(d(i,j)-d(i-1,j))/(q-1.0)
                deriv_richardson_complex = d(i,j+1)
                q=4.0*q
            end do
            h=h/2.0
        end do
    end function deriv_richardson_complex



    function romberg_trap ( f, a, b, tol, maxM, param )

        !*****************************************************************************80
        !
        !! ROMBERG_TRAP approximates an integral by extrapolating the trapezoidal rule.
        !
        !  Discussion:
        !
        !    This routine computes a sequence of integral estimates involving the
        !    trapezoid rule.  Extrapolation of successive trapezoidal estimates
        !    produces an estimate with a higher rate of convergence.
        !
        !  Licensing:
        !
        !    This code is distributed under the GNU LGPL license.
        !
        !  Modified:
        !
        !    28 February 2013
        !
        !  Author:
        !
        !    John Burkardt
        !
        !  Parameters:
        !
        !    Input, real  A, B, the endpoints of the interval.
        !
        !    Input, external real  F, the name of the function
        !    which evaluates the integrand, and whose form is:
        !      function f ( x )
        !      real  f
        !      real  x
        !      f = ...
        !
        !    Input, real  TOL, the tolerance.  When two successive integral
        !    estimates differ by less than this value, the iteration will halt.
        !
        !    Output, integer  M, the number of trapezoidal estimates
        !    that were required.
        !
        !    Output, real  ROMBERG_TRAP, the integral estimate.
        !
        implicit none

        real  a
        real b
        !real , external :: f
        integer  m, maxM
        real  q
        real  q_new
        real  r
        real  r_new
        real  romberg_trap
        real  tol
        !real  trapezoid_refine
        TYPE(ParamsCommonInf) :: param

        interface

            real function f(theta, param)

                USE COM_VAR
                TYPE(ParamsCommonInf) :: param
                REAL, INTENT(IN) :: theta

            end function f


        end interface

        q_new = 0.0D+00
        r_new = 0.0D+00
        m = 1

        do

            q = q_new
            q_new = trapezoid_refine ( a, b, m, f, q, param )

            if ( m == 1 ) then
                r_new = q_new
            else
                r = r_new
                r_new = ( 4.0D+00 * q_new - q ) / 3.0D+00
                if ( abs ( r_new - r ) .lt. tol * ( 1.0 + abs ( r_new ) ) ) then
                    exit
                end if
            end if

            if ( maxM <= m ) then
                !write ( *, '(a)' ) ' '
                !write ( *, '(a)' ) 'ROMBERG_TRAP!'
                !write ( *, '(a)' ) '  No convergence in 20 iterations.'
                !write ( *, '(a)' ) '  The algorithm is halting.'
                return
            end if

            m = m + 1

        end do

        romberg_trap = r_new

        return
    end function romberg_trap

    function trapezoid_refine ( a, b, m, f, q, param)

        !*****************************************************************************80
        !
        !! TRAPEZOID_REFINE carries out a step of trapezoidal refinement.
        !
        !  Discussion:
        !
        !    This routine is designed to be an efficient way to carry out a
        !    sequence of integral estimates, using the trapezoidal rule
        !    and a nested sequence of evaluation points, in such a way that
        !    the function is not re-evaluated unnecessarily.
        !
        !    The user calls first with M = 1 and Q = 0 to get a 2 point
        !    integral estimate.  On the second call, the user sets M = 2,
        !    and the input value of Q should be the integral estimate returned
        !    on the previous call.  By incrementing M and passing the previous
        !    estimate, the user gets a sequence of increasingly improved
        !    integral estimates:
        !
        !    q = 0.0
        !    m = 1
        !    do
        !      q_new = trapezoid_refine ( a, b, m, f, q )
        !      if ( satisfied ) then
        !        exit
        !      end if
        !      q = q_new
        !      m = m + 1
        !    end do
        !
        !    The number of points used on each step of the iteration is:
        !
        !    M   N
        !    1   2
        !    2   3
        !    3   5
        !    4   9
        !    5  17
        !    6  33
        !    m   2^(m-1)+1
        !
        !  Licensing:
        !
        !    This code is distributed under the GNU LGPL license.
        !
        !  Modified:
        !
        !    27 February 2013
        !
        !  Author:
        !
        !    John Burkardt
        !
        !  Parameters:
        !
        !    Input, real  A, B, the endpoints of the interval.
        !
        !    Input, integer  M, the current level of refinement.
        !    On the first call, M should be 1.  On each subsequent call,
        !    the user should increment M by 1.
        !
        !    Input, external real  F, the name of the function
        !    which evaluates the integrand, and whose form is:
        !      function f ( x )
        !      real  f
        !      real  x
        !      f = ...
        !
        !    Input, real  Q, the integral estimate return on
        !    the previous call.  But on the first call, set Q to zero.
        !
        !    Output, real  TRAPEZOID_REFINE, the improved
        !    integral estimate.
        !
        implicit none

        real  a
        real b
        !real , external :: f
        integer  i
        integer k
        integer  m
        real  q
        real  trapezoid_refine
        real  value
        real  x
        TYPE(ParamsCommonInf) :: param

        interface

            real function f(theta, param)

                USE COM_VAR
                TYPE(ParamsCommonInf) :: param
                REAL, INTENT(IN) :: theta

            end function f


        end interface

        if ( m < 1 ) then

            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'TRAPEZOID_REFINE - Fatal error!'
            write ( *, '(a)' ) '  Illegal input value of M.'
            stop

        else if ( m == 1 ) then

            value = ( b - a ) / 2.0D+00 * ( f ( a, param ) + f ( b, param ) )

        else if ( 2 <= m ) then

            k = 2 ** ( m - 2 )
            value = 0.0D+00
            do i = 1, k
                x = ( real ( 2 * k - 2 * i + 1 ) * a   &
                    + real (         2 * i - 1 ) * b ) &
                    / real ( 2 * k )
                value = value + f ( x, param )
            end do

            value = 0.5D+00 * q + ( b - a ) * value / real ( 2 * k )

        end if

        trapezoid_refine = value

        return
    end function trapezoid_refine





    complex function romberg_trap_complex ( f, a, b, tol, maxM, param )

        !*****************************************************************************80
        !
        !! ROMBERG_TRAP approximates an integral by extrapolating the trapezoidal rule.
        !
        !  Discussion:
        !
        !    This routine computes a sequence of integral estimates involving the
        !    trapezoid rule.  Extrapolation of successive trapezoidal estimates
        !    produces an estimate with a higher rate of convergence.
        !
        !  Licensing:
        !
        !    This code is distributed under the GNU LGPL license.
        !
        !  Modified:
        !
        !    28 February 2013
        !
        !  Author:
        !
        !    John Burkardt
        !
        !  Parameters:
        !
        !    Input, real  A, B, the endpoints of the interval.
        !
        !    Input, external real  F, the name of the function
        !    which evaluates the integrand, and whose form is:
        !      function f ( x )
        !      real  f
        !      real  x
        !      f = ...
        !
        !    Input, real  TOL, the tolerance.  When two successive integral
        !    estimates differ by less than this value, the iteration will halt.
        !
        !    Output, integer  M, the number of trapezoidal estimates
        !    that were required.
        !
        !    Output, real  ROMBERG_TRAP, the integral estimate.
        !
        implicit none

        real  a
        real  b
        complex, external :: f
        integer  m, maxM
        complex q
        complex q_new
        complex r
        complex r_new
        real  tol
        TYPE(ParamsCommonInf) :: param
        COMPLEX,PARAMETER :: ZERO=CMPLX(0.,0.)

        q_new = ZERO
        r_new = ZERO
        m = 1

        do

            q = q_new
            q_new = trapezoid_refine_complex ( a, b, m, f, q, param )

            if ( m == 1 ) then
                r_new = q_new
            else
                r = r_new
                r_new = ( 4.0D+00 * q_new - q ) / 3.0D+00
                if ( abs ( r_new - r ) .lt. tol * ( 1.0 + abs ( r_new ) ) ) then
                    exit
                end if
            end if

            if ( maxM <= m ) then
                !write ( *, '(a)' ) ' '
                !write ( *, '(a)' ) 'ROMBERG_TRAP!'
                !write ( *, '(a)' ) '  No convergence in 20 iterations.'
                !write ( *, '(a)' ) '  The algorithm is halting.'
                return
            end if

            m = m + 1

        end do

        romberg_trap_complex = r_new

        return
    end function romberg_trap_complex

    complex function trapezoid_refine_complex ( a, b, m, f, q, param)

        !*****************************************************************************80
        !
        !! TRAPEZOID_REFINE carries out a step of trapezoidal refinement.
        !
        !  Discussion:
        !
        !    This routine is designed to be an efficient way to carry out a
        !    sequence of integral estimates, using the trapezoidal rule
        !    and a nested sequence of evaluation points, in such a way that
        !    the function is not re-evaluated unnecessarily.
        !
        !    The user calls first with M = 1 and Q = 0 to get a 2 point
        !    integral estimate.  On the second call, the user sets M = 2,
        !    and the input value of Q should be the integral estimate returned
        !    on the previous call.  By incrementing M and passing the previous
        !    estimate, the user gets a sequence of increasingly improved
        !    integral estimates:
        !
        !    q = 0.0
        !    m = 1
        !    do
        !      q_new = trapezoid_refine ( a, b, m, f, q )
        !      if ( satisfied ) then
        !        exit
        !      end if
        !      q = q_new
        !      m = m + 1
        !    end do
        !
        !    The number of points used on each step of the iteration is:
        !
        !    M   N
        !    1   2
        !    2   3
        !    3   5
        !    4   9
        !    5  17
        !    6  33
        !    m   2^(m-1)+1
        !
        !  Licensing:
        !
        !    This code is distributed under the GNU LGPL license.
        !
        !  Modified:
        !
        !    27 February 2013
        !
        !  Author:
        !
        !    John Burkardt
        !
        !  Parameters:
        !
        !    Input, real  A, B, the endpoints of the interval.
        !
        !    Input, integer  M, the current level of refinement.
        !    On the first call, M should be 1.  On each subsequent call,
        !    the user should increment M by 1.
        !
        !    Input, external real  F, the name of the function
        !    which evaluates the integrand, and whose form is:
        !      function f ( x )
        !      real  f
        !      real  x
        !      f = ...
        !
        !    Input, real  Q, the integral estimate return on
        !    the previous call.  But on the first call, set Q to zero.
        !
        !    Output, real  TRAPEZOID_REFINE, the improved
        !    integral estimate.
        !
        implicit none

        real  a
        real  b
        complex, external :: f
        integer  i
        integer  k
        integer  m
        complex q
        complex value
        real  x
        TYPE(ParamsCommonInf) :: param
        COMPLEX,PARAMETER :: ZERO=CMPLX(0.,0.)

        if ( m < 1 ) then

            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'TRAPEZOID_REFINE - Fatal error!'
            write ( *, '(a)' ) '  Illegal input value of M.'
            stop

        else if ( m == 1 ) then

            value = ( b - a ) / 2.0D+00 * ( f ( a, param ) + f ( b, param ) )

        else if ( 2 <= m ) then

            k = 2 ** ( m - 2 )
            value = ZERO
            do i = 1, k
                x = ( real ( 2 * k - 2 * i + 1) * a   &
                    + real (         2 * i - 1) * b ) &
                    / real ( 2 * k)
                value = value + f ( x, param )
            end do

            value = 0.5D+00 * q + ( b - a ) * value / real ( 2 * k)

        end if

        trapezoid_refine_complex = value

        return
    end function trapezoid_refine_complex


    SUBROUTINE COMPUTE_NEXP(param)
    ! This is a utility to compute precompute values needed to integrate green function
    ! and to compute next
    ! param: common variables that contains useful variables
    !        See param definition for explanation of the member used here

    USE COMPUTE_GREEN_FD
    USE COMPUTE_GREEN_INFD
    USE COM_VAR

        REAL:: GM, BJ, DIJ, AKK
        TYPE(ParamsCommonInf), TARGET :: param
        INTEGER:: NP1, N,  JJ, I, J, I1, NJJ
        INTEGER, POINTER:: NQ
        REAL, POINTER:: CQ(:),QQ(:),AMBDA(:),AR(:)
        REAL, POINTER :: FSP,FSM,VSXP,VSYP,VSZP,VSXM,VSYM,VSZM
        TYPE(TempVar), POINTER :: SolverVar


        SolverVar => param%SolverVar



        AMBDA => SolverVar%AMBDA
        AR => SolverVar%AR
        NQ => SolverVar%NQ
        CQ => SolverVar%CQ
        QQ => SolverVar%QQ
        VSXP => SolverVar%VSXP
        VSYP => SolverVar%VSYP
        VSZP => SolverVar%VSZP
        VSXM => SolverVar%VSXM
        VSYM => SolverVar%VSYM
        VSZM => SolverVar%VSZM

         !--------------Initilizations---------------
        VSXP=0.
        VSXM=0.
        VSYP=0.
        VSYM=0.
        VSZP=0.
        VSZM=0.
        !--------------------------------------------

        NJJ=NSYMY+1


        IF(ProblemSavedAt(SolverVar%ProblemNumber) < 0 ) THEN
            CQ=0.0
            QQ=0.0
            AMBDA=0.0
            AR=0.0

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
            AKK=param%AM0*GM

        END IF


        IF(param%is_infinite == 1) THEN

            IF(ProblemSavedAt(SolverVar%ProblemNumber) < 0 ) THEN

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

        ELSE
            IF(ProblemSavedAt(SolverVar%ProblemNumber) < 0 ) THEN

                CALL CINT_FD(AKK,N)
                !N=NQ Bug Fixes
                NQ = N
                CALL LISC(param%AKH,param%AMH,param%NEXP, SolverVar)

                IF(ProblemSavedAt(SolverVar%ProblemNumber) == -1 ) THEN
                    DO I=1,31
                        ARCache(ProblemToSaveLocation(SolverVar%ProblemNumber), I)  = AR(I)
                        AMBDACache(ProblemToSaveLocation(SolverVar%ProblemNumber), I) =  AMBDA(I)

                    END DO
                    NEXPCache(ProblemToSaveLocation(SolverVar%ProblemNumber)) = param%NEXP

                END IF

            ELSE

                DO I=1,31
                    AR(I) = ARCache(ProblemSavedAt(SolverVar%ProblemNumber), I)
                    AMBDA(I) = AMBDACache(ProblemSavedAt(SolverVar%ProblemNumber), I)
                END DO

                param%NEXP = NEXPCache(ProblemSavedAt(SolverVar%ProblemNumber))

            END IF

        END IF

    END SUBROUTINE COMPUTE_NEXP


end module UTILITY
