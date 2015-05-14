!--------------------------------------------------------------------------------------
!
!Copyright (C) 2014 TopCoder Inc., All Rights Reserved.
!
!--------------------------------------------------------------------------------------

!*
!*  Copyright (C) CERFACS 1998
!*
!*  SOFTWARE LICENSE AGREEMENT NOTICE - THIS SOFTWARE IS BEING PROVIDED TO
!*  YOU BY CERFACS UNDER THE FOLLOWING LICENSE. BY DOWN-LOADING, INSTALLING
!*  AND/OR USING THE SOFTWARE YOU AGREE THAT YOU HAVE READ, UNDERSTOOD AND
!*  WILL COMPLY WITH THESE FOLLOWING TERMS AND CONDITIONS.
!*
!*  1 - This software program provided in source code format ("the " Source
!*  Code ") and any associated documentation (the " Documentation ") are
!*  licensed, not sold, to you.
!*
!*  2 - CERFACS grants you a personal, non-exclusive, non-transferable and
!*  royalty-free right to use, copy or modify the Source Code and
!*  Documentation, provided that you agree to comply with the terms and
!*  restrictions of this agreement. You may modify the Source Code and
!*  Documentation to make source code derivative works, object code
!*  derivative works and/or documentation derivative Works (called "
!*  Derivative Works "). The Source Code, Documentation and Derivative
!*  Works (called " Licensed Software ") may be used by you for personal
!*  and non-commercial use only. " non-commercial use " means uses that are
!*  not or will not result in the sale, lease or rental of the Licensed
!*  Software and/or the use of the Licensed Software in any commercial
!*  product or service. CERFACS reserves all rights not expressly granted
!*  to you. No other licenses are granted or implied.
!*
!*  3 - The Source Code and Documentation are and will remain the sole
!*  property of CERFACS. The Source Code and Documentation are copyrighted
!*  works. You agree to treat any modification or derivative work of the
!*  Licensed Software as if it were part of the Licensed Software itself.
!*  In return for this license, you grant CERFACS a non-exclusive perpetual
!*  paid-up royalty-free license to make, sell, have made, copy, distribute
!*  and make derivative works of any modification or derivative work you
!*  make of the Licensed Software.
!*
!*  4- The licensee shall acknowledge the contribution of the Source Code
!*  (using the reference [1]) in any publication of material dependent upon
!*  upon the use of the Source Code. The licensee shall use reasonable
!*  endeavours to notify the authors of the package of this publication.
!*
!*  [1] V. Frayssé, L. Giraud, S. Gratton, and J. Langou, A set of GMRES
!*    routines for real and complex arithmetics on high performance
!*    computers, CERFACS Technical Report TR/PA/03/3, public domain software
!*    available on www.cerfacs/algor/Softs, 2003
!*
!*  5- CERFACS has no obligation to support the Licensed Software it is
!*  providing under this license.
!*
!*  THE LICENSED SOFTWARE IS PROVIDED " AS IS " AND CERFACS MAKE NO
!*  REPRESENTATIONS OR WARRANTIES, EXPRESS OR IMPLIED. BY WAY OF EXAMPLE,
!*  BUT NOT LIMITATION, CERFACS MAKE NO REPRESENTATIONS OR WARRANTIES OF
!*  MERCHANTIBILY OR FITNESS FOR ANY PARTICULAR PURPOSE OR THAT THE USE OF
!*  THE LICENSED SOFTWARE OR DOCUMENTATION WILL NOT INFRINGE ANY THIRD
!*  PARTY PATENTS, COPYRIGHTS, TRADEMARKS OR OTHER RIGHTS. CERFACS WILL NOT
!*  BE LIABLE FOR ANY CONSEQUENTIAL, INCIDENTAL, OR SPECIAL DAMAGES, OR ANY
!*  OTHER RELIEF, OR FOR ANY CLAIM BY ANY THIRD PARTY, ARISING FROM YOUR
!*  USE OF THE LICENSED SOFTWARE.
!*
!*  6- For information regarding a commercial license for the Source Code
!*  and Documentation, please contact Mrs Campassens (campasse@cerfacs.fr)
!*
!*  7- This license is effective until terminated. You may terminate this
!*  license at any time by destroying the Licensed Software.
!*
!*    I agree all the terms and conditions of the above license agreement
!*


!   This module contains an iterative solver to find a solution for a system of complex linear equations
!   It implements the GMRES algorithm.
!
!   @author TCSASSEMBLER
!   @version 1.1
MODULE GMRES

! This type contains variables that need to be saved between gmres iteration.
! It will be passed to every subroutine of this module. It was introduced so that
! gmres can be used independently to solve multiple equation in same program
TYPE GMRESVAR

integer :: icheck
integer :: iterOut, jH, retlbl, j, nOrtho, compRsd
real ::  beta, bn, dnormres, sA, sb, sPA, sPb, dnormx, trueNormRes, bea, be, dloo

END TYPE GMRESVAR

CONTAINS

    subroutine drive_cgmres(n,nloc,m,lwork,work, &
    irc,icntl,cntl,info,rinfo, env)

!  Purpose
!  =======
!    drive_cgmres is the driver routine for solving the linear system
!  Ax = b using the *  Generalized Minimal Residual iterative method
!  with preconditioning.
!  This solver is implemented with a reverse communication scheme: control
!  is returned to the user for computing the (preconditioned)
!  matrix-vector product.
!  See the User's Guide for an example of use.


! Written : June 1996
! Authors : Luc Giraud, Serge Gratton, V. Fraysse
!             Parallel Algorithms - CERFACS

! Updated : April 1997
! Authors :  Valerie Fraysse, Luc Giraud, Serge Gratton
!             Parallel Algorithms - CERFACS

! Updated : June 1998
! Authors : Valerie Fraysse, Luc Giraud, Serge Gratton
!             Parallel Algorithms - CERFACS
! Purpose : Make clear that the warning and error messages come from the
!           cgmres modules.

! Updated : December 2002 - L. Giraud, J.Langou
! Purpose : Add the capability to avoid explicit residual calculation at restart

! Updated : March 2005 - L. Giraud
! Purpose : small bug in the computation of the restart if too large vs the
!           size of the workspace

!  Arguments
!  =========

!  n      (input) INTEGER.
!          On entry, the dimension of the problem.
!          Unchanged on exit.

!  nloc   (input) INTEGER.
!          On entry, the dimension of the local problem.
!          In a parallel distributed envirionment, this corresponds
!          to the size of the subset of entries of the right hand side
!          and solution allocated to the calling process.
!          Unchanged on exit.


!  m      (input) INTEGER
!          Restart parameter, <= N. This parameter controls the amount
!          of memory required for matrix H (see WORK and H).
!          Unchanged on exit.

!  lwork  (input) INTEGER
!          size of the workspace
!          if (icntl(5) = 0 or 1 )
!            lwork >= m*m + m*(n+5) + 5*n+2, if icntl(8) = 1
!            lwork >= m*m + m*(n+5) + 6*n+2, if icntl(8) = 0
!          if (icntl(5) = 2 or 3 )
!            lwork >= m*m + m*(n+5) + 5*n+m+1, if icntl(8) = 1
!            lwork >= m*m + m*(n+5) + 6*n+m+1, if icntl(8) = 0

!  work   (workspace) real/complex array, length lwork
!          work contains the required vector and matrices stored in the
!          following order :
!            x  (n,1)       : computed solution.
!            b  (n,1)       : right hand side.
!            r0 (n,1)       : vector workspace.
!            w  (n,1)       : vector workspace.
!            V  (n,m)       : Krylov basis.
!            H  (m+1,m+1)   : Hessenberg matrix (full storage).
!            yCurrent (m,1) : solution of the current LS
!            xCurrent (n,1) : current iterate
!            rotSin (m,1)   : Sine of the Givens rotation
!            rotCos (m,1)   : Cosine of the Givens rotation

!  irc     (input/output) INTEGER array. length 5
!            irc(1) : REVCOM   used for reverse communication
!                             (type of external operation)
!            irc(2) : COLX     used for reverse communication
!            irc(3) : COLY     used for reverse communication
!            irc(4) : COLZ     used for reverse communication
!            irc(5) : NBSCAL   used for reverse communication

!  icntl   (input) INTEGER array. length 7
!            icntl(1) : stdout for error messages
!            icntl(2) : stdout for warnings
!            icntl(3) : stdout for convergence history
!            icntl(4) : 0 - no preconditioning
!                       1 - left preconditioning
!                       2 - right preconditioning
!                       3 - double side preconditioning
!                       4 - error, default set in Init
!            icntl(5) : 0 - modified Gram-Schmidt
!                       1 - iterative modified Gram-Schmidt
!                       2 - classical Gram-Schmidt
!                       3 - iterative classical Gram-Schmidt
!            icntl(6) : 0 - default initial guess x_0 = 0 (to be set)
!                       1 - user supplied initial guess
!            icntl(7) : maximum number of iterations
!            icntl(8) : 0 - use recurence formula at restart
!                       1 - default compute the true residual at each restart

!  cntl    (input) real array, length 5
!            cntl(1) : tolerance for convergence
!            cntl(2) : scaling factor for normwise perturbation on A
!            cntl(3) : scaling factor for normwise perturbation on b
!            cntl(4) : scaling factor for normwise perturbation on the
!                      preconditioned matrix
!            cntl(5) : scaling factor for normwise perturbation on
!                      preconditioned right hand side

!  info    (output) INTEGER array, length 3
!            info(1) :  0 - normal exit
!                      -1 - n < 1
!                      -2 - m < 1
!                      -3 - lwork too small
!                      -4 - convergence not achieved after icntl(7) iterations
!                      -5 - precondition type not set by user
!            info(2) : if info(1)=0 - number of iteration to converge
!                      if info(1)=-3 - minimum workspace size necessary
!            info(3) : optimal size for the workspace

!  rinfo   (output) real array, length 2
!            if info(1)=0
!              rinfo(1) : backward error for the preconditioned system
!              rinfo(2) : backward error for the unpreconditioned system

! Input variables
! ---------------
    integer :: n, nloc, lwork, icntl(*)
    real ::   cntl(*)
    real ::   sA, sb, sPA, sPb
! Output variables
! ----------------
    integer ::  info(*)
    real ::    rinfo(*)
! Input/Output variables
! ----------------------
    integer ::  m, irc(*)
    complex :: work(*)
    TYPE(GMRESVAR), TARGET :: env
! Local variables
! ---------------
    integer :: xptr, bptr, wptr, r0ptr, Vptr, Hptr, dotptr
    integer :: yCurrent,rotSin, rotCos, xCurrent
    integer :: sizeWrk, newRestart
    integer :: iwarn, ierr, ihist, compRsd
    real ::    rn,rx,rc
    real :: DZRO
    parameter (DZRO = 0.0e0)

    integer, pointer :: icheck
    !save icheck
    !DATA icheck /0/

    intrinsic ifix, float

!       Executable statements :

    icheck => env%icheck

    ierr    = icntl(1)
    iwarn   = icntl(2)
    ihist   = icntl(3)
    compRsd = icntl(8)

    if (ierr < 0) ierr = 6

    if (compRsd == 1) then
        sizeWrk  = m*m + m*(nloc+5) + 5*nloc+ 1
    else
        sizeWrk  = m*m + m*(nloc+5) + 6*nloc+ 1
    endif

    if (icheck == 0) then
    ! Check the value of the arguments
        if ((n < 1) .OR. (nloc < 1)) then
            write(ierr,*)
            write(ierr,*)' ERROR GMRES : '
            write(ierr,*)'     N < 1 '
            write(ierr,*)
            info(1) = -1
            irc(1)  = 0
            return
        endif
        if (m < 1) then
            write(ierr,*)
            write(ierr,*)' ERROR GMRES :'
            write(ierr,*)'     M < 1 '
            write(ierr,*)
            info(1) = -2
            irc(1)  = 0
            return
        endif
        if ((icntl(4) /= 0) .AND. (icntl(4) /= 1) .AND. &
        (icntl(4) /= 2) .AND. (icntl(4) /= 3)) then
            write(ierr,*)
            write(ierr,*)' ERROR GMRES : '
            write(ierr,*)'     Undefined preconditioner '
            write(ierr,*)
            info(1) = -5
            irc(1)  = 0
            return
        endif

        if ((icntl(5) < 0) .OR. (icntl(5) > 3)) then
            icntl(5) = 0
            if (iwarn /= 0) then
                write(iwarn,*)
                write(iwarn,*) ' WARNING  GMRES : '
                write(iwarn,*) '       Undefined orthogonalisation '
                write(iwarn,*) '       Default MGS '
                write(iwarn,*)
            endif
        endif

        if ((icntl(5) == 2) .OR. (icntl(5) == 3)) then
        ! the workspace should be large enough to store the m dot-products
            sizeWrk  = sizeWrk  + m
        else
            sizeWrk  = sizeWrk  + 1
        endif

        if (iwarn /= 0) then
            write(iwarn,*)
            write(iwarn,*) ' WARNING GMRES : '
            write(iwarn,*) '       For M = ',m,' optimal value '
            write(iwarn,*) '       for LWORK =  ', sizeWrk
            write(iwarn,*)
        endif

        if ((icntl(6) /= 0) .AND. (icntl(6) /= 1)) then
            icntl(6) = 0
            if (iwarn /= 0) then
                write(iwarn,*)
                write(iwarn,*) ' WARNING GMRES : '
                write(iwarn,*) '       Undefined intial guess '
                write(iwarn,*) '       Default x0 = 0 '
                write(iwarn,*)
            endif
        endif
        if (icntl(7) <= 0) then
            icntl(7) = n
            if (iwarn /= 0) then
                write(iwarn,*)
                write(iwarn,*) ' WARNING GMRES :'
                write(iwarn,*) '       Negative max number of iterations'
                write(iwarn,*) '       Default N '
                write(iwarn,*)
            endif
        endif
        if ((icntl(8) /= 0) .AND. (icntl(8) /= 1)) then
            icntl(8) = 1
            write(iwarn,*)
            write(iwarn,*) ' WARNING GMRES :'
            write(iwarn,*) '       Undefined strategy for the residual'
            write(iwarn,*) '       at restart'
            write(iwarn,*) '       Default 1 '
            write(iwarn,*)
        endif
    ! Check if the restart parameter is correct and if the size of the
    !  workspace is big enough for the restart.
    ! If not try to fix correctly the parameters

        if ((m > n) .OR. (lwork < sizeWrk)) then
            if (m > n) then
                m = n
                if (iwarn /= 0) then
                    write(iwarn,*)
                    write(iwarn,*) ' WARNING GMRES : '
                    write(iwarn,*) '       Parameter M bigger than N'
                    write(iwarn,*) '       New value for M ',m
                    write(iwarn,*)
                endif
                if (compRsd == 1) then
                    sizeWrk = m*m + m*(nloc+5) + 5*nloc+1
                else
                    sizeWrk = m*m + m*(nloc+5) + 6*nloc+1
                endif
                if ((icntl(5) == 2) .OR. (icntl(5) == 3)) then
                ! the workspace should be large enough to store the m dot-products
                    sizeWrk  = sizeWrk  + m
                else
                    sizeWrk  = sizeWrk  + 1
                endif
            endif
            if ((lwork < sizeWrk) .AND. (n == nloc)) then
            ! Compute the maximum size of the m according to the memory space
                rn         = float(n)
                rx         = rn + 5.0
                rc         = 5.0*rn + 1 - float(lwork)

            ! Update the linear part of the second order equation to be solved
                if ((icntl(5) == 2) .OR. (icntl(5) == 3)) then
                    rx = rx + 1
                endif
            ! Update the constant part of the second order equation to be solved

                if (icntl(8) == 0) then
                    rc = rc + rn
                endif
            ! Solution of the the second order equation
                newRestart = ifix((-rx+sqrt(rx**2-4.0*rc))/2.0)
                if (newRestart > 0) then
                    m = newRestart
                    if (iwarn /= 0) then
                        write(iwarn,*)
                        write(iwarn,*)' WARNING GMRES : '
                        write(iwarn,*)'       Workspace too small for M'
                        write(iwarn,*)'       New value for M ',m
                        write(iwarn,*)
                    endif
                else
                    write(ierr,*)
                    write(ierr,*)' ERROR GMRES : '
                    write(ierr,*)'     Not enough space for the problem'
                    write(ierr,*)'     the space does not permit any m'
                    write(ierr,*)
                    info(1) = -3
                    irc(1)  = 0
                    return
                endif
            endif
            if ((lwork < sizeWrk) .AND. (n /= nloc)) then
                write(ierr,*)
                write(ierr,*)' ERROR GMRES : '
                write(ierr,*)'     Not enough space for the problem'
                write(ierr,*)
                info(1) = -3
                irc(1)  = 0
                return
            endif
        endif

        info(3) = sizeWrk
        icheck = 1

    ! save the parameters the the history file

        if (ihist /= 0) then
            write(ihist,'(10x,A39)') 'CONVERGENCE HISTORY FOR GMRES'
            write(ihist,*)
            write(ihist,'(A30,I2)') 'Errors are displayed in unit: ',ierr
            if (iwarn == 0) then
                write(ihist,'(A27)') 'Warnings are not displayed:'
            else
                write(ihist,'(A32,I2)') 'Warnings are displayed in unit: ', &
                iwarn
            endif
            write(ihist,'(A13,I7)') 'Matrix size: ',n
            write(ihist,'(A19,I7)') 'Local matrix size: ',nloc
            write(ihist,'(A9,I7)') 'Restart: ',m
            if (icntl(4) == 0) then
                write(ihist,'(A18)') 'No preconditioning'
            elseif (icntl(4) == 1) then
                write(ihist,'(A20)') 'Left preconditioning'
            elseif (icntl(4) == 2) then
                write(ihist,'(A21)') 'Right preconditioning'
            elseif (icntl(4) == 3) then
                write(ihist,'(A30)') 'Left and right preconditioning'
            endif
            if (icntl(5) == 0) then
                write(ihist,'(A21)') 'Modified Gram-Schmidt'
            elseif (icntl(5) == 1) then
                write(ihist,'(A31)') 'Iterative modified Gram-Schmidt'
            elseif (icntl(5) == 2) then
                write(ihist,'(A22)') 'Classical Gram-Schmidt'
            else
                write(ihist,'(A32)') 'Iterative classical Gram-Schmidt'
            endif
            if (icntl(6) == 0) then
                write(ihist,'(A29)') 'Default initial guess x_0 = 0'
            else
                write(ihist,'(A27)') 'User supplied initial guess'
            endif
            if (icntl(8) == 1) then
                write(ihist,'(A33)') 'True residual computed at restart'
            else
                write(ihist,'(A30)') 'Recurrence residual at restart'
            endif
            write(ihist,'(A30,I5)') 'Maximum number of iterations: ', &
            icntl(7)
            write(ihist,'(A27,E9.2)') 'Tolerance for convergence: ', &
            cntl(1)

            write(ihist,'(A53)') &
            'Backward error on the unpreconditioned system Ax = b:'
            sA       = cntl(2)
            sb       = cntl(3)
            if ((sA == DZRO) .AND. (sb == DZRO)) then
                write(ihist,'(A39)') &
                '    the residual is normalised by ||b||'
            else
                write(ihist,1) sA,sb
                1 format('    the residual is normalised by         ',E9.2, &
                ' * ||x|| + ',E9.2)
            endif
            sPA      = cntl(4)
            sPb      = cntl(5)
            write(ihist,2)
            2 format('Backward error on the preconditioned system', &
            ' (P1)A(P2)y = (P1)b:')
            if ((sPA == DZRO) .AND. (sPb == DZRO)) then
                write(ihist,3)
                3 format('    the preconditioned residual is normalised ', &
                'by ||(P1)b||')
            else
                write(ihist,4) sPA, sPb
                4 format('    the preconditioned residual is normalised by ', E9.2, &
                ' * ||(P2)y|| + ',E9.2)
            endif

            write(ihist,5) info(3)
            5 format('Optimal size for the local workspace:',I7)
            write(ihist,*)
            write(ihist,6)
            6 format('Convergence history: b.e. on the preconditioned system')
            write(ihist,7)
            7 format(' Iteration   Arnoldi b.e.    True b.e.')
        endif

    endif
! setup some pointers on the workspace
    xptr     = 1
    bptr     = xptr + nloc
    r0ptr    = bptr + nloc
    wptr     = r0ptr + nloc
    Vptr     = wptr + nloc
    if (compRsd == 1) then
        Hptr     = Vptr + m*nloc
    else
        Hptr     = Vptr + (m+1)*nloc
    endif
    dotptr   = Hptr + (m+1)*(m+1)
    if ((icntl(5) == 2) .OR. (icntl(5) == 3)) then
        yCurrent = dotptr + m
    else
        yCurrent = dotptr + 1
    endif
    xCurrent = yCurrent + m
    rotSin   = xCurrent + nloc
    rotCos   = rotSin + m

    call cgmres(nloc,m,work(bptr),work(xptr), &
    work(Hptr),work(wptr),work(r0ptr),work(Vptr), &
    work(dotptr),work(yCurrent),work(xCurrent), &
    work(rotSin),work(rotCos),irc,icntl,cntl,info,rinfo, env)

    if (irc(1) == 0) then
        icheck = 0
    endif

    return
    end subroutine drive_cgmres

    subroutine cgmres(n,m,b,x,H,w,r0,V,dot,yCurrent,xCurrent,rotSin, &
    rotCos,irc,icntl,cntl,info,rinfo, env)


!  Purpose
!  =======
!  cgmres solves the linear system Ax = b using the
!  Generalized Minimal Residual iterative method

! When preconditioning is used we solve :
!     M_1^{-1} A M_2^{-1} y = M_1^{-1} b
!     x = M_2^{-1} y

!   Convergence test based on the normwise backward error for
!  the preconditioned system

! Written : June 1996
! Authors : Luc Giraud, Serge Gratton, V. Fraysse
!             Parallel Algorithms - CERFACS

! Updated : April 1997
! Authors :  Valerie Fraysse, Luc Giraud, Serge Gratton
!             Parallel Algorithms - CERFACS

! Updated : March 1998
! Purpose : Pb with F90 on DEC ws
!           cure : remove "ZDSCAL" when used to initialize vectors to zero

! Updated : May 1998
! Purpose : r0(1) <-- r0'r0 : pb when used with DGEMV for the dot product
!           cure : w(1) <--  r0'r0

! Updated : June 1998
! Purpose : Make clear that the warning and error messages come from the
!           cgmres modules.

! Updated : February 2001 - L. Giraud
! Purpose : In complex version, initializations to zero performed  in complex
!           arithmetic to avoid implicit conversion by the compiler.

! Updated : July 2001 - L. Giraud, J. Langou
! Purpose : Avoid to compute the approximate solution at each step of
!           the Krylov space construction when spA is zero.

! Updated : November 2002 - S. Gratton
! Purpose : Use Givens rotations conform to the classical definition.
!           No impact one the convergence history.

! Updated : November 2002 - L. Giraud
! Purpose : Properly handle the situation when the convergence is obtained
!           exactly at the "IterMax" iteration

! Updated : December 2002 - L. Giraud, J.Langou
! Purpose : Add the capability to avoid explicit residual calculation at restart

! Updated : January  2003 - L. Giraud, S. Gratton
! Purpose : Use Givens rotations from BLAS.

! Updated : March    2003 - L. Giraud
! Purpose : Set back retlbl to zero, if initial guess is solution
!           or right-hand side is zero

! Updated : September 2003 - L. Giraud
! Purpose : Include room in the workspace to store the results of the dot products
!           Fix the bugs that appeared when M > Nloc

!  Arguments
!  =========

!  n       (input) INTEGER.
!           On entry, the dimension of the problem.
!           Unchanged on exit.

!  m        (input) INTEGER
!           Restart parameter, <= N. This parameter controls the amount
!           of memory required for matrix H (see WORK and H).
!           Unchanged on exit.

!  b        (input) real/complex
!           Right hand side of the linear system.

!  x        (output) real/complex
!           Computed solution of the linear system.

!  H        (workspace)  real/complex
!           Hessenberg matrix built within dgmres

!  w        (workspace)  real/complex
!           Vector used as temporary storage

!  r0       (workspace)  real/complex
!           Vector used as temporary storage

!  V        (workspace)  real/complex
!           Basis computed by the Arnoldi's procedure.

!  dot      (workspace) real/complex
!           Store the results of the dot product calculation

!  yCurrent (workspace) real/complex
!           solution of the current LS

!  xCurrent (workspace) real/complex
!           current iterate

!  rotSin   (workspace) real/complex
!           Sine of the Givens rotation

!  rotCos   (workspace) real
!           Cosine of the Givens rotation

!  irc      (input/output) INTEGER array. length 3
!             irc(1) : REVCOM   used for reverse communication
!                              (type of external operation)
!             irc(2) : COLX     used for reverse communication
!             irc(3) : COLY     used for reverse communication
!             irc(4) : COLZ     used for reverse communication
!             irc(5) : NBSCAL   used for reverse communication

!  icntl    (input) INTEGER array. length 7
!             icntl(1) : stdout for error messages
!             icntl(2) : stdout for warnings
!             icntl(3) : stdout for convergence history
!             icntl(4) : 0 - no preconditioning
!                        1 - left preconditioning
!                        2 - right preconditioning
!                        3 - double side preconditioning
!                        4 - error, default set in Init
!             icntl(5) : 0 - modified Gram-Schmidt
!                        1 - iterative modified Gram-Schmidt
!                        2 - classical Gram-Schmidt
!                        3 - iterative classical Gram-Schmidt
!             icntl(6) : 0 - default initial guess x_0 = 0 (to be set)
!                        1 - user supplied initial guess
!             icntl(7) : maximum number of iterations
!             icntl(8) : 1 - default compute the true residual at each restart
!                        0 - use recurence formula at restart

!  cntl     (input) real array, length 5
!             cntl(1) : tolerance for convergence
!             cntl(2) : scaling factor for normwise perturbation on A
!             cntl(3) : scaling factor for normwise perturbation on b
!             cntl(4) : scaling factor for normwise perturbation on the
!                       preconditioned matrix
!             cntl(5) : scaling factor for normwise perturbation on
!                       preconditioned right hand side

!  info     (output) INTEGER array, length 2
!             info(1) :  0 - normal exit
!                       -1 - n < 1
!                       -2 - m < 1
!                       -3 - lwork too small
!                       -4 - convergence not achieved after icntl(7) iterations
!                       -5 - precondition type not set by user
!             info(2) : if info(1)=0 - number of iteration to converge
!                       if info(1)=-3 - minimum workspace size necessary
!             info(3) : optimal size for the workspace

! rinfo     (output) real array, length 2
!             if info(1)=0
!               rinfo(1) : backward error for the preconditioned system
!               rinfo(2) : backward error for the unpreconditioned system

! Input variables
! ---------------
    integer ::  n, m, icntl(*)
    complex :: b(*)
    real ::    cntl(*)

! Output variables
! ----------------
    integer ::  info(*)
    real ::    rinfo(*)

! Input/Output variables
! ----------------------
    integer ::  irc(*)
    complex :: x(*), H(m+1,*), w(*), r0(*)
    complex :: V(n,*),dot(*),yCurrent(*)
    complex :: xCurrent(*), rotSin(*)
    complex :: rotCos(*)
    TYPE(GMRESVAR), TARGET :: env

! Local variables
! ---------------
    integer ::  iterMax, initGuess, iOrthog
    integer ::  xptr, bptr, wptr, r0ptr, Vptr, Hptr, yptr, xcuptr
    integer ::  dotptr
    integer ::  typePrec, leftPrec, rightPrec, dblePrec, noPrec
    integer ::  iwarn, ihist
    real ::     temp
    real ::     dnormw
    complex :: dVi, aux
    complex :: auxHjj, auxHjp1j

    parameter (noPrec = 0, leftPrec = 1)
    parameter (rightPrec = 2, dblePrec = 3)

    complex :: ZERO, ONE
    parameter (ZERO = (0.0e0, 0.0e0), ONE = (1.0e0, 0.0e0))
    real :: DZRO,DONE
    parameter (DZRO = 0.0e0, DONE = 1.0e0)

    integer, pointer :: iterOut, jH, retlbl, j, nOrtho, compRsd
    real, pointer ::  beta, bn, dnormres, sA, sb, sPA, sPb, dnormx, trueNormRes, bea, be, dloo


! External functions
! ------------------
    real ::    scnrm2
    external scnrm2

! Saved variables
! ---------------
    !save iterOut, jH, beta, bn, dnormres, retlbl, j
    !save sA, sb, sPA, sPb, dnormx, trueNormRes, bea, be
    !save dloo, nOrtho, compRsd

! Intrinsic function
! ------------------
    intrinsic abs, sqrt, real, conjg


! Reverse communication variables
! -------------------------------
    integer :: matvec, precondLeft, precondRight, prosca
    parameter(matvec=1, precondLeft=2, precondRight=3, prosca=4)
    !integer :: retlbl
    !DATA retlbl /0/

!       Executable statements
iterOut => env%iterOut
jH => env%jH
retlbl => env%retlbl
j => env%j
nOrtho => env%nOrtho
compRsd => env%compRsd
beta => env%beta
bn => env%bn
dnormres => env%dnormres
sA => env%sA
sb => env%sb
sPA => env%sPA
sPb => env%sPb
dnormx => env%dnormx
trueNormRes => env%trueNormRes
bea => env%bea
be => env%be
dloo => env%dloo

! setup some pointers on the workspace
    xptr     = 1
    bptr     = xptr + n
    r0ptr    = bptr + n
    wptr     = r0ptr + n
    Vptr     = wptr + n
    if (icntl(8) == 1) then
        Hptr     = Vptr + m*n
    else
        Hptr     = Vptr + (m+1)*n
    endif
    dotptr   = Hptr + (m+1)*(m+1)
    if ((icntl(5) == 2) .OR. (icntl(5) == 3)) then
        yptr = dotptr + m
    else
        yptr = dotptr + 1
    endif
    xcuptr   = yptr + m

    iwarn      = icntl(2)
    ihist      = icntl(3)
    typePrec   = icntl(4)
    iOrthog    = icntl(5)
    initGuess  = icntl(6)
    iterMax    = icntl(7)

    if (retlbl == 0) then
        compRsd    = icntl(8)
    endif

    if (retlbl /= 0) then
        if (retlbl == 5) then
            goto 5
        else if (retlbl == 6) then
            goto 6
        else if (retlbl == 8) then
            goto 8
        else if (retlbl == 11) then
            goto 11
        else if (retlbl == 16) then
            goto 16
        else if (retlbl == 18) then
            goto 18
        else if (retlbl == 21) then
            goto 21
        else if (retlbl == 26) then
            goto 26
        else if (retlbl == 31) then
            goto 31
        else if (retlbl == 32) then
            goto 32
        else if (retlbl == 33) then
            goto 33
        else if (retlbl == 34) then
            goto 34
        else if (retlbl == 36) then
            goto 36
        else if (retlbl == 37) then
            goto 37
        else if (retlbl == 38) then
            goto 38
        else if (retlbl == 41) then
            goto 41
        else if (retlbl == 43) then
            goto 43
        else if (retlbl == 46) then
            goto 46
        else if (retlbl == 48) then
            goto 48
        else if (retlbl == 51) then
            goto 51
        else if (retlbl == 52) then
            goto 52
        else if (retlbl == 61) then
            goto 61
        else if (retlbl == 66) then
            goto 66
        else if (retlbl == 68) then
            goto 68
        endif
    endif


! intialization of various variables

    iterOut  = 0
    beta     = DZRO

    if (initGuess == 0) then
        do j=1,n
            x(j) = ZERO
        enddo
    endif

!        bn = scnrm2(n,b,1)

    irc(1) = prosca
    irc(2) = bptr
    irc(3) = bptr
    irc(4) = dotptr
    irc(5) = 1
    retlbl = 5
    return
    5 continue
    bn = sqrt(real(dot(1)))

    if (bn == DZRO) then
        do j=1,n
            x(j) = ZERO
        enddo
        if (iwarn /= 0) then
            write(iwarn,*)
            write(iwarn,*) ' WARNING GMRES : '
            write(iwarn,*) '       Null right hand side'
            write(iwarn,*) '       solution set to zero'
            write(iwarn,*)
        endif
        jH = 0
        bea = DZRO
        be  = DZRO
        if(ihist /= 0) then
            write(ihist,'(I5,11x,E9.2)',advance='NO') jH,bea
            write(ihist,'(7x,E9.2)') be
        end if
        info(1)  = 0
        info(2)  = 0
        rinfo(1) = DZRO
        rinfo(2) = DZRO
        irc(1)   = 0
        retlbl = 0
        return
    endif

! Compute the scaling factor for the backward error on the
!  unpreconditioned sytem

    sA       = cntl(2)
    sb       = cntl(3)
    if ((sA == DZRO) .AND. (sb == DZRO)) then
        sb = bn
    endif
! Compute the scaling factor for the backward error on the
!  preconditioned sytem

    sPA      = cntl(4)
    sPb      = cntl(5)
    if ((sPA == DZRO) .AND. (sPb == DZRO)) then
        if ((typePrec == noPrec) .OR. (typePrec == rightPrec)) then
            sPb = bn
        else
            irc(1) = precondLeft
            irc(2) = bptr
            irc(4) = r0ptr
            retlbl = 6
            return
        endif
    endif
    6 continue
    if ((sPA == DZRO) .AND. (sPb == DZRO)) then
        if ((typePrec == dblePrec) .OR. (typePrec == leftPrec)) then

        !           sPb = scnrm2(n,r0,1)

            irc(1) = prosca
            irc(2) = r0ptr
            irc(3) = r0ptr
            irc(4) = dotptr
            irc(5) = 1
            retlbl = 8
            return
        endif
    endif
    8 continue
    if ((sPA == DZRO) .AND. (sPb == DZRO)) then
        if ((typePrec == dblePrec) .OR. (typePrec == leftPrec)) then
            sPb = sqrt(real(dot(1)))

        endif
    endif


! Compute the first residual
!           Y = AX : r0 <-- A x

! The residual is computed only if the initial guess is not zero

    if (initGuess /= 0) then
        irc(1) = matvec
        irc(2) = xptr
        irc(4) = r0ptr
        retlbl = 11
        return
    endif
    11 continue
    if (initGuess /= 0) then
        do j=1,n
            r0(j) = b(j)-r0(j)
        enddo
    else
        call ccopy(n,b,1,r0,1)
    endif

! Compute the preconditioned residual if necessary
!      M_1Y = X : w <-- M_1^{-1} r0

    if ((typePrec == noPrec) .OR. (typePrec == rightPrec)) then
        call ccopy(n,r0,1,w,1)
    else
        irc(1) = precondLeft
        irc(2) = r0ptr
        irc(4) = wptr
        retlbl = 16
        return
    endif
    16 continue


!       beta = scnrm2(n,w,1)


    irc(1) = prosca
    irc(2) = wptr
    irc(3) = wptr
    irc(4) = dotptr
    irc(5) = 1
    retlbl = 18
    return
    18 continue
    beta = sqrt(real(dot(1)))

    if (beta == DZRO) then
    !  The residual is exactly zero : x is the exact solution
        info(1) = 0
        info(2) = 0
        rinfo(1) = DZRO
        rinfo(2) = DZRO
        irc(1)   = 0
        retlbl = 0
        jH = 0
        bea = DZRO
        be  = DZRO
        write(ihist,'(I5,11x,E9.2)',advance='NO') jH,bea
        write(ihist,'(7x,E9.2)') be
        if (iwarn /= 0) then
            write(iwarn,*)
            write(iwarn,*) ' WARNING GMRES : '
            write(iwarn,*) '       Intial residual is zero'
            write(iwarn,*) '       initial guess is solution'
            write(iwarn,*)
        endif
        return
    endif

    aux = ONE/beta
    do j=1,n
        V(j,1) = ZERO
    enddo
    call caxpy(n,aux,w,1,V(1,1),1)

!       Most outer loop : cgmres iteration

!       REPEAT
    7 continue


    H(1,m+1)=beta
    do j=1,m
        H(j+1,m+1) = ZERO
    enddo

!        Construction of the hessenberg matrix WORK and of the orthogonal
!        basis V such that AV=VH

    jH = 1
    10 continue
! Remark : this  do loop has been written with a while do
!          because the
!               " do jH=1,restart "
!         fails with the reverse communication.
!      do  jH=1,restart


! Compute the preconditioned residual if necessary

    if ((typePrec == rightPrec) .OR. (typePrec == dblePrec)) then

    !           Y = M_2^{-1}X : w <-- M_2^{-1} V(1,jH)

        irc(1) = precondRight
        irc(2) = vptr + (jH-1)*n
        irc(4) = wptr
        retlbl = 21
        return
    else
        call ccopy(n,V(1,jH),1,w,1)
    endif
    21 continue

!           Y = AX : r0 <-- A w

    irc(1) = matvec
    irc(2) = wptr
    irc(4) = r0ptr
    retlbl = 26
    return
    26 continue

!      MY = X : w <-- M_1^{-1} r0

    if ((typePrec == noPrec) .OR. (typePrec == rightPrec)) then
        call ccopy(n,r0,1,w,1)
    else
        irc(1) = precondLeft
        irc(2) = r0ptr
        irc(4) = wptr
        retlbl = 31
        return
    endif
    31 continue

! Orthogonalization using either MGS or IMGS

! initialize the Hessenberg matrix to zero in order to be able to use
!     IMGS as orthogonalization procedure.
    do j=1,jH
        H(j,jH) = ZERO
    enddo
    nOrtho = 0
    19 continue
    nOrtho = nOrtho +1
    dloo   = DZRO

    if ((iOrthog == 0) .OR. (iOrthog == 1)) then
    ! MGS

    !           do j=1,jH

        j = 1
    !           REPEAT
    endif
    23 continue
    if ((iOrthog == 0) .OR. (iOrthog == 1)) then

    !             dVi     = cdotc(n,V(1,j),1,w,1)

        irc(1) = prosca
        irc(2) = vptr + (j-1)*n
        irc(3) = wptr
        irc(4) = dotptr
        irc(5) = 1
        retlbl = 32
        return
    endif
    32 continue
    if ((iOrthog == 0) .OR. (iOrthog == 1)) then
        dVi     = dot(1)
        H(j,jH) = H(j,jH) + dVi
        dloo    = dloo + abs(dVi)**2
        aux = -ONE*dVi
        call caxpy(n,aux,V(1,j),1,w,1)
        j = j + 1
        if (j <= jH) goto 23
    !          enddo_j
    else
    ! CGS
    ! Gathered dot product calculation

    !           call cgemv('C',n,jH,ONE,V(1,1),n,w,1,ZERO,r0,1)

        irc(1) = prosca
        irc(2) = vptr
        irc(3) = wptr
        irc(4) = dotptr
        irc(5) = jH
        retlbl = 34
        return
    endif
    34 continue
    if ((iOrthog == 2) .OR. (iOrthog == 3)) then

        call caxpy(jH,ONE,dot,1,H(1,jH),1)
        call cgemv('N',n,jH,-ONE,V(1,1),n,dot,1,ONE,w,1)
        dloo = scnrm2(jH,dot,1)**2
    endif

!         dnormw = scnrm2(n,w,1)

    irc(1) = prosca
    irc(2) = wptr
    irc(3) = wptr
    irc(4) = dotptr
    irc(5) = 1
    retlbl = 33
    return
    33 continue
    dnormw = sqrt(real(dot(1)))

    if ((iOrthog == 1) .OR. (iOrthog == 3)) then
    ! IMGS / CGS orthogonalisation
        dloo = sqrt(dloo)
    ! check the orthogonalization quality
        if (((2.0*dnormw) <= dloo) .AND. (nOrtho < 3)) then
            goto 19
        endif
    endif

    H(jH+1,jH) = dnormw
    if ((jH < m) .OR. (icntl(8) == 0)) then
        aux = ONE/dnormw
        do j=1,n
            V(j,jH+1) = ZERO
        enddo
        call caxpy(n,aux,w,1,V(1,jH+1),1)
    endif
! Apply previous Givens rotations to the new column of H
    do j=1,jH-1
        call crot(1, H(j,jH), 1, H(j+1,jH), 1,real(rotCos(j)), &
        rotSin(j))
    enddo
    auxHjj = H(jH,jH)
    auxHjp1j= H(jH+1,jH)
    call crotg(auxHjj, auxHjp1j,temp,rotSin(jH))
    rotCos(jH)= cmplx(temp, DZRO)
! Apply current rotation to the rhs of the least squares problem
    call crot(1, H(jH,m+1), 1, H(jH+1,m+1), 1, real(rotCos(jH)), &
    rotSin(jH))

! zabs(H(jH+1,m+1)) is the residual computed using the least squares
!          solver
! Complete the QR factorisation of the Hessenberg matrix by apply the current
! rotation to the last entry of the collumn
    call crot(1, H(jH,jH), 1, H(jH+1,jH), 1, real(rotCos(jH)), &
    rotSin(jH))
    H(jH+1,jH) = ZERO

! Get the Least square residual

    dnormres = abs(H(jH+1,m+1))
    if (sPa /= DZRO) then

    ! Compute the solution of the current linear least squares problem

        call ccopy(jH,H(1,m+1),1,yCurrent,1)
        call ctrsv('U','N','N',jH,H,m+1,yCurrent,1)

    ! Compute the value of the new iterate

        call cgemv('N',n,jH,ONE,v,n, &
        yCurrent,1,ZERO,xCurrent,1)

        if ((typePrec == rightPrec) .OR. (typePrec == dblePrec)) then

        !         Y = M_2^{-1}X : r0 <-- M_2^{-1} xCurrent

            irc(1) = precondRight
            irc(2) = xcuptr
            irc(4) = r0ptr
            retlbl = 36
            return
        else
            call ccopy(n,xCurrent,1,r0,1)
        endif
    endif
    36 continue


    if (sPa /= DZRO) then
    ! Update the current solution
        call ccopy(n,x,1,xCurrent,1)
        call caxpy(n,ONE,r0,1,xCurrent,1)

    !         dnormx = scnrm2(n,xCurrent,1)

        irc(1) = prosca
        irc(2) = xcuptr
        irc(3) = xcuptr
        irc(4) = dotptr
        irc(5) = 1
        retlbl = 38
        return
    else
        dnormx    = DONE
    endif
    38 continue
    if (sPa /= DZRO) then
        dnormx = sqrt(real(dot(1)))
    endif

    bea = dnormres/(sPA*dnormx+sPb)

! Check the convergence based on the Arnoldi Backward error for the
! preconditioned system
    if ((bea <= cntl(1)) .OR. (iterOut*m+jH >= iterMax)) then

    ! The Arnoldi Backward error indicates that cgmres might have converge
    ! enforce the calculation of the true residual at next restart
        compRsd = 1

    !  If the update of X has not yet been performed
        if (sPA == DZRO) then

        ! Compute the solution of the current linear least squares problem

            call ccopy(jH,H(1,m+1),1,yCurrent,1)
            call ctrsv('U','N','N',jH,H,m+1,yCurrent,1)

        ! Compute the value of the new iterate

            call cgemv('N',n,jH,ONE,v,n, &
            yCurrent,1,ZERO,xCurrent,1)

            if ((typePrec == rightPrec) .OR. (typePrec == dblePrec)) then

            !         Y = M_2^{-1}X : r0 <-- M_2^{-1} xCurrent

                irc(1) = precondRight
                irc(2) = xcuptr
                irc(4) = r0ptr
                retlbl = 37
                return
            else
                call ccopy(n,xCurrent,1,r0,1)
            endif
        endif
    endif
    37 continue
    if ((bea <= cntl(1)) .OR. (iterOut*m+jH >= iterMax)) then
        if (sPA == DZRO) then
        ! Update the current solution
            call ccopy(n,x,1,xCurrent,1)
            call caxpy(n,ONE,r0,1,xCurrent,1)
        endif

        call ccopy(n,xCurrent,1,r0,1)
    ! Compute the true residual, the Arnoldi one may be unaccurate

    !           Y = AX : w  <-- A r0

        irc(1) = matvec
        irc(2) = r0ptr
        irc(4) = wptr
        retlbl = 41
        return
    endif
    41 continue
    if ((bea <= cntl(1)) .OR. (iterOut*m+jH >= iterMax)) then

        do j=1,n
            w(j) = b(j) - w(j)
        enddo
    ! Compute the norm of the unpreconditioned residual

    !        trueNormRes = scnrm2(n,w,1)

        irc(1) = prosca
        irc(2) = wptr
        irc(3) = wptr
        irc(4) = dotptr
        irc(5) = 1
        retlbl = 43
        return
    endif
    43 continue
    if ((bea <= cntl(1)) .OR. (iterOut*m+jH >= iterMax)) then
        trueNormRes = sqrt(real(dot(1)))

        if ((typePrec == leftPrec) .OR. (typePrec == dblePrec)) then

        !      MY = X : r0 <-- M_1^{-1} w

            irc(1) = precondLeft
            irc(2) = wptr
            irc(4) = r0ptr
            retlbl = 46
            return
        else
            call ccopy(n,w,1,r0,1)
        endif
    endif
    46 continue
    if ((bea <= cntl(1)) .OR. (iterOut*m+jH >= iterMax)) then

    !        dnormres = scnrm2(n,r0,1)

        irc(1) = prosca
        irc(2) = r0ptr
        irc(3) = r0ptr
        irc(4) = dotptr
        irc(5) = 1
        retlbl = 48
        return
    endif
    48 continue
    if ((bea <= cntl(1)) .OR. (iterOut*m+jH >= iterMax)) then
        dnormres = sqrt(real(dot(1)))

        be = dnormres/(sPA*dnormx+sPb)
    ! Save the backward error on a file if convergence history requested
        if (ihist /= 0) then
            write(ihist,1000)iterOut*m+jH,bea,be
            1000 format(I5,11x,E9.2,7x,E9.2)
        endif

    endif


! Check again the convergence
    if ((bea <= cntl(1)) .OR. (iterOut*m+jH >= iterMax)) then
        if ((be <= cntl(1)) .OR. (iterOut*m+jH >= iterMax)) then
        ! The convergence has been achieved, we restore the solution in x
        ! and compute the two backward errors.
            call ccopy(n,xCurrent,1,x,1)

            if (sA /= DZRO) then

            !            dnormx = scnrm2(n,x,1)

                irc(1) = prosca
                irc(2) = xptr
                irc(3) = xptr
                irc(4) = dotptr
                irc(5) = 1
                retlbl = 51
                return
            endif
        endif
    endif
    51 continue
    if ((bea <= cntl(1)) .OR. (iterOut*m+jH >= iterMax)) then
        if ((be <= cntl(1)) .OR. (iterOut*m+jH >= iterMax)) then
            if (sA /= DZRO) then
                dnormx = sqrt(real(dot(1)))

            else
                dnormx = DONE
            endif
        ! Return the backward errors
            rinfo(1) = be
            rinfo(2) = trueNormRes/(sA*dnormx+sb)
            if (be <= cntl(1)) then
                info(1) = 0
                if (ihist /= 0) then
                    write(ihist,*)
                    write(ihist,'(A20)') 'Convergence achieved'
                endif
            else if (be > cntl(1)) then
                if (iwarn /= 0) then
                    write(iwarn,*)
                    write(iwarn,*) ' WARNING GMRES : '
                    write(iwarn,*) '       No convergence after '
                    write(iwarn,*) iterOut*m+jH,' iterations '
                    write(iwarn,*)
                endif
                if (ihist /= 0) then
                    write(ihist,*)
                    write(ihist,*) ' WARNING GMRES :'
                    write(ihist,*) '       No convergence after '
                    write(ihist,*) iterOut*m+jH,' iterations '
                    write(ihist,*)
                endif
                info(1) = -4
            endif
            if (ihist /= 0) then
                write(ihist,1010) rinfo(1)
                write(ihist,1011) rinfo(2)
                1010 format('B.E. on the preconditioned system:   ',E9.2)
                1011 format('B.E. on the unpreconditioned system: ',E9.2)
            endif
            info(2) = iterOut*m+jH
            if (ihist /= 0) then
                write(ihist,'(A10,I2)') 'info(1) = ',info(1)
                write(ihist,'(A32,I5)') &
                'Number of iterations (info(2)): ',info(2)
            endif
            irc(1)  = 0
            retlbl  = 0
            return
        endif
    else
    ! Save the backward error on a file if convergence history requested
        if (ihist /= 0) then
            write(ihist,1001)iterOut*m+jH,bea
            1001 format(I5,11x,E9.2,7x,'--')
        endif

    endif

    jH = jH + 1
    if (jH <= m) then
        goto 10
    endif

    iterOut = iterOut + 1

! we have completed the Krylov space construction, we restart if
! we have not yet exceeded the maximum number of iterations allowed.

    if ((sPa == DZRO) .AND. (bea > cntl(1))) then

    ! Compute the solution of the current linear least squares problem

        jH = jH - 1
        call ccopy(jH,H(1,m+1),1,yCurrent,1)
        call ctrsv('U','N','N',jH,H,m+1,yCurrent,1)

    ! Compute the value of the new iterate

        call cgemv('N',n,jH,ONE,v,n, &
        yCurrent,1,ZERO,xCurrent,1)

        if ((typePrec == rightPrec) .OR. (typePrec == dblePrec)) then

        !         Y = M_2^{-1}X : r0 <-- M_2^{-1} xCurrent

            irc(1) = precondRight
            irc(2) = xcuptr
            irc(4) = r0ptr
            retlbl = 52
            return
        else
            call ccopy(n,xCurrent,1,r0,1)
        endif
    endif
    52 continue
    if ((sPa == DZRO) .AND. (bea > cntl(1))) then
    ! Update the current solution
        call ccopy(n,x,1,xCurrent,1)
        call caxpy(n,ONE,r0,1,xCurrent,1)
    endif

    call ccopy(n,xCurrent,1,x,1)

    if (compRsd == 1) then

    ! Compute the true residual

        call ccopy(n,x,1,w,1)
        irc(1) = matvec
        irc(2) = wptr
        irc(4) = r0ptr
        retlbl = 61
        return
    endif
    61 continue
    if (compRsd == 1) then
        do j=1,n
            r0(j) = b(j) - r0(j)
        enddo

    ! Precondition the new residual if necessary

        if ((typePrec == leftPrec) .OR. (typePrec == dblePrec)) then

        !      MY = X : w <-- M_1^{-1} r0

            irc(1) = precondLeft
            irc(2) = r0ptr
            irc(4) = wptr
            retlbl = 66
            return
        else
            call ccopy(n,r0,1,w,1)
        endif
    endif
    66 continue

!           beta = scnrm2(n,w,1)

    if (compRsd == 1) then
        irc(1) = prosca
        irc(2) = wptr
        irc(3) = wptr
        irc(4) = dotptr
        irc(5) = 1
        retlbl = 68
        return
    endif
    68 continue
    if (compRsd == 1) then
        beta = sqrt(real(dot(1)))

    else
    ! Use recurrence to approximate the residual at restart
        beta = abs(H(m+1,m+1))
    ! Apply the Givens rotation is the reverse order
        do j=m,1,-1
            H(j,m+1)   = ZERO
            call crot(1, H(j,m+1), 1, H(j+1,m+1), 1, &
            real(rotCos(j)), -rotSin(j))
        enddo

    ! On applique les vecteurs V

        call cgemv('N',n,m+1,ONE,v,n,H(1,m+1),1,ZERO,w,1)

    endif
    do j=1,n
        V(j,1) = ZERO
    enddo
    aux = ONE/beta
    call caxpy(n,aux,w,1,V(1,1),1)

    goto 7

    end subroutine cgmres


    subroutine init_cgmres(icntl,cntl, env)

!  Purpose
!  =======
!    Set default values for the parameters defining the characteristics
! of the Gmres algorithm.
!  See the User's Guide for an example of use.


! Written : April 1997
! Authors :  Valerie Fraysse, Luc Giraud, Serge Gratton
!             Parallel Algorithms - CERFACS


!  Arguments
!  =========

! icntl    (input) INTEGER array. length 6
!            icntl(1) : stdout for error messages
!            icntl(2) : stdout for warnings
!            icntl(3) : stdout for convergence history
!            icntl(4) : 0 - no preconditioning
!                       1 - left preconditioning
!                       2 - right preconditioning
!                       3 - double side preconditioning
!                       4 - error, default set in Init
!            icntl(5) : 0 - modified Gram-Schmidt
!                       1 - iterative modified Gram-Schmidt
!                       2 - classical Gram-Schmidt
!                       3 - iterative classical Gram-Schmidt
!            icntl(6) : 0 - default initial guess x_0 = 0 (to be set)
!                       1 - user supplied initial guess
!            icntl(7) : maximum number of iterations
!            icntl(8) : 1 - default compute the true residual at each restart
!                       0 - use recurence formaula at restart

! cntl     (input) real array, length 5
!            cntl(1) : tolerance for convergence
!            cntl(2) : scaling factor for normwise perturbation on A
!            cntl(3) : scaling factor for normwise perturbation on b
!            cntl(4) : scaling factor for normwise perturbation on the
!                      preconditioned matrix
!            cntl(5) : scaling factor for normwise perturbation on
!                      preconditioned right hand side

! Output variables
! ----------------
    integer :: icntl(*)
    real ::   cntl(*)
    TYPE(GMRESVAR), TARGET :: env

    icntl(1) = 6
    icntl(2) = 6
    icntl(3) = 0
    icntl(4) = 4
    icntl(5) = 0
    icntl(6) = 0
    icntl(7) = -1
    icntl(8) = 1

    cntl(1) = 1.0e-5
    cntl(2) = 0.0e0
    cntl(3) = 0.0e0
    cntl(4) = 0.0e0
    cntl(5) = 0.0e0

    env%icheck = 0
    env%retlbl = 0

    return
    end subroutine init_cgmres


    END MODULE
