program gmresvsgauss
    USE GMRES
    USE M_SOLVER
    USE IFPORT
    implicit none


    COMPLEX, DIMENSION(:,:), ALLOCATABLE :: random_matrix, gauss_matrix, gmres_matrix
    COMPLEX, DIMENSION(:), ALLOCATABLE :: random_sol, gmres_sol

    integer:: arg_count, num_iterations,i, n,j,iter
    character(200) :: arg
    real:: cpu_time_begin, cpu_time_end, total_cpu_gauss, total_cpu_gmres
    DOUBLE PRECISION :: WALL_START_TIME, WALL_STOP_TIME, total_wall_gauss, total_wall_gmres
    complex :: ONE, ZERO

    n = 2000
    ONE = CMPLX(1., 0.)
    ZERO = CMPLX(0., 0.)

    arg_count = command_argument_count()

    num_iterations = 1

    if(arg_count >=1 ) then
        call get_command_argument(1, arg)
        read(arg,*) num_iterations

    end if

    !CALL RANDOM_SEED()
    total_cpu_gmres = 0
    total_cpu_gauss = 0
    total_wall_gmres = 0
    total_wall_gauss = 0

    do iter=1,num_iterations

        write(*,*) 'Iteration ', iter
        write(*,*) 'Generates ', n , ' x ', n, ' matrix'
        ALLOCATE(random_matrix(n,n))
        ALLOCATE(random_sol(n))
        ALLOCATE(gauss_matrix(n,n+1), gmres_matrix(n,n), gmres_sol(n))


        random_matrix = cmplx(0.,0.)
        do i=1,n
            random_matrix(i,i) = (4.0e0,  0.0e0)
        end do
        do i = 1,n-1
            random_matrix(i,i+1) = (-2.0e0, 1.0e0)
            random_matrix(i+1,i) = (-1.0e0, 1.0e0)
        enddo

        gmres_sol = ONE

        call cgemv('N',n,n,ONE,random_matrix,n,gmres_sol,1,ZERO,random_sol,1)

        do i=1,n
            do j=1,n
                gauss_matrix(i, j) = random_matrix(i,j)
                gmres_matrix(i, j) = random_matrix(i,j)
            end do
            gauss_matrix(i, n+1) = random_sol(i)
            gmres_sol(i) = random_sol(i)
        end do

        write(*,*) 'Solving with gmres started'

        CALL CPU_TIME ( cpu_time_begin )
        WALL_START_TIME = DCLOCK()

        CALL compute_gmres(gmres_matrix, gmres_sol, n)

        WALL_STOP_TIME = DCLOCK()
        CALL CPU_TIME ( cpu_time_end )

        write(*, *) 'GMRES completed with following statistics:'
        write(*,*) 'cpu time(s) ' , cpu_time_end - cpu_time_begin
        write(*,*) 'wall elapsed time(s) ', WALL_STOP_TIME-WALL_START_TIME
        write(*,*)
        total_cpu_gmres = total_cpu_gmres + cpu_time_end - cpu_time_begin
        total_wall_gmres = total_wall_gmres + WALL_STOP_TIME-WALL_START_TIME

        write(*,*) 'Solving with gauss started'
        CALL CPU_TIME ( cpu_time_begin )
        WALL_START_TIME = DCLOCK()

        CALL compute_gauss(gauss_matrix, n)

        WALL_STOP_TIME = DCLOCK()
        CALL CPU_TIME ( cpu_time_end )

        write(*, *) 'GAUSS completed with following statistics:'
        write(*,*) 'cpu time(s) ' , cpu_time_end - cpu_time_begin
        write(*,*) 'wall elapsed time(s) ', WALL_STOP_TIME-WALL_START_TIME
        write(*,*)
        total_cpu_gauss = total_cpu_gauss + cpu_time_end - cpu_time_begin
        total_wall_gauss = total_wall_gauss + WALL_STOP_TIME-WALL_START_TIME

        DEALLOCATE(random_matrix)
        DEALLOCATE(random_sol)
        DEALLOCATE(gauss_matrix, gmres_matrix, gmres_sol)


    end do

    if (num_iterations >= 1) then

        total_cpu_gauss = total_cpu_gauss/real(num_iterations)
        total_cpu_gmres = total_cpu_gmres/real(num_iterations)
        total_wall_gauss = total_wall_gauss/ dble(num_iterations)
        total_wall_gmres = total_wall_gmres/ dble(num_iterations)

        write(*,*) '----------  Average Summary -------------'
        write(*, *) 'GMRES completed with following statistics:'
        write(*,*) 'avg cpu time(s) ' , total_cpu_gmres
        write(*,*) 'avg wall elapsed time(s) ', total_wall_gmres
        write(*,*)
        write(*, *) 'GAUSS completed with following statistics:'
        write(*,*) 'avg cpu time(s) ' , total_cpu_gauss
        write(*,*) 'avg wall elapsed time(s) ', total_wall_gauss
        write(*,*)


    end if




contains
    subroutine compute_gmres(A, b, N)

        COMPLEX,DIMENSION(:,:):: A
        COMPLEX, DIMENSION(:):: b
        INTEGER::N
        integer :: lwork
        integer :: revcom, colx, coly, colz, nbscal
        integer :: irc(5), icntl(8), info(3)
        TYPE(GMRESVAR):: gmres_env
        integer :: matvec, precondLeft, precondRight, dotProd
        parameter (matvec=1, precondLeft=2, precondRight=3, dotProd=4)
        complex, dimension(:), allocatable ::  work
        real ::  cntl(5), rinfo(2)
        complex :: ONE, ZZERO
        INTEGER:: restartParam, maxIterations
        real:: TOL_GMRES

        restartParam = 20
        maxIterations = 100
        TOL_GMRES = 5.E-07

        lwork = (restartParam+10)**2 + (restartParam+10)*(N+5) + 5*N + 1
        ALLOCATE(work(lwork))
        ONE = CMPLX(1., 0.)
        ZZERO = CMPLX(0., 0.)


        !******************************************************
        !* Initialize the control parameters to default value
        !******************************************************

        call init_cgmres(icntl,cntl, gmres_env)

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

        work(1:N) = ONE/2.0

        DO I=1,N

            work(N+I) = b(i)

        END DO

        DO

            call drive_cgmres(N,N,restartParam,lwork,work, &
                irc,icntl,cntl,info,rinfo, gmres_env)
            revcom = irc(1)
            colx   = irc(2)
            coly   = irc(3)
            colz   = irc(4)
            nbscal = irc(5)

            if (revcom == matvec) then
                ! perform the matrix vector product
                !        work(colz) <-- A * work(colx)
                call cgemv('N',N,N,ONE,A,N,work(colx),1, &
                    ZZERO,work(colz),1)

            else if (revcom == precondLeft) then
                ! perform the left preconditioning
                !         work(colz) <-- M^{-1} * work(colx)
                call ccopy(N,work(colx),1,work(colz),1)
                call ctrsm('L','L','N','N',N,1,ONE,A,N,work(colz),N)

            else if (revcom == precondRight) then
                ! perform the right preconditioning
                call ccopy(N,work(colx),1,work(colz),1)
                call ctrsm('L','U','N','N',N,1,ONE,A,N,work(colz),N)

            else if (revcom == dotProd) then
                !      perform the scalar product
                !      work(colz) <-- work(colx) work(coly)

                call cgemv('C',N,nbscal,ONE,work(colx),N, &
                    work(coly),1,ZZERO,work(colz),1)

            ELSE
                exit
            endif

        END DO

    end subroutine compute_gmres

    subroutine compute_gauss(A, N)
        integer:: n
        COMPLEX,DIMENSION(:,:):: A
        CALL GAUSSZ(A,N,N,N)
    end subroutine compute_gauss





end program gmresvsgauss
