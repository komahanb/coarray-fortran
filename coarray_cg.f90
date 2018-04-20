module conjugate_gradient

  ! import dependencies
  use iso_fortran_env, only : dp => REAL64

  implicit none

contains

  subroutine dparcg(A, b, max_it, max_tol, x, iter, tol, flag)

    real(dp), intent(in) :: A(:,:)
    real(dp), intent(in) :: b(:)
    integer , intent(in) :: max_it
    real(dp), intent(in) :: max_tol

    real(dp), intent(inout) :: x(:)
    integer , intent(out)   :: iter
    real(dp), intent(out)   :: tol
    integer , intent(out)   :: flag

    ! create local data
    real(dp), allocatable :: p(:), r(:), w(:)
    real(dp), allocatable :: rho(:)
    real(dp) :: alpha, beta
    real(dp) :: bnorm, rnorm

    ! Memory allocations
    allocate(r, p, w, mold=x)
    allocate(rho(max_it))

    ! Start the iteration counter
    iter = 1

    ! Norm of the right hand side
    bnorm = co_norm2(b)

    ! Norm of the initial residual
    r         = b - co_matmul(A, x)
    rnorm     = co_norm2(r)
    tol       = rnorm/bnorm
    rho(iter) = rnorm*rnorm

    open(10, file='cg.log', action='write', position='append')

    ! Apply Iterative scheme until tolerance is achieved
    do while ((tol .gt. max_tol) .and. (iter .lt. max_it))

       ! step (a) compute the descent direction
       if ( iter .eq. 1) then
          ! steepest descent direction p
          p = r
       else
          ! take a conjugate direction
          beta = rho(iter)/rho(iter-1)
          p = r + beta*p
       end if

       ! step (b) compute the solution update
       w = co_matmul(A,p)

       ! step (c) compute the step size for update
       alpha = rho(iter)/dot_product(p, w)

       ! step (d) Add dx to the old solution
       x = x + alpha*p

       ! step (e) compute the new residual
       r = r - alpha*w
       !r = b - matmul(A, x)

       ! step(f) update values before next iteration
       rnorm = co_norm2(r)
       tol = rnorm/bnorm

       write(10,*) iter, tol
       print *, iter, tol

       iter = iter + 1

       rho(iter) = rnorm*rnorm

    end do

    close(10)

    deallocate(r, p, w, rho)

    flag = 0

  end subroutine dparcg

  !===================================================================!
  ! Function to compute the norm of a distributed vector
  !===================================================================!

  function co_norm2(x) result(norm)

    real(8), intent(in) :: x(:)    
    real(8) :: xdot, norm

    ! find dot product, sum over processors, take sqrt and return
    xdot = dot_product(x,x)  
    call co_sum (xdot)
    norm = sqrt(xdot)

  end function co_norm2

  !===================================================================!
  ! Function that computes the matrix vector product in a distributed
  ! fashion for columnwise decomposition of matrix
  ! ===================================================================!

  function co_matmul(A, x) result(b)

    ! Arguments
    real(8), intent(in) :: A(:,:)
    real(8), intent(in) :: x(:)
    real(8) :: b(size(x))

    ! Local variables
    integer :: nimages
    integer :: me, local_size
    integer :: stat
    character(10) :: msg

    ! Create a local vector of global sizse (optimize this!)
    real(8), allocatable :: work(:)
    allocate(work, mold=A(:,1))

    ! Determine partition
    nimages = num_images()
    me = this_image()
    local_size = size(x)

    ! Multiply, sum and distrbute
    work = matmul(A,x)
    call co_sum(work, stat=stat, errmsg=msg)
    b = work((me-1)*local_size+1:me*local_size)

    deallocate(work)

  end function co_matmul

end module conjugate_gradient

module system

  implicit none

contains

  !-------------------------------------------------------------------!
  ! Assemble -U_xx = 2x - 0.5, U(0) = 1;  U(1)= 0, x in [0,1]
  !-------------------------------------------------------------------!

  subroutine assemble_system_dirichlet(a, b, npts, V, rhs, u, P)

    implicit none

    real(8), intent(in)  :: a, b ! bounds of the domain
    integer, intent(in)  :: npts ! number of interior points
    real(8), intent(out) :: V(npts,npts) ! banded matrix
    real(8), intent(out) :: rhs(npts)
    real(8), intent(out) :: u(npts)
    real(8), intent(out) :: P(npts,npts)
    real(8) :: S(npts,npts), D(npts, npts)

    real(8), parameter :: PI = 3.141592653589793d0
    real(8)            :: h, alpha
    integer            :: M, N
    integer            :: i, j, k

    ! h = width / num_interior_pts + 1
    h = (b-a)/dble(npts+1)
    V = 0.0d0

    ! Size of the linear system = unknowns (interior nodes)
    M = npts ! nrows
    N = npts ! ncols

    ! Set the inner block
    rows: do i = 1, M
       cols: do j = 1, N
          if (j .eq. i-1) then
             ! lower triangle
             V(j,i) = -1.0d0
          else if (j .eq. i+1) then           
             ! upper triangle
             V(j,i) = -1.0d0
          else if (j .eq. i) then           
             ! diagonal
             V(j,i) = 2.0d0
          else
             ! skip
          end if
       end do cols
    end do rows

    ! Assemble the RHS
    do i = 1, M
       rhs(i) = h*h*(2.0d0*dble(i)*h - 0.5d0)
    end do
    rhs(1) = rhs(1) + 1.0d0
    rhs(M) = rhs(M)

    ! Initial solution profile use sin function as a first guess
    do i = 1, M
       u(i) =  sin(dble(i)*h*PI)
    end do

    ! Find the sine transform matrix
    alpha = sqrt(2.0d0/dble(npts+1))
    do j = 1, M
       do k = 1, N
          S(k,j) = alpha*sin(PI*dble(j*k)/dble(npts+1))
       end do
    end do

    ! Find the diagonal matrix
    D = matmul(S, matmul(V, S))

    ! Invert the digonal matrix easily
    do j = 1, M
       D(j,j) = 1.0d0/D(j,j)
    end do

    ! Define the preconditioner
    p = matmul(S, matmul(D, S))

  end subroutine assemble_system_dirichlet

end module system

program main

  use conjugate_gradient, only: dparcg
  use system, only : assemble_system_dirichlet

!!$  serial : block
!!$    
!!$    integer, parameter :: npts = 64
!!$    real(8), parameter :: max_tol = 1.0d-8
!!$    integer, parameter :: max_it = 100000
!!$    real(8) :: x(npts,3), b(npts), A(npts,npts), P(npts, npts)
!!$    integer :: iter, flag, i, j
!!$    real(8) :: tol
!!$
!!$    ! solve using CG
!!$    call assemble_system_dirichlet(0.0d0, 1.0d0, npts, A, b, x(:,2), P) 
!!$    call dparcg(A, b, max_it, max_tol, x(:,2), iter, tol, flag)
!!$    print *, 'cg', tol, iter
!!$
!!$  end block serial

  parallel : block

    integer, parameter :: npts = 8
    integer, parameter :: max_it = 100000
    real(8), parameter :: max_tol = 1.0d-8

    real(8), allocatable :: x(:)[:], b(:)[:], A(:,:)[:], P(:,:)[:]
    real(8), allocatable :: xtmp(:), btmp(:), Atmp(:,:)

    integer :: iter, flag, i, j
    real(8) :: tol   
    integer :: nimages
    integer :: me, local_size

    allocate(A(npts,npts)[*])
    allocate(P(npts,npts)[*])
    allocate(x(npts)[*])
    allocate(b(npts)[*])

    ! Determine partition
    nimages    = num_images()
    me         = this_image()
    local_size = npts/nimages

    ! Assemble system on master
    if (me .eq. 1) then
       call assemble_system_dirichlet(0.0d0, 1.0d0, npts, A, b, x, P)
       print *, A, shape(A)
    end if

    sync all

    ! Split A, b, x into pieces
    allocate(Atmp(npts,local_size))
    allocate(xtmp(local_size))
    allocate(btmp(local_size))

    ! Copy from proc 1
    Atmp = A(1:npts,(me-1)*local_size+1:me*local_size)[1]
    xtmp = x((me-1)*local_size+1:me*local_size)[1]
    btmp = b((me-1)*local_size+1:me*local_size)[1]

    ! clearup memory
    deallocate(A,b,x,P)

    ! Distribute the work to processors
    call dparcg(Atmp, btmp, max_it, max_tol, xtmp, iter, tol, flag)
    print *, 'cg', tol, iter

  end block parallel

end program main
