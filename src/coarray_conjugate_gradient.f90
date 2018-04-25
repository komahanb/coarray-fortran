!======================================================================!
! 
!======================================================================!
module physics_class

  implicit none

  private
  public :: physics

  type :: physics

     integer :: nvars

   contains
     
     !procedure :: add_residual
     !procedure :: add_jacobian_vector_pdt

  end type physics

contains

!!$  subroutine add_residual(this, res)
!!$
!!$    res % values =
!!$  end subroutine add_residual
  
 end module physics_class

!======================================================================!
! Module that contains parallel nonlinear solvers
!
! Author : Komahan Boopathy (komahan@gatech.edu)
!======================================================================!

module parallel_nonlinear_solvers

  implicit none
 
  private
  
  public :: parcg

  type :: parcg
     
     integer :: max_iterations
     real(8) :: tolerance
     integer :: error_flag

   contains
     
     procedure :: solve
     
  end type parcg
  
  interface nonlinear_solver
     module procedure parcg
     !module procedure pargmres
  end interface nonlinear_solver
  
contains

  !-------------------------------------------------------------------!
  ! Solve a linear system Ax = b using matrix-free conjugate gradient
  ! method.
  !-------------------------------------------------------------------!

  subroutine solve(this, system, x)

    use vector_class, only : vector
    use physics_class, only : physics

    ! input and output variables
    class(parcg) , intent(inout) :: this
    type(physics), intent(in)    :: system
    type(vector), intent(inout) :: x

    ! temporary variables
    type(vector)         :: r, p, w, b
    real(8), allocatable :: rho(:)
    real(8)              :: alpha, beta
    real(8)              :: bnorm, rnorm
    real(8)              :: tolerance
    integer              :: iteration

    ! create local data
    ! real(8), allocatable :: p(:), r(:), w(:)
    ! real(8), allocatable :: rho(:)

    this % error_flag = 0
    !this % error_msg = ''

    ! Memory allocations
    allocate(rho(this % max_iterations))
    r = vector(system % nvars/num_images())
    p = vector(system % nvars/num_images())
    w = vector(system % nvars/num_images())
    b = vector(system % nvars/num_images())

    ! Start the iteration counter
    iteration = 1

    ! Norm of the right hand side
!!    call system % assemble_rhs(b)
    bnorm = b % norm()

    ! Norm of the initial residual
    !r  = b - matmul(A, x)
!!    call system % add_jacobian_vector_pdt(x, r)
    call r % axpby(1.0d0, b, -1.0d0)
    !rnorm = norm2(r)
    rnorm = r % norm()

    tolerance = rnorm/bnorm
    rho(iteration) = rnorm*rnorm

    open(10, file='cg.log', action='write', position='append')

    ! Apply Iterative scheme until tolerance is achieved
    do while ( &
         & (tolerance .gt. this % tolerance) &
         & .and. &
         & (iteration .lt. this % max_iterations) &
         & )

       ! step (a) compute the descent direction
       if ( iteration .eq. 1) then
          ! steepest descent direction p
          call p % copy_values(r)
       else
          ! take a conjugate direction
          beta = rho(iteration)/rho(iteration-1)
          !p = r + beta*p
          call p % axpby(1.0d0, r, beta)
       end if

       ! step (b) compute the solution update
       !w = matmul(A,p)
!!       call system % add_jacobian_vector_pdt(p, w)

       ! step (c) compute the step size for update
       alpha = rho(iteration)/(p % dot(w))

       ! step (d) Add dx to the old solution
       !x = x + alpha*p
       call x % axpy(alpha, p)

       ! step (e) compute the new residual
       !r = r - alpha*w
       call r % axpy(-alpha, w)
       !r = b - matmul(A, x)

       ! step(f) update values before next iteration
       !rnorm = norm2(r)
       rnorm = r % norm()
       tolerance = rnorm/bnorm

       write(10,*) iteration, tolerance
       print *, iteration, tolerance

       iteration = iteration + 1

       rho(iteration) = rnorm*rnorm

    end do

    close(10)

    if (tolerance .gt. this % tolerance) this % error_flag = 1

    deallocate(rho)

  end subroutine solve

end module parallel_nonlinear_solvers

program test_nonlinear_solvers

  use parallel_nonlinear_solvers, only : parcg

  implicit none

  !type(system) :: laplace

end program test_nonlinear_solvers
