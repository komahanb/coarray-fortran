!=====================================================================!
! Program to compute the norm of a large vector in a distributed
! fashion using fortran coarrays and collective routines.

! Author : Komahan Boopathy (komahan@gatech.edu)
!
!=====================================================================!

program test_norm

  integer, parameter :: global_size = 1000000000
  integer :: nimages = 1
  integer :: local_size

  ! Distributed work data
  real(8), allocatable :: x(:)[:]
  real(8) :: xdot

  ! Determine partition
  nimages = num_images()
  local_size = global_size/nimages

  ! Create smaller vectors
  allocate(x(local_size)[*])
  
  ! Fill the local arrays with random values or known entries
  ! CALL RANDOM_NUMBER(x)
  x = 1.0d0

  ! Compute norm = sqrt(x.x) and reduce across all processors
  xdot = dot_product(x,x)  
  call co_reduce (xdot, result_image=1, operator=sum)
  if (this_image() == 1) then
     write(*,*) "Norm of the vector is", sqrt(xdot) , this_image()
  end if

contains

  ! MPI_SUM
  pure function sum(a, b)
    real(8), intent(in) :: a, b
    real(8) :: sum
    sum = a + b
  end function sum
  
end program test_norm
