program main

  use clock_class, only : clock

  type(clock) :: timer

  call timer % start()
  call test_norm()
  call timer % stop()

  if (this_image() .eq. 1) then
     write(*, '("model run time:",f8.3," seconds")') timer % getelapsed()
  end if

end program main

!=====================================================================!
! program to compute the norm of a large vector in a distributed
! fashion using fortran coarrays and collective routines.
!
! author : komahan boopathy (komahan@gatech.edu)
!=====================================================================!

subroutine test_norm

  use vector_class, only : vector
  use coarray_util, only: co_norm2

  implicit none

  integer, parameter :: global_size = 2e7
  integer :: nimages = 1
  integer :: local_size

  ! distributed work data
  real(8), allocatable :: x(:)
  type(vector) :: xvec
  real(8) :: xdot
  real(8) :: xnorm

  ! determine partition
  nimages = num_images()
  local_size = global_size/nimages

  ! create smaller vectors
  allocate(x(local_size))

  !-------------------------------------------------------------------!
  ! test norm computations of a large vector with random values on each rank
  !-------------------------------------------------------------------!

  print *, nimages, this_image()

  test_random: block

    ! set random values into the vector
    call random_number(x)

    ! using co_norm2 function
    xnorm = co_norm2(x)
    write(*,*) "norm of the vector using co_norm2 function is", xnorm, this_image()

    ! using derived datatype that uses direct procedure internally
    xvec % values = x
    xnorm = xvec % norm()
    write(*,*) "norm of the vector datatype is", xnorm , this_image()

  end block test_random

  !-------------------------------------------------------------------!
  ! test norm computations of a large vector with same values on each rank
  !-------------------------------------------------------------------!

  test_deterministic: block

    ! set deterministic values into the vector
    x = 1.0d0

    ! using co_norm2 function
    xnorm = co_norm2(x)
    write(*,*) "norm of the vector using co_norm2 function is", xnorm, this_image()

    ! using derived datatype that uses direct procedure internally
    xvec % values = x
    xnorm = xvec % norm()

    write(*,*) "norm of the vector datatype is", xnorm , this_image()

  end block test_deterministic

  deallocate(x)

contains

end subroutine test_norm
