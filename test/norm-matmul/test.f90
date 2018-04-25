program main

  use clock_class, only : clock
  type(clock) :: timer

  call timer % start()
  call test_norm()
  call timer % stop()

  if (this_image() .eq. 1) then
     write(*, '("Model run time:",F8.3," seconds")') timer % getelapsed()
  end if

end program main

!=====================================================================!
! Program to compute the norm of a large vector in a distributed
! fashion using fortran coarrays and collective routines.
!
! Author : Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

subroutine test_norm
  
  use vector_class, only : vector

  implicit none

  integer, parameter :: global_size = 2000000
  integer :: nimages = 1
  integer :: local_size

  ! Distributed work data
  real(8), allocatable :: x(:)
  type(vector) :: xvec
  real(8) :: xdot
  real(8) :: xnorm
    
  ! Determine partition
  nimages = num_images()
  local_size = global_size/nimages

  ! Create smaller vectors
  allocate(x(local_size))

  !-------------------------------------------------------------------!
  ! Test norm computations of large vector with random values
  !-------------------------------------------------------------------!

  test_random: block

    ! Set random values into the vector
    call random_number(x)

    ! Using direct procedure
    xdot = dot_product(x,x)  
    call co_sum (xdot)
    if (this_image() == 1) then
       write(*,*) "Norm of the array using direct procedure is", sqrt(xdot) , this_image()
    end if

    ! Using CO_NORM2 function
    xnorm = co_norm2(x)
    if (this_image() == 1) then
       write(*,*) "Norm of the vector using co_norm2 function is", xnorm, this_image()
    end if

    ! Using derived datatype that uses direct procedure internally
    xvec % values = x
    xnorm = xvec % norm()
    if (this_image() == 1) then
       write(*,*) "Norm of the vector datatype is", xnorm , this_image()
    end if
    
  end block test_random

  !-------------------------------------------------------------------!
  ! Test norm computations of large vector with same values
  !-------------------------------------------------------------------!

  test_deterministic: block
    
    ! Set deterministic values into the vector
    x = 1.0d0

    ! Using direct procedure
    xdot = dot_product(x,x)  
    call co_sum (xdot)
    if (this_image() == 1) then
       write(*,*) "Norm of the array using direct procedure is", sqrt(xdot) , this_image()
    end if

    ! Using CO_NORM2 function
    xnorm = co_norm2(x)
    if (this_image() == 1) then
       write(*,*) "Norm of the vector using co_norm2 function is", xnorm, this_image()
    end if
    
    ! Using derived datatype that uses direct procedure internally
    xvec % values = x
    xnorm = xvec % norm()
    if (this_image() == 1) then
       write(*,*) "Norm of the vector datatype is", xnorm , this_image()
    end if

  end block test_deterministic

  deallocate(x)

contains
  
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

  ! MPI_SUM
  pure function sum(a, b)
    real(8), intent(in) :: a, b
    real(8) :: sum
    sum = a + b
  end function sum

end subroutine test_norm
