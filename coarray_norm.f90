!=====================================================================!
! Vector datatype module that computes norm in a distributed fashion
! 
! Author : Komahan Boopathy (komahan@gaech.edu)
!=====================================================================!

module vector_class

  implicit none

  ! Everything is private to module
  private

  ! Only the datatype is exposed to public
  public :: vector

  ! Vector datatype
  type :: vector

     ! Attributes
     real(8), allocatable :: values(:)
     integer              :: size

   contains

     ! Type-bound procedures
     procedure :: zero_values, copy_values, random_values
     procedure :: scale, axpy, axpby
     procedure :: dot, norm

  end type vector

  ! Constructor interface
  interface vector
     module procedure construct_vector
  end interface vector

contains

  ! MPI_SUM
  pure function sum(a, b)
    
    real(8), intent(in) :: a, b
    real(8) :: sum
    sum = a + b
  end function sum

  !===================================================================!
  ! Constructor for vector
  !===================================================================!

  pure type(vector) function construct_vector(size) result (this)

    integer, intent(in) :: size

    this % size = size

    allocate(this % values(this % size))

    call this % zero_values()

  end function construct_vector

  !===================================================================!
  ! Set the values of the vector as zero
  !===================================================================!

  pure subroutine zero_values(this)

    class(vector), intent(inout) :: this

    this % values = 0.0d0

  end subroutine zero_values
  
  !===================================================================!
  ! Copies the values of one vector into another
  !===================================================================!

  subroutine copy_values(this, a)

    class(vector), intent(inout) :: this
    type(vector), intent(in) :: a

    if (a % size .ne. this % size) &
         & stop "dimension-mismatch-for-copy"
    
    this % values = a % values

  end subroutine copy_values

  !===================================================================!
  ! Set random values into the vector
  !===================================================================!

  subroutine random_values(this)

    class(vector), intent(inout) :: this

    call random_number(this % values)

  end subroutine random_values

  !===================================================================!
  ! y <- y + alpha*x where alpha is a scalar
  !===================================================================!

  subroutine axpy(this, alpha, x)

    class(vector) , intent(inout) :: this
    real(8)       , intent(in)    :: alpha
    type(vector)  , intent(in)    :: x

    this % values = this % values + alpha * x % values

  end subroutine axpy
  
  !===================================================================!
  ! y <- beta*y + alpha*x where alpha is a scalar
  !===================================================================!

  subroutine axpby(this, alpha, x, beta)

    class(vector) , intent(inout) :: this
    real(8)       , intent(in)    :: alpha, beta
    type(vector)  , intent(in)    :: x

    this % values = beta * this % values + alpha * x % values

  end subroutine axpby

  !===================================================================!
  ! Scale the values by alpha: y <- alpha*y
  !===================================================================!

  subroutine scale(this, alpha)

    class(vector) , intent(inout) :: this
    real(8)       , intent(in)    :: alpha

    this % values = alpha * this % values

  end subroutine scale
  
  !===================================================================!
  ! Return the dot product of this vector with another vector 'b'
  !===================================================================!

  real(8) function dot(this, b) result(xdot)

    class(vector), intent(in) :: this
    type(vector), intent(in) :: b

    xdot = dot_product(this % values, b % values)
    !call co_reduce (xdot, operator=sum)
    call co_sum(xdot)

  end function dot

  !===================================================================!
  ! Norm of the vector
  !===================================================================!

  real(8) function norm(this)

    class(vector), intent(in) :: this
    real(8) :: xdot

    xdot = dot_product(this % values, this % values)
    !call co_reduce (xdot, operator=sum)
    call co_sum (xdot)
    norm = sqrt(xdot)

  end function norm

end module vector_class

!=====================================================================!
! Program to compute the norm of a large vector in a distributed
! fashion using fortran coarrays and collective routines.
!
! Author : Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

subroutine test_norm
  
  use vector_class, only : vector

  implicit none

  integer, parameter :: global_size = 2000000000
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

!!$    ! Set random values into the vector
!!$    call random_number(x)
!!$
!!$    ! Using direct procedure
!!$    xdot = dot_product(x,x)  
!!$    call co_sum (xdot)
!!$    if (this_image() == 1) then
!!$       write(*,*) "Norm of the array using direct procedure is", sqrt(xdot) , this_image()
!!$    end if

!!$    ! Using CO_NORM2 function
!!$    xnorm = co_norm2(x)
!!$    if (this_image() == 1) then
!!$       write(*,*) "Norm of the vector using co_norm2 function is", xnorm, this_image()
!!$    end if
    
!!$    ! Using derived datatype that uses direct procedure internally
!!$    xvec % values = x
!!$    xnorm = xvec % norm()
!!$    if (this_image() == 1) then
!!$       write(*,*) "Norm of the vector datatype is", xnorm , this_image()
!!$    end if
    
  end block test_random

  !-------------------------------------------------------------------!
  ! Test norm computations of large vector with same values
  !-------------------------------------------------------------------!

  test_deterministic: block
    
    ! Set deterministic values into the vector
    x = 1.0d0

    ! Using direct procedure
!!$    xdot = dot_product(x,x)  
!!$    call co_sum (xdot)
!!$    if (this_image() == 1) then
!!$       write(*,*) "Norm of the array using direct procedure is", sqrt(xdot) , this_image()
!!$    end if

    ! Using CO_NORM2 function
    xnorm = co_norm2(x)
    if (this_image() == 1) then
       write(*,*) "Norm of the vector using co_norm2 function is", xnorm, this_image()
    end if
    
    ! Using derived datatype that uses direct procedure internally
!!$    xvec % values = x
!!$    xnorm = xvec % norm()
!!$    if (this_image() == 1) then
!!$       write(*,*) "Norm of the vector datatype is", xnorm , this_image()
!!$    end if

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
