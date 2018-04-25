!=====================================================================!
! Vector datatype module that computes norm in a distributed fashion
! 
! Author : Komahan Boopathy (komahan@gaech.edu)
!=====================================================================!

module vector_class

  use iso_fortran_env, only : dp => REAL64
  use opencoarrays, only : co_sum
  
  implicit none

  ! Everything is private to module
  private

  ! Only the datatype is exposed to public
  public :: vector

  ! Vector datatype
  type :: vector

     ! Attributes
     real(dp), allocatable :: values(:)
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
  pure function reduce_sum(a, b)
    
    real(dp), intent(in) :: a, b
    real(dp) :: reduce_sum
    reduce_sum = a + b
  end function reduce_sum

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
    real(dp)       , intent(in)    :: alpha
    type(vector)  , intent(in)    :: x

    this % values = this % values + alpha * x % values

  end subroutine axpy
  
  !===================================================================!
  ! y <- beta*y + alpha*x where alpha is a scalar
  !===================================================================!

  subroutine axpby(this, alpha, x, beta)

    class(vector) , intent(inout) :: this
    real(dp)       , intent(in)    :: alpha, beta
    type(vector)  , intent(in)    :: x

    this % values = beta * this % values + alpha * x % values

  end subroutine axpby

  !===================================================================!
  ! Scale the values by alpha: y <- alpha*y
  !===================================================================!

  subroutine scale(this, alpha)

    class(vector) , intent(inout) :: this
    real(dp)       , intent(in)    :: alpha

    this % values = alpha * this % values

  end subroutine scale
  
  !===================================================================!
  ! Return the dot product of this vector with another vector 'b'
  !===================================================================!

  real(dp) function dot(this, b) result(xdot)

    class(vector), intent(in) :: this
    type(vector), intent(in) :: b

    xdot = dot_product(this % values, b % values)
    call co_sum(xdot)

  end function dot

  !===================================================================!
  ! Norm of the vector
  !===================================================================!

  real(dp) function norm(this)

    class(vector), intent(in) :: this
    real(dp) :: xdot

    xdot = dot_product(this % values, this % values)
    call co_sum (xdot)
    norm = sqrt(xdot)

  end function norm

end module vector_class
