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
    call co_reduce (xdot, operator=sum)

  end function dot

  !===================================================================!
  ! Norm of the vector
  !===================================================================!

  real(8) function norm(this)

    class(vector), intent(in) :: this
    real(8) :: xdot

    xdot = dot_product(this % values, this % values)
    call co_reduce (xdot, operator=sum)
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

  implicit none

  integer, parameter :: global_size = 20000
  integer :: nimages = 1
  integer :: local_size

  ! Distributed work data
  real(8), allocatable :: x(:)
  real(8) :: xdot

  ! Determine partition
  nimages = num_images()
  local_size = global_size/nimages

  ! Create smaller vectors
  allocate(x(local_size))

  !-------------------------------------------------------------------!
  ! Test primitive array norm using coarrays
  !-------------------------------------------------------------------!

  test_array : block

    ! Set known values and compute norm
    x = 1.0d0
    xdot = dot_product(x,x)  
    call co_reduce (xdot, operator=sum)
    if (this_image() == 1) then
       write(*,*) "Norm of the array is", sqrt(xdot) , this_image()
    end if

    ! Compute norm with random values in the vector
    CALL RANDOM_NUMBER(x)
    xdot = dot_product(x,x)  
    call co_reduce (xdot, operator=sum)
    if (this_image() == 1) then
       write(*,*) "Norm of the array is", sqrt(xdot) , this_image()
    end if

  end block test_array

  !-------------------------------------------------------------------!
  ! Test derived data type for distributed vector
  !-------------------------------------------------------------------!

  test_datatype : block

    use vector_class, only : vector

    type(vector) :: xvec
    real(8) :: xnorm

    ! create a vector datatype
    xvec = vector(local_size)

    ! set known values and compute norm
    xvec % values = 1.0d0
    xnorm = xvec % norm()
    if (this_image() == 1) then
       write(*,*) "Norm of the vector is", xnorm , this_image()
    end if

    ! set random values and compute norm
    !call xvec % random_values()
    xvec % values(:) = x(:)
    xnorm = xvec % norm()
    if (this_image() == 1) then
       write(*,*) "Norm of the vector is", xnorm , this_image()
    end if

    call xvec % axpy(1.0d0, xvec)
    xnorm = xvec % norm()
    if (this_image() == 1) then
       write(*,*) "Norm of the vector is", xnorm , this_image()
    end if

  end block test_datatype

  deallocate(x)

contains

  ! MPI_SUM
  pure function sum(a, b)
    real(8), intent(in) :: a, b
    real(8) :: sum
    sum = a + b
  end function sum

end subroutine test_norm

program main
  call test_norm
end program main
