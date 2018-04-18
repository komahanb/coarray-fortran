!=====================================================================!
! Program to compute the matrix vector product using fortran coarrays.
!
! Author : Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program test_norm

  implicit none
  
  integer, parameter :: global_size = 4

  ! Distributed work data
  real(8), allocatable :: A(:,:)
  real(8), allocatable :: x(:)
  real(8), allocatable :: b(:)

  integer :: nimages = 1
  integer :: local_size
  integer :: i, me
  
  ! Determine partition
  nimages    = num_images()
  local_size = global_size/nimages
  me         = this_image()

  ! Decompose matrix columnwise and vector rowwise
  allocate(A(global_size,local_size))
  allocate(x(local_size))

!!$  print *, shape(A), this_image()
!!$  print *, shape(x), this_image()
  
  ! Fill the local arrays with random values or known entries
  call random_number(A)
  call random_number(x)

!!$  if (this_image() == 1) then
!!$     A(:,1) = [1.0d0, 2.0d0, 3.0d0, 4.0d0]
!!$     A(:,2) = 2.0d0*[1.0d0, 2.0d0, 3.0d0, 4.0d0]
!!$     x(1:2) = [1.0d0, 2.0d0]
!!$  else
!!$     A(:,1) = 3.0d0*[1.0d0, 2.0d0, 3.0d0, 4.0d0]
!!$     A(:,2) = 4.0d0*[1.0d0, 2.0d0, 3.0d0, 4.0d0]
!!$     x(1:2) = [3.0d0, 4.0d0]
!!$  end if

  ! Matrix vector multiplication
  b = matmul(A, x)
  print *, "A=", A, "from", me
  print *, "x=", x, "from", me
  print *, "b=", b, "from", me
  
  ! Compiler does not support array reduction yet!
   do i = 1, global_size
     call co_reduce (b(i), operator=sum)
  end do

  ! print b = Ax vector
  if (this_image() == 1) then
     write(*,*) "the new vector is", b , this_image()
  end if
  
contains

  ! MPI_SUM
  pure function sum(a, b)
    real(8), intent(in) :: a, b
    real(8) :: sum
    sum = a + b
  end function sum
  
end program test_norm
