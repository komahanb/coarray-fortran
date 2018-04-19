!=====================================================================!
! Program to compute the matrix vector product using fortran coarrays.
!
! Author : Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program test_norm

  implicit none
  
  integer, parameter :: global_size = 40

  ! Distributed work data
  real(8), allocatable :: A(:,:)
  real(8), allocatable :: x(:)
  real(8), allocatable :: b(:)
  real(8), allocatable :: btmp(:)

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
  allocate(b(local_size))
  allocate(btmp(global_size))

  ! Fill the local arrays with random values or known entries
  call random_number(A)
  call random_number(x)

  ! Matrix vector multiplication
  btmp = matmul(A, x)

  ! Sum btmp across all processors
  call co_sum(btmp)

  b = btmp((me-1)*local_size+1:me*local_size)
  !write(*,*) "the new vector is", b , this_image()

  deallocate(A,x,b,btmp)
  
contains

  ! MPI_SUM
  pure function sum(a, b)
    real(8), intent(in) :: a, b
    real(8) :: sum
    sum = a + b
  end function sum
  
end program test_norm
