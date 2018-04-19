!=====================================================================!
! Program to compute the matrix vector product using fortran coarrays.
!
! Author : Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program test_matmul

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

  ! Multiply as you would normally
  b = co_matmul(A, x)

  ! call comatmul(A, x, btmp, b)
  
  write(*,*) "the new vector is", b , this_image()

  deallocate(A,x,b,btmp)
  
contains
  
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

  !===================================================================!
  ! Subroutine that computes the matrix vector product in a
  ! distributed fashion for columnwise decomposition of matrix
  !===================================================================!
  
  subroutine comatmul(A, x, work, b)

    ! Arguments
    real(8), intent(in) :: A(:,:)
    real(8), intent(in) :: x(:)
    real(8), intent(inout) :: work(:)
    real(8), intent(inout) :: b(:)
    
    ! Local variables
    integer :: nimages
    integer :: me, local_size
    integer :: stat
    character(10) :: msg

    ! Determine partition
    nimages = num_images()
    me = this_image()
    local_size = size(x)

    ! Multiply, sum and distrbute
    work = matmul(A,x)
    call co_sum(work, stat=stat, errmsg=msg)
    b = work((me-1)*local_size+1:me*local_size)

  end subroutine comatmul

end program test_matmul
