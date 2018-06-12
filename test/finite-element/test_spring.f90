!=====================================================================!
! 
!=====================================================================!

program main

  use element_class   , only : spring
  use assembler_class , only : assembler
  use linear_algebra  , only : eigvals, solve

  implicit none

  ! Domain and meshing
  real(8), parameter        :: a = 0.0d0, b = 1.0d0
  integer, parameter        :: nelems = 10
  integer                   :: nnodes = nelems + 1
  integer                   :: me, nimages
  
  real(8)                   :: dx
  integer                   :: lnelems
  integer                   :: e
  integer                   :: leid
  integer                   :: lnode_ids(2)
  real(8)                   :: lnodes(2, 1) = 0.0d0
  real(8), allocatable      :: x(:), rhs(:), mat(:,:), jac(:,:)
  
  type(spring), allocatable :: springs(:)
  type(assembler)           :: spring_assembler
  
  integer :: i,j 

  !-------------------------------------------------------------------!
  ! Mesh the one dimensional geometry and partition between procs
  !-------------------------------------------------------------------!
  
  mesh_partition : block

    ! Mesh sizse (element length)
    dx = (b-a)/dble(nelems)

    ! Identify processors
    me      = this_image()
    nimages = num_images()
    if (nimages .gt. nelems) STOP "Too many processors for the problem size"
    if (mod(nelems, nimages) .ne. 0) STOP "Unequal partition. Use different number of processors."

    ! Decompose domain based on procs
    lnelems = nelems/nimages

  end block mesh_partition

  !-------------------------------------------------------------------!
  ! Create elements associated with the mesh
  !-------------------------------------------------------------------!
  
  elements : block
    
    ! Create elements locally on each processor
    allocate(springs(lnelems))
    do e = 1, lnelems
       leid        = e
       lnode_ids   = [leid, leid + 1]
       lnodes(:,1) = [dble(lnode_ids - 1)*dx]
       springs(e)  = spring(leid, lnode_ids, lnodes, 1.0d0, me)
       call springs(e) % to_string()  
    end do

  end block elements

  !-------------------------------------------------------------------!
  ! Assemble linear system
  !-------------------------------------------------------------------!
  
  assembly : block

    ! Create an assembler for these elements
    spring_assembler = assembler(springs)
    print *, "created assembler"

    ! Create solution, right hand side, jacobian arrays
    allocate(x(spring_assembler % ndof))
    allocate(rhs(spring_assembler % ndof))
    allocate(jac(spring_assembler % ndof, spring_assembler % ndof))

    call random_number(x) ! assumed solution
    x = 1.0d0
    call spring_assembler % residual(x, rhs)
    call spring_assembler % jacobian(x, jac)
    
    print *, "res=", rhs
    do i = 1, spring_assembler % ndof
       print *,  ( jac(i,j), j = 1, spring_assembler % ndof )
    enddo

!!$    print *, "jac-vec pdt 1:", matmul(jac, x)
!!$    call spring_assembler % jacobian_vector_product(x, rhs)
!!$    print *, "jac-vec pdt 2:", rhs

    !print *, "exact sol", solve(jac, rhs)
    
  end block assembly
  
  !-------------------------------------------------------------------!
  ! Solve the linear system
  !-------------------------------------------------------------------!
  
  solver: block

    integer, parameter :: max_it = 100000
    real(8), parameter :: max_tol = 1.0d-8
    integer :: iter, flag, i, j
    real(8) :: tol   
    real(8), allocatable :: x(:)[:]

    allocate(x(spring_assembler % ndof)[*])
    x = 1.0d0

    call spring_assembler % solve(max_it, max_tol, x, iter, tol, flag)

    print *, "solution", x, "from", this_image()
    
    !rhs = 0.0d0
    !call spring_assembler % jacobian_vector_product(x, rhs)
    !print *, rhs

  end block solver

  deallocate(springs)

end program main
