module element_class

  implicit none

  private
  
  public :: spring
  
  type :: spring

     ! Unique global identifier for this element
     integer :: element_id

     ! Global node identifiers associated with this element
     integer :: nnodes
     integer, allocatable :: node_ids(:)
     real(8), allocatable :: nodes(:,:)

     ! DOF calculations
     integer :: ndof_per_node = 1
     integer :: ndof = 0
     real(8) :: k
     integer :: proc_id
     
   contains

     procedure :: add_residual
     procedure :: add_jacobian     
     procedure :: add_jacobian_vector_pdt
     procedure :: to_string
     
  end type spring

  ! Constructor
  interface spring
     module procedure construct_spring
  end interface spring
  
contains

  !===================================================================!
  ! Constructor for spring
  !===================================================================!

  type(spring) function construct_spring( element_id, node_ids, &
       & nodes, spring_constant, proc_id ) result (this)

    integer, intent(in) :: element_id
    integer, intent(in) :: node_ids(:)
    real(8), intent(in) :: nodes(:,:)
    real(8), intent(in) :: spring_constant
    integer, intent(in) :: proc_id

    this % element_id = element_id
    this % node_ids = node_ids
    this % nodes = nodes
    this % k = spring_constant
    this % proc_id = proc_id

    !print *, size(node_ids)
    
    this % ndof = this % ndof_per_node * size(node_ids)
    
  end function construct_spring
  
  !===================================================================!
  ! Print object info
  !===================================================================!
  
  subroutine to_string(this)

    class(spring) :: this

    print *, "eid", this % element_id,  "nids", this % node_ids, "x", this % nodes, "proc", this % proc_id
    
  end subroutine to_string

  !===================================================================!
  ! Add the residual 
  !===================================================================!

  subroutine add_residual(this, x, res)

    class(spring) :: this
    real(8), intent(in) :: x(:)
    real(8), intent(inout) :: res(:)

    !res(1) = res(1) + 1.0d0
    res(2) = res(2) + 1.0d0

    ! use forces may be?
    !res(1) = res(1) + 1.0d0 !this % k * (x(1) - x(2)) 
    !res(2) = res(2) + 1.0d0 !this % k * (x(2) - x(1)) 
    
  end subroutine add_residual

  !===================================================================!
  ! Add the jacobian
  !===================================================================!
  
  subroutine add_jacobian(this, x, jac, bc)

    class(spring) :: this
    real(8), intent(in) :: x(:)
    real(8), intent(inout) :: jac(:,:)
    integer, intent(in) :: bc
    
    if ( bc .eq. 1 ) then
       jac(1,1) = jac(1,1) + this % k + 100.0d0
    else
       jac(1,1) = jac(1,1) + this % k
    end if
       jac(2,1) = jac(2,1) - this % k
    
       jac(1,2) = jac(1,2) - this % k
       jac(2,2) = jac(2,2) + this % k
    
  end subroutine add_jacobian

  !===================================================================!
  ! Add the jacobian-vector product
  !===================================================================!
  
  subroutine add_jacobian_vector_pdt(this, x, pdt, bc)

    class(spring) :: this
    real(8), intent(in) :: x(:)
    real(8), intent(inout) :: pdt(:)
    integer :: bc
    real(8) :: A(this % ndof, this % ndof)
    
    A = 0.0d0
    call this % add_jacobian(x, A, bc)
    pdt = pdt + matmul(A, x)

    !pdt(1) = pdt(1) + this % k * x(1) - this % k * x(2)
    !pdt(2) = pdt(2) - this % k * x(2) + this % k * x(2)
    !print *, "out", pdt
    ! Any communication with neighbor?
    
  end subroutine add_jacobian_vector_pdt
  
end module element_class

!=====================================================================!
! Assembler class for spring elements
!=====================================================================!

module assembler_class

  use element_class

  implicit none

  type :: assembler

     type(spring), allocatable :: elems(:)
     integer :: nelems
     integer :: ndof

   contains

     procedure :: residual
     procedure :: jacobian
     procedure :: jacobian_vector_product
     procedure :: solve

  end type assembler

  ! 
  interface assembler
     module procedure create_assembler
  end interface assembler

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
  
  type(assembler) function create_assembler(elems) result (this)

    type(spring), intent(in) :: elems(:)
    integer :: i

    allocate(this % elems, source = elems)
    this % nelems = size(elems)
    
!!$    do i = 1, this % nelems
!!$       this % ndof = this % ndof + this % elems(i) % ndof
!!$    end do

    ! removing common dofs
    this % ndof = this % nelems + 1
    
    print *, "NDOF=",this % ndof, "proc", this_image()
    
  end function create_assembler
  
  !===================================================================!
  ! Add the residual 
  !===================================================================!
  
  subroutine residual(this, x, res)

    class(assembler) :: this
    real(8), intent(in) :: x(:)
    real(8), intent(inout) :: res(:)
    integer :: nids(2)
    integer :: e

    res = 0.0d0
    
    do e = 1, this % nelems
       nids = this % elems(e) % node_ids
       call this % elems(e) % add_residual( &
            & x(nids(1):nids(2)), &
            & res(nids(1):nids(2)) &
            & )
    end do

    !call this % apply_bc(res)
    
  end subroutine residual

  !===================================================================!
  ! Add the jacobian
  !===================================================================!
  
  subroutine jacobian(this, x, jac)

    class(assembler)       :: this
    real(8), intent(in)    :: x(:)
    real(8), intent(inout) :: jac(:,:)
    integer                :: nids(2)
    integer                :: e

    ! Zero the jacobian matrix entries
    jac = 0.0d0

    ! Add individual blocks to the matrix
    do e = 1, this % nelems      
       nids = this % elems(e) % node_ids     
       associate ( &
            & jmat => jac(nids(1):nids(2), nids(1):nids(2)), &
            & xvec => x(nids(1):nids(2)) ) 
         call this % elems(e) % add_jacobian(xvec, jmat, e)         
       end associate
  end do

  end subroutine jacobian
  
  !===================================================================!
  ! Add the jacobian-vector product
  !===================================================================!
  
  subroutine jacobian_vector_product(this, x, Jx)

    class(assembler) :: this
    real(8), intent(in) :: x(:)
    real(8), intent(inout) :: Jx(:)
    integer :: e
    integer :: nids(this % ndof)

    Jx = 0.0d0
    
    do e = 1, this % nelems

          nids = this % elems(e) % node_ids
          
          call this % elems(e) % add_jacobian_vector_pdt( &
               & x(nids(1):nids(2)), &
               & Jx(nids(1):nids(2)), e &
               & )
       
    end do

    ! 
    !call this % apply_bc(res)
    
    ! Any communication with neighbor?

  end subroutine jacobian_vector_product
  
  subroutine solve(this, max_it, max_tol, x, iter, tol, flag)

    class(assembler) :: this
    integer , intent(in) :: max_it
    real(8), intent(in) :: max_tol

    real(8), intent(inout) :: x(:)
    integer , intent(out)   :: iter
    real(8), intent(out)   :: tol
    integer , intent(out)   :: flag

    ! create local data
    real(8), allocatable :: p(:), r(:), w(:), b(:)
    real(8), allocatable :: rho(:)
    real(8) :: alpha, beta
    real(8) :: bnorm, rnorm

    ! Memory allocations
    allocate(r, p, w, b, mold=x)
    allocate(rho(max_it))

    ! Start the iteration counter
    iter = 1

    call this % residual(x, b)

    ! Norm of the right hand side
    bnorm = co_norm2(b)

    ! Norm of the initial residual
    call this % jacobian_vector_product(x, r)
    r         = b - r !co_matmul(A, x)
    rnorm     = co_norm2(r)
    tol       = rnorm/bnorm
    rho(iter) = rnorm*rnorm

    if (this_image() .eq. 1) then
       open(10, file='cg.log', action='write', position='append')
    end if
    
    ! Apply Iterative scheme until tolerance is achieved
    do while ((tol .gt. max_tol) .and. (iter .lt. max_it))

       ! step (a) compute the descent direction
       if ( iter .eq. 1) then
          ! steepest descent direction p
          p = r
       else
          ! take a conjugate direction
          beta = rho(iter)/rho(iter-1)
          p = r + beta*p
       end if

       ! step (b) compute the solution update
       call this % jacobian_vector_product(p, w)
       !w = co_matmul(A,p)

       ! step (c) compute the step size for update
       alpha = rho(iter)/co_dot_product(p, w)

       ! step (d) Add dx to the old solution
       x = x + alpha*p

       ! step (e) compute the new residual
       r = r - alpha*w
       !r = b - matmul(A, x)

       ! step(f) update values before next iteration
       rnorm = co_norm2(r)
       tol = rnorm/bnorm

       if (this_image() .eq. 1) then
          write(10,*) iter, tol
          print *, iter, tol
       end if
       
       iter = iter + 1

       rho(iter) = rnorm*rnorm

    end do

    close(10)

    deallocate(r, p, w, rho)

    flag = 0

  end subroutine solve

  
  !===================================================================!
  ! Function to compute the dot product of two distributed vector
  !===================================================================!

  function co_dot_product(a, b) result(dot)

    real(8), intent(in) :: a(:), b(:)    
    real(8) :: dot

    ! find dot product, sum over processors, take sqrt and return
    dot = dot_product(a, b)  
    call co_sum (dot)

  end function co_dot_product

end module assembler_class
