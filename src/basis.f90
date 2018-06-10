module basis_class

  implicit none

  private
  public :: basis

  ! Derived type for basis
  type :: basis
     
   contains
     
     procedure :: basis_vector_product
     procedure :: transposed_basis_vector_product
     
  end type basis

contains

   subroutine basis_vector_product(this, q, qprime)

    class(basis), intent(in) :: this
    real(8), intent(in)  :: q(3)
    real(8), intent(out) :: qprime(3)

    real(8) :: B(3,3)
    
    associate(r => 1.0d0, theta => 0.0d0, phi => 0.0d0)
      B(:,1) = [cos(theta)*sin(phi), sin(theta)*sin(phi), cos(phi)]
      B(:,2) = [-sin(theta), cos(theta), 0.0d0]
      B(:,3) = [cos(theta)*cos(phi), sin(theta)*cos(phi), -sin(phi)]
    end associate

    print *, B
    
    qprime = matmul(B, q)
    
!!$    associate(r => q(1), theta => q(2), z => q(3))
!!$      B(:,1) = [cos(theta), sin(theta), 0.0d0]
!!$      B(:,2) = [-sin(theta), cos(theta), 0.0d0]
!!$      B(:,3) = [0.0d0, 0.0d0, 1.0d0]
!!$    end associate

  end subroutine basis_vector_product

  
  pure subroutine transposed_basis_vector_product(this, q, qprime)

    class(basis), intent(in) :: this
    real(8), intent(in)  :: q(3)
    real(8), intent(out) :: qprime(3)

    real(8) :: B(3,3)
    
    associate(r => 1.0d0, theta => 0.0d0, phi => 0.0d0)
      B(:,1) = [cos(theta)*sin(phi), sin(theta)*sin(phi), cos(phi)]
      B(:,2) = [-sin(theta), cos(theta), 0.0d0]
      B(:,3) = [cos(theta)*cos(phi), sin(theta)*cos(phi), -sin(phi)]
    end associate

    qprime = matmul(transpose(B), q)
    
!!$    associate(r => q(1), theta => q(2), z => q(3))
!!$      B(:,1) = [cos(theta), sin(theta), 0.0d0]
!!$      B(:,2) = [-sin(theta), cos(theta), 0.0d0]
!!$      B(:,3) = [0.0d0, 0.0d0, 1.0d0]
!!$    end associate

  end subroutine transposed_basis_vector_product
  
end module basis_class
