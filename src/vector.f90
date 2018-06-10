module vector_class

  use basis_class, only : basis
  implicit none
  
  type :: vector

     type(basis) :: B ! rectangular, spherical, cylindrical
     real(8)  :: coordinates(3)
     
   contains
     
     procedure :: decompose
     
  end type vector

  ! Constructor
  interface vector
     module procedure create_vector
  end interface vector

contains

  function decompose(this, B) result(q)
    
    class(vector), intent(in) :: this
    type(basis), intent(in) :: B
    real(8) :: q(3)
    real(8) :: qtmp(3)
    
    call this % B % basis_vector_product(this % coordinates, qtmp)    
    call B % transposed_basis_vector_product(qtmp, q)

  end function decompose
  
  type(vector) function create_vector(coordinates, base) result (this)

    real(8), intent(in) :: coordinates(3)
    type(basis), intent(in) :: base

    this % coordinates = coordinates
    this % B = base

  end function create_vector
  
end module vector_class

!=====================================================================!
!
!=====================================================================!

program test_vector

  use vector_class, only : vector
  use basis_class, only : basis

  type(basis) :: rect_basis
  type(basis) :: spherical
  type(basis) :: cylindrical
  type(vector) :: v

  cylindrical = basis() 
  v = vector([1.0d0, 2.0d0, 3.0d0], cylindrical)
  print *, v % decompose(cylindrical)
  
end program test_vector
