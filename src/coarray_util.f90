module coarray_util

  implicit none

  public :: co_norm2, co_dot_product

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

end module coarray_util
