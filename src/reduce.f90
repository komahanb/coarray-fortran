program test
  real(8) :: val
  val = DBLE(this_image ())
  call co_reduce (val, operator=sum)! result_image=1, 
  if (this_image() == 1) then
     write(*,*) "Product value", val , this_image() ! prints num_images() factorial
  end if
contains
  pure function sum(a, b)
    real(8), intent(in) :: a, b
    real(8) :: sum
    myprod = a + b
  end function sum
end program test
