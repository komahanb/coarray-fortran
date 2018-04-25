program main

  implicit none

  real :: a[*]
  type(real) :: b[*]

  integer :: image
  
  a = this_image()
  b = this_image()*2.0


  print *, a, b
  if (this_image() .eq. 1) then

     do image  = 1, num_images()
        print *, 'Image', this_image(), a[image], b[image]
     end do
     
  end if
  
end program main
