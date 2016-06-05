program convTest
  use UtilModule
  implicit none

  integer :: i
  integer, parameter :: s1 = 3
  integer, parameter :: s2 = 3
  real(kind = PRECISION), dimension(s1) :: poly1
  real(kind = PRECISION), dimension(s2) :: poly2

  read *, poly1
  read *, poly2


  print '(f10.5)', conv(poly1,poly2,s1,s2)


end program
