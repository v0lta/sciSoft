program utilTest
  use UtilModule
  implicit none

  real(kind = PRECISION), dimension(3) :: poly
  real (kind = PRECISION) :: res

  poly = (/ 3, 2, 1 /)
  print *,  polyval(poly,2.5)


end program
