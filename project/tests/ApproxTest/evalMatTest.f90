program evalMatTest
  use GramSchmidtModule
  implicit none

  integer, parameter :: M = 6
  real, dimension(M,M) :: polyMat = 1
  real, dimension(3,M) :: evalMat
  real, dimension(3)   :: desValues = (/ -PI/16, 0.0 , PI/16 /);


  call evalMatrix(polyMat, desValues, M, evalMat)
  print '(3f6.3)', evalMat


end program
