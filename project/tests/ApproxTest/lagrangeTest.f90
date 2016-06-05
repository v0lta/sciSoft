
program LagrangeTest
  use UtilModule
  use LagrangeModule

  integer, parameter :: M = 3
  real(kind = PRECISION), dimension(:) :: testarray1(M)
  real(kind = PRECISION), dimension(M,M) :: testarray2

  print *, 'start'

  testarray2 = 1
  !print *, M, shape(testarray2)
  call genLagrangePols(M,testarray2)
  print *, '---------------------------------done'
  call printRows(testarray2,'(5f8.4)   ')

  !print *, conv((/ 1.0,2.0,3.0,4.0,5.0,6.0,7.0 /), (/ 1.0,2.0 /), 7,2)








end program
