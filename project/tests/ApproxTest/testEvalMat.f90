program gramSchmidtTest
  use GramSchmidtModule
  use UtilModule


  integer, parameter :: M = 6
  real(kind = PRECISION) :: alpha, beta
  real(kind = PRECISION), dimension(M,M) :: polyMat
  real(kind = PRECISION), dimension(:,:) :: evalMat(2,M)
  integer :: i

  polyMat = reshape((/ 0.00, 0.00, 0.00, 0.00, 0.00, 0.64,&
                       0.00, 0.00, 0.00, 0.00, 1.32, 0.55,&
                       0.00, 0.00, 0.00, 2.67, 1.25,-0.57,&
                       0.00, 0.00, 5.37, 2.61,-2.41,-0.57,&
                       0.00, 10.77, 5.34,-7.44,-2.44,0.56,&
                       10.48, 5.27, -9.81, -3.68, 1.69, 0.28 /), shape(polyMat))
  print *, (/ ('-',i=1,70) /), 'evalMat'
  call printRows(polyMat,'(6f6.2)   ')

  call evalMatrix(polyMat, (/ 1.0, -1.0 /), M, evalMat)

  print *, (/ ('-',i=1,70) /), 'evalMat'
  call printRows(evalMat,'(6f6.2)   ')

  !print *, (/ ('-',i=1,70) /), 'polyVal'
  !print *, polyval((/ 1.0,2.0 /), -1.0)


  !P =  0.6400    1.8700    3.3500    5.0000    6.7900
  !     0.6400   -0.7700    0.8500   -0.9200    0.9900
  !That works.


end program
