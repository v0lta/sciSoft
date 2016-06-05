program projTest
  use LagrangeModule
  use UtilModule
  use FunctionModule
  use ProjectionModule
  use ReadModule

  integer, parameter :: M = 3
  integer :: intPoints, i
  real(kind = PRECISION) :: gAlpha, gBeta, wAlpha, wBeta, curEps
  real(kind = PRECISION), dimension(M,M) :: polyMat, eyeMat
  real(kind = PRECISION), dimension(:,:),allocatable :: evalMat
  type(FunctionType)  :: f
  type(FunctionType)  :: solFunction
  real(kind = PRECISION) :: transA,transB
  real(kind = PRECISION),dimension(M) :: c
  real(kind = PRECISION), dimension(:), allocatable :: result

  gAlpha = 0.4
  walpha = 0
  gBeta = -0.13
  wbeta = -0

  !call genLagrangePols(M,polyMat)

  !Read in data
  f = readFunctionData()
  call sortFunction(f)
  allocate(result(f%desSize))

  !Generate a non orthonormal monomial basis.
  call eye(polyMat)
  call flipud(polyMat)
  call gramSchmidt(f%funSize,gAlpha,gBeta,wAlpha,wBeta,polyMat,2)
  call printRows(polyMat)

  call transformX(f%xValues,f%desValues, transA, transB)
  !curEps = epsilon(f%xValues)
  !result = sqrt(curEps) - 1 + ((f%desValues - minval(f%desValues)) &
  !        / (maxval(f%desValues) - minval(f%desValues))) &
  !        * (2 - 2 * sqrt(curEps))

  c = projC(gAlpha, gBeta, wAlpha, wBeta, polyMat,f,M,0)

  print *, (/ ('-',i=1,70) /), 'c'
  print '(6f10.4)', c
  print *, (/ ('-',i=1,70) /)

  allocate(evalMat(f%desSize,M))

  call evalMatrix(polyMat, f%desValues, M, evalMat)
  print *, size(result)

  result = matmul(evalMat,c)
  !print *, (/ ('-',i=1,70) /), 'result'
  print '(f10.4)', result
  !print *, (/ ('-',i=1,70) /), 'result'
  deallocate(result)
  deallocate(evalMat)

end program
