program lagApproxTest
  use GramSchmidtModule
  use UtilModule
  use FunctionModule
  use ReadModule
  use LagrangeModule

  integer, parameter :: M = 4
  integer :: intPoints, i
  real(kind = PRECISION) :: gAlpha, gBeta, wAlpha, wBeta
  real(kind = PRECISION), dimension(M,M) :: polyMat
  real(kind = PRECISION), dimension(M,M) :: eyeMat
  real(kind = PRECISION), dimension(:,:),allocatable :: evalMat
  type(FunctionType)  :: f
  type(FunctionType)  :: solFunction
  real(kind = PRECISION) :: transA,transB
  real(kind = PRECISION),dimension(:) :: c(M)
  real(kind = PRECISION), dimension(:), allocatable :: result

  !Generate non orthonormal lagrange polynomials.
  polyMat = 1
  call genLagrangePols(M,polyMat)

  walpha = 0.1
  wbeta = 1.5
  galpha = 0
  gbeta = -0.4

  !Read in data
  f = readFunctionData()
  call sortFunction(f)
  allocate(evalMat(f%funSize,M))
  allocate(result(f%desSize))

  print *, f%funSize,gAlpha,gBeta,wAlpha,wBeta,polyMat,1
  call gramSchmidt(f%funSize,gAlpha,gBeta,wAlpha,wBeta,polyMat,1)
  call transformX(f%xValues,f%desValues, transA, transB)

  !Transfomrmation for polyval...
  call evalMatrix(polyMat, f%xValues, M, evalMat)
  c = compC(evalMat,f,f%funSize, M)
  print '(6f8.4)', c

  deallocate(evalMat)
  allocate(evalMat(f%desSize,M))
  call evalMatrix(polyMat, f%desValues, M, evalMat)

  call reverseTransformX(f%xValues, f%desValues, transA, transB)


  result = matmul(evalMat,c)
  !print *, (/ ('-',i=1,70) /), 'result'
  print '(f10.4)', result
  !print *, (/ ('-',i=1,70) /), 'result'
  deallocate(result)
  deallocate(evalMat)
end program
