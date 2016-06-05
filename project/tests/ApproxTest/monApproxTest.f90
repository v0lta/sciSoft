program gramSchmidtTest
  use GramSchmidtModule
  use UtilModule
  use FunctionModule
  use ReadModule

  integer, parameter :: M = 4
  !trapz points.
  integer, parameter :: N = 1000
  integer :: i
  real(kind = PRECISION) :: galpha, gbeta, walpha, wbeta
  real(kind = PRECISION), dimension(M,M) :: polyMat
  type(FunctionType)  :: f
  real(kind = PRECISION) :: transA,transB
  real(kind = PRECISION),dimension(:) :: c(M)
  real(kind = PRECISION), dimension(:), allocatable :: result
  real(kind = PRECISION), dimension(:,:),allocatable :: evalMat

  walpha = -0.3
  wbeta = 0.3
  galpha = 0
  gbeta = -0.4

  !Read in data
  f = readFunctionData()
  call sortFunction(f)
  allocate(evalMat(f%funSize,M))
  allocate(result(f%desSize))

  !Generate a non orthonormal monomial basis.
  call eye(polyMat)
  call flipud(polyMat)

  !call printRows(polyMat)

  !call gramSchmidt(N,gAlpha,gBeta,wAlpha,wBeta,polyMat,0)
  call gramSchmidt(f%funSize,gAlpha,gBeta,wAlpha,wBeta,polyMat,1)
  print '(4f6.1)', polyMat

  call transformX(f%xValues,f%desValues, transA, transB)
  !Transfomrmation for polyval...

  call evalMatrix(polyMat, f%xValues, M, evalMat)
  !call printRows(evalMat,'(6f8.4)   ')

  c = compC(evalMat,f,f%funSize, M)

  !print *, (/ ('-',i=1,70) /), 'c'
  !print '(6f10.4)', c
  !print *, (/ ('-',i=1,70) /)

  deallocate(evalMat)
  allocate(evalMat(f%desSize,M))
  call evalMatrix(polyMat, f%desValues, M, evalMat)

  !print *, f%desValues
  !call printRows(evalMat,'(6f8.4)   ')

  call reverseTransformX(f%xValues, f%desValues, transA, transB)


  result = matmul(evalMat,c)
  !print *, (/ ('-',i=1,70) /), 'result'
  print '(9f10.4)', result
  !print *, (/ ('-',i=1,70) /), 'result'
  deallocate(result)
  deallocate(evalMat)
end program
