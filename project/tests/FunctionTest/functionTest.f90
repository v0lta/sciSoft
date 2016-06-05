program functionTest
  use FunctionModule
  use UtilModule
  implicit none

  real(kind = PRECISION),dimension(100) :: testX
  real(kind = PRECISION),dimension(100) :: testY
  integer :: i
  type(FunctionType) :: testFunction



  do i = 1,100
    testX(i) = i * 3.14
  end do

  testY = sin(testY)

  print *, 'it works.', PRECISION


  testFunction = FunctionType(testX,testY)

  call printFunction(testFunction)

end program functionTest
