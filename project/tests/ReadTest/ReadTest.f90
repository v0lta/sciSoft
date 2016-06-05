program valueInTest
  use ReadModule
  use UtilModule
  use FunctionModule
  implicit none
  integer :: i

  type(FunctionType) :: testFunction

  testFunction = readFunctionData()

  print *, size(testFunction%functionValues)
  print '(es24.16)', testFunction%functionValues

  if (size(testFunction%desValues) .ne. 0) then
    print '(es24.16)', testFunction%xValues
    print *, size(testFunction%desValues)
    print '(es24.15)', testFunction%desValues
  end if

end program
