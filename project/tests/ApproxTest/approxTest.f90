program approxTest
  use LagrangeModule
  use UtilModule
  use FunctionModule

  type(FunctionType) :: testFunction
  real(kind = PRECISION), dimension(:) :: desValues(63)
  real(kind = PRECISION), dimension(:) :: functionValues(26)
  real(kind = PRECISION), dimension(:) :: xValues(26)
  integer :: funSize = 26
  integer :: desSize = 63

  call linspace(-2*PI,2*PI,xValues)
  functionValues = cos(xValues)
  call linspace(-4*PI,4*PI,desValues)

  testFunction = FunctionType(funSize,functionValues,xValues,desSize,desValues)
  call sortFunction(testFunction)
  print '(es24.16)', lagrange(testFunction)


end program
