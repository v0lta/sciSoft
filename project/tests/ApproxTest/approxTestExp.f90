program approxTest
  use LagrangeModule
  use UtilModule
  use FunctionModule

  type(FunctionType) :: testFunction
  real(kind = PRECISION), dimension(:) :: desValues(1)
  real(kind = PRECISION), dimension(:) :: functionValues(26)
  real(kind = PRECISION), dimension(:) :: xValues(26)
  integer :: funSize = 26
  integer :: desSize = 1

  call linspace(-2*PI,2*PI,xValues)
  functionValues = exp(xValues)
  desValues = -PI/16;

  testFunction = FunctionType(funSize,functionValues,xValues,desSize,desValues)

  print '(f12.5)', lagrange(testFunction)


end program
