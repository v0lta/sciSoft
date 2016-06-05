program GaussTest
  use ReadModule
  use IntegrationModule
  use UtilModule
  use FunctionModule
  implicit none

  integer, parameter :: N = 61
  real(kind = PRECISION), dimension(N+1,2) :: ab
  real(kind = PRECISION) :: alpha, beta, sum
  integer :: i = 1
  type(FunctionType) :: testFunction

  testFunction = readFunctionData()

  alpha = 0.4_PRECISION
  beta = -0.13_PRECISION


  !function integral(alpha, beta, N, funp, x, intFunction) result(sum)
  sum = integral(alpha,beta,N, funValues = testFunction%functionValues )
  print '(f16.12)', sum

end program
