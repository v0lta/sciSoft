program GaussTest
  use GramSchmidtModule
  use IntegrationModule
  use UtilModule
  use FunctionModule
  implicit none

  integer, parameter :: N = 100
  real(kind = PRECISION), dimension(N+1,2) :: ab
  real(kind = PRECISION) :: alpha, beta, sum
  real(kind = PRECISION), dimension(N) :: x
  integer :: i = 1

  procedure (funPointer), pointer :: funp => null()

  funp => x2
  !print *, funp(2.0)
  print *, trapz(funp,-1.0,1.0,N)

  !use the gauss integrator
  alpha = 0
  beta = 0

  call linspace(-1.0,1.0,x)
  print *, integral(alpha, beta, N, funp, x)

contains


  function x2(x) result(res)
    real, intent(in) :: x
    real :: res

    res = x**2
  end function
end program
