program GaussTest
  use IntegrationModule
  use UtilModule
  implicit none

  integer, parameter :: N = 100
  real(kind = PRECISION), dimension(N+1,2) :: ab
  real(kind = PRECISION) :: alpha, beta, sum
  real(kind = PRECISION), dimension(:) :: fun(N), x(N)
  integer :: i

  alpha = 0.0
  beta = 0.0

  call linspace(-1.0,1.0,x)

  do i = 1,N
    fun(i) = x2(x(i))
  end do

  !function integral(alpha, beta, N, funp, x, intFunction) result(sum)
  sum = integral(alpha,beta,N, funValues = fun )
  print '(f16.12)', sum

  sum = trapz(N, x, fun )
  print '(f16.12)', sum

contains

  function x2(x) result(res)
    real(kind = PRECISION) :: x,res
    res = x**2
  end function x2

end program
