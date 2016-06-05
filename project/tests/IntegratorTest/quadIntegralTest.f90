program quadIntegralTest
  use UtilModule
  use IntegrationModule
  implicit none


  real :: res
  procedure (funPointer), pointer :: f => omega
  real(kind = PRECISION) :: a,b

  !set the public variables.
  pAlpha = 0
  pBeta = 0
  !integration limits a,b
  a = -1.0
  b = 1.0
  res = quadIntegral(f,a,b)
  print '(f10.8)',res

contains

  function omega(x) result(res)
    real(kind = PRECISION), intent(in) :: x
    real(kind = PRECISION) :: res
    res = (1+x)**pAlpha * (1-x)**pBeta
  end function

  function x2(x) result(res)
    real(kind = PRECISION), intent(in) :: x
    real(kind = PRECISION) :: res
    res = x**2.0
  end function
end program quadIntegralTest
