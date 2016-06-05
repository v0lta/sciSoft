program trapzTest
  use IntegrationModule
  use UtilModule
  implicit none

  integer, parameter :: N = 1000
  integer :: i
  real(kind = PRECISION),dimension(N) :: x,f


  call linspace(-1.0,1.0,x)
  do i = 1,N
    f(i) = x2(x(i))
  end do
  print *, trapz(N,x,f)

contains

  function x2(x) result(res)
    real(kind = PRECISION) :: x,res
    res = x**2
  end function x2



end program trapzTest
