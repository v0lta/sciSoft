program quadTest
  implicit none
  abstract interface
    function funp (x)
      real :: funp
      real, intent (in) :: x
    end function funp
  end interface

  procedure (funp), pointer :: F => null ()

  real :: A,ABSERR,B,EPSABS,EPSREL,RESULT
  integer :: IER,IWORK,LAST,NEVAL

  !integration limits a,b
  A = -1.0
  B = 1.0
	!absolute and relative accuracy requested.
  EPSABS = 0.0E0
  EPSREL = 1.0E-3

  F => x2

  CALL QAGS(F,A,B,EPSABS,EPSREL,RESULT,ABSERR, &
            NEVAL,IER)!,LIMIT,LENW,LAST,IWORK,WORK)

  print '(f10.8)', RESULT

contains

  function x2(x) result(res)
    real, intent(in) :: x
    real :: res
    res = x**2
  end function
end program
