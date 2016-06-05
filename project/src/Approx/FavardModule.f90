module FavardModule
  use UtilModule
  use IntegrationModule
  implicit none

contains

    subroutine favard(N, walpha, wbeta, polyMat)
      integer, intent(in) :: N
      real(kind = PRECISION), intent(in) :: wbeta, walpha
      real(kind = PRECISION), dimension(:,:), intent(inout) :: polyMat
      real(kind = PRECISION), dimension(N+1) :: convRes
      real(kind = PRECISION), dimension(N+1,2) :: ab
      integer :: i
      procedure (funPointer), pointer :: funp => null()
      
      !nagfor and ifort cannot initilaize pointers at creation.
      funp => normFun

      ab = rJacobi(N,walpha,wbeta)
      convRes = 0
      polyMat = 0
      polyMat(N,1) = 1
      convRes = conv((/ 1.0, ab(1,1) /),polyMat(:,1),2,N)
      polyMat(:,2) = convRes(2:N+1)

      !generate the orthogonal polynomial.
      do i = 3,N
        convRes = conv((/ 1.0, ab(i-1,1) /), polyMat(:,i-1),2,N)
        convRes = convRes - (/ 0.0, ab(i-1,2)*polyMat(:,i-2) /)
        polyMat(:,i) = convRes(2:N+1)
      end do

      pAlpha = wBeta
      pBeta = wBeta
      !normalize the polynomial to get an orthonormal basis.
      do i = 1,N
        pol1 = polyMat(:,i)
        polyMat(:,i) = polyMat(:,i) / sqrt(quadIntegral(funp,-1.0,1.0))
      end do
    end subroutine favard

    function normFun(x) result(res)
      real(kind = PRECISION), intent(in) :: x
      real(kind = PRECISION) :: res
      res = omega(x) * polyval(pol1,x)**2
    end function normFun

    function omega(x) result(res)
      real(kind = PRECISION), intent(in) :: x
      real(kind = PRECISION) :: res
      !use the public alpha and beta for quadpack interface conformance.
      res = (1-x)**pAlpha * (1+x)**pBeta
    end function omega

end module FavardModule
