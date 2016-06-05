program lensTest
  use UtilModule
  use AsyModule
  implicit none


  integer, parameter :: maxOrder = 4
  integer, parameter :: pointNumber = 4
  integer, parameter :: mo = nint(maxOrder-2/2.0)
  real(kind = PRECISION), parameter :: alpha = 0.2_PRECISION
  real(kind = PRECISION), parameter :: beta = -0.2_PRECISION

  real(kind = PRECISION) :: test, Dinf = 0
  integer :: i,j,k,n = 0
  complex(kind = PRECISION), dimension(2,2,1,1) :: WV
  complex(kind = PRECISION), dimension(4) :: c
  complex(kind = PRECISION), dimension(4) :: d
  complex(kind = PRECISION), dimension(2,2,maxOrder-1,ceiling((maxOrder-1.0)/2.0)) :: Uleft = 0
  complex(kind = PRECISION), dimension(2,2,maxOrder-1,ceiling((maxOrder-1.0)/2.0)) :: Uright = 0
  complex(kind = PRECISION), dimension(pointNumber, maxOrder) :: polyMat
  complex(kind = PRECISION) :: ipiVal
  real(kind = PRECISION), dimension(pointNumber) :: z

  c = 0
  d = 0
  Dinf = 2**(-alpha/2.0-beta/2.0);
  call compUExplicit(alpha, beta, Dinf, c , d, Uright, Uleft)

  n = 6
  call linspace(-0.9,0.9,z)
  Dinf = 2**(-alpha/2.0-beta/2.0);
  do i = 1,pointNumber
        do j = 1,maxOrder
        call asyLens(n,j,cmplx(z(i),0.0),alpha,beta,Dinf,Uright,Uleft,ipiVal)
        polyMat(i,j) = ipiVal
        end do
  end do
  !call printRows(cmplxMatrix=polyMat)
  print '(4f10.4)', real(polyMat)
  !do i = 1,size(polyMat,1)
  !      do j = 1,size(polyMat,2)
  !        print '(f10.4)', real(polyMat(i,j))
  !        print '(f10.4)', aimag(polyMat(i,j))
  !      end do
  !end do

end program lensTest
