program asyTest
  use UtilModule
  use AsyModule
  implicit none


  integer, parameter :: maxOrder = 2
  integer, parameter :: maxOrderE = 4
  integer, parameter :: mo = nint(maxOrder-2/2.0)

  real(kind = PRECISION) :: test, alpha, beta, Dinf = 0
  integer :: n,k = 0
  complex(kind = PRECISION), dimension(2,2,1,1) :: WV
  complex(kind = PRECISION), dimension(4) :: c
  complex(kind = PRECISION), dimension(4) :: d
  complex(kind = PRECISION), dimension(2,2,maxOrder-1,ceiling((maxOrder-1.0)/2.0)) :: Uleft
  complex(kind = PRECISION), dimension(2,2,maxOrder-1,ceiling((maxOrder-1.0)/2.0)) :: Uright
  complex(kind = PRECISION), dimension(2,2,maxOrderE-1,ceiling((maxOrderE-1.0)/2.0)) :: UleftE
  complex(kind = PRECISION), dimension(2,2,maxOrderE-1,ceiling((maxOrderE-1.0)/2.0)) :: UrightE

  c = 0
  d = 0
  Dinf = 1

  call compU(alpha, beta, Dinf, c ,d, maxOrder, Uleft, Uright)
  call compUExplicit(alpha, beta, Dinf, c, d, UrightE, UleftE)

  print *, '------Uleft error first order-----'
  print *, abs(norm2(abs(UleftE(:,:,1,1))-abs(UleftE(:,:,1,1))))
  print *, '------Uright error first order---'
  print *, abs(norm2(abs(UrightE(:,:,1,1))-abs(UrightE(:,:,1,1))))

end program asyTest
