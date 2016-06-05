program asyTest
  use UtilModule
  use AsyModule
  implicit none


  integer, parameter :: maxOrder = 4
  integer, parameter :: mo = nint(maxOrder-2/2.0)
  real(kind = PRECISION), parameter :: alpha = 0
  real(kind = PRECISION), parameter :: beta = 0
  real(kind = PRECISION) :: test, Dinf = 0
  integer :: n,k = 0
  complex(kind = PRECISION), dimension(2,2,1,1) :: WV
  complex(kind = PRECISION), dimension(4) :: c
  complex(kind = PRECISION), dimension(4) :: d
  complex(kind = PRECISION), dimension(2,2,maxOrder-1,ceiling((maxOrder-1.0)/2.0)) :: Uleft
  complex(kind = PRECISION), dimension(2,2,maxOrder-1,ceiling((maxOrder-1.0)/2.0)) :: Uright

  c = 0
  d = 0
  Dinf = 2**(-alpha/2.0-beta/2.0);
  call compUExplicit(alpha, beta, Dinf, c, d, Uright, Uleft)

  print *, '------Uright---(:,:,1,1)'
  call printRows(cmplxMatrix=Uright(:,:,1,1))
  print *, '------Uleft----(:,:,1,1)'
  call printRows(cmplxMatrix=Uleft(:,:,1,1))
  print *, '------Uright---(:,:,2,1)'
  call printRows(cmplxMatrix=Uright(:,:,2,1))
  print *, '------Uleft----(:,:,2,1)'
  call printRows(cmplxMatrix=Uleft(:,:,2,1))
  print *, '------Uright---(:,:,3,1)'
  call printRows(cmplxMatrix=Uright(:,:,3,1))
  print *, '------Uleft----(:,:,3,1)'
  call printRows(cmplxMatrix=Uleft(:,:,3,1))
  print *, '------Uright---(:,:,3,2)'
  call printRows(cmplxMatrix=Uright(:,:,3,2))
  print *, '------ Uleft---(:,:,3,2)'
  call printRows(cmplxMatrix=Uleft(:,:,3,2))


  !call printRows(transpose(reshape((/ 1.0,2.0,3.0,4.0 /), (/2,2/))))

end program asyTest
