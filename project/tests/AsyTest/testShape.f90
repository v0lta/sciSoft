program testShape
  use UtilModule
  use AsyModule
  implicit none

  real   (kind = PRECISION), dimension(1,1,2) ::  g
  real   (kind = PRECISION), dimension(2,2,2) ::  h
  real   (kind = PRECISION), dimension(2,2,2) ::  h2

  print *, (/1,2 /)
  g = reshape((/ 1, 2/), (/ 1,1,2 /))
  print *, g(1,1,1), g(1,1,2)
  print * , '---------------------------'
  h = repmat(g,2)
  print *, h(:,:,1)
  print *, h(:,:,2)
  print *, sum(h,3)

end program testShape
