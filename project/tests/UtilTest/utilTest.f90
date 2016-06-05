program utilTest
  use UtilModule
  implicit none

  real(kind = PRECISION), dimension(100) :: toSort

  read *, toSort
  call sort(toSort)
  print '(f14.6)',toSort


end program
