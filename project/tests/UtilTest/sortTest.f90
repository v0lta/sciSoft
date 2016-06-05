program sortTest
  use UtilModule
  real(kind = PRECISION), dimension(5) :: test

  test = (/ 2.0 ,-1.5, 0.2, -0.327, 5.0 /)
  call sort(test)
  print *, test



end program sortTest
