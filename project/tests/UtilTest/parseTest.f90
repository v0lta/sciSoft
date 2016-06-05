program parseTest
  use ParseModule

  !part 1
  integer :: i, N, M
  character(len=32) :: arg
  logical  :: hlp, tst, lag, tme, int, intGaus, quad, pNoS, ws, degS, orth, proj
  real  :: wAlpha, wBeta, gAlpha, gBeta, x
  !part 2
  logical :: evOrth, isGM, isFav, isAsy
  integer :: T, h


  call readCmdLine(hlp, tst, lag, tme, int, intGaus,quad, pNoS, ws, degS, orth, &
                  proj, wAlpha, wBeta, gAlpha, gBeta, i, N, M, evOrth, h, x, &
                  isGM, isFav, isAsy, T)


end program parseTest
