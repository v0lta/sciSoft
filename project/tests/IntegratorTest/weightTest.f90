program weightTest
  use IntegrationModule
  use UtilModule
  implicit none

  integer, parameter :: N = 61
  real(kind = PRECISION), dimension(N+1,2) :: ab
  real(kind = PRECISION) :: alpha, beta
  integer :: i = 1

  alpha = 0.4
  beta = -0.13

  ab = rJacobi(N,alpha,beta)

  do while (i < N+2)

    print '(4es12.4)', ab(i,1) , ab(i,2)
    i = i + 1
  end do

end program
