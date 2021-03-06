program checkOrthogonality
  use UtilModule
  use GramSchmidtModule
  use IntegrationModule
  use LagrangeModule

  integer, parameter :: M = 3
  integer, parameter :: N = 500
  integer :: i
  real(kind = PRECISION) :: walpha, wbeta, galpha, gbeta
  real(kind = PRECISION), dimension(M,M) :: polyMat
  real(kind = PRECISION), dimension(:)  :: row(M)
  logical :: pass

  walpha = 0.0
  wbeta  = 0.0
  galpha = 0.0
  gbeta  = 0.0

  !A monomial basis!
  call eye(polyMat)
  call flipud(polyMat)

  !Lagrange polynomials.
  !call genLagrangePols(M,polyMat)


  !call gramSchmidt(N,galpha,gbeta,walpha,wbeta,polyMat,0)
  call gramSchmidt(N,galpha,gbeta,walpha,wbeta,polyMat,2)
  !print '(3f10.5)',polyMat
  !print *, (/ ('-',i=1,70) /), 'done'
  !call printRows(polyMat)

  call checkBasis(walpha,wbeta,polyMat,N,pass)
  print *, pass
  !print *, '-------------------------- done'
  !print '(6f10.4)', polyMat


contains


!Also compare with octave file ./tests/ApproxTest/matlab/checkOrth.m
  subroutine checkBasis(alpha,beta,polyMat,N,pass)
    logical,intent(out) :: pass
    integer :: N
    real(kind = PRECISION), dimension(:,:) :: polyMat
    real(kind = PRECISION), dimension(:) :: x(N), funVals(N)
    real(kind = PRECISION) :: res, alpha, beta, eps
    integer :: i,j,k

    pass = .true.
    call linspace(-1.0,1.0,x)
    eps = sqrt(epsilon(x))
    !Transform
    x = (x - minval(x))/(maxval(x) - minval(x)) * (2 - 2*eps) + eps - 1

    do i = 1,size(polyMat,1)
      do j = 1,size(polyMat,2)
        !Check orthogonality.
        if (i .ne. j) then
          do k = 1,N
            funVals(k) = orthFunTest(x(k),alpha,beta,polyMat(:,i),polyMat(:,j))
          end do
          res = trapz(N,x,funVals)
          if (res .gt. 1.0/N) then
            print *, 'othogonality check failed!', res, i,j
            pass = .false.
          end if
        end if
      end do
      !check normality.
      do k = 1,N
        funVals(k) = normFunTest(x(k),alpha,beta,polyMat(:,i))
      end do
      res = trapz(N,x,funVals)
!      print *, res
      if (abs(res-1) .gt. 1.0/N) then
        print *, 'normality check failed!', res, i
        pass = .false.
      end if
    end do
  end subroutine checkBasis

  function orthFunTest(x,alpha,beta,poly1,poly2) result(res)
    real(kind = PRECISION) :: res,x,alpha,beta
    real(kind = PRECISION), dimension(:) :: poly1, poly2

    res = omegaTest(x,alpha,beta) * polyvaL(poly1,x) * polyVal(poly2,x)
  end function orthFunTest

  function normFunTest(x,alpha,beta,poly) result(res)
    real(kind = PRECISION) :: res,x,alpha,beta
    real(kind = PRECISION), dimension(:) :: poly

    res = omegaTest(x,alpha,beta) * polyval(poly,x)**2
  end function normFunTest

  function omegaTest(x,alpha,beta) result(res)
    real(kind = PRECISION) :: res,x,alpha,beta
    if (x == -1) then
      res = 0
    else
      res = (1.0 - x)**alpha * (1.0 + x)**beta
    end if
  end function omegaTest

end program
