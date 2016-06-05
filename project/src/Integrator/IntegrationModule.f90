!A module featuring a function which is able to do numerical integration according
!to gauss' rule.
module IntegrationModule
  use UtilModule
  use FunctionModule
  implicit none

  !this abstrcat Interface allows function pointers.
  abstract interface
     function funPointer (z)
        real  :: funPointer
        real, intent (in) :: z
     end function funPointer
  end interface
contains

  !The well known trapeziodal rule.
  function trapz(N,xValues,funValues) result(res)
    integer :: N,i
    real(kind = PRECISION) :: res,h
    real(kind = PRECISION), dimension(:) :: xValues,funValues

    res = 0
    do i = 2,N
      h = (xvalues(i) - xvalues(i-1))
      res = res + h/2 * (funValues(i) + funValues(i-1))
    end do
  end function trapz

  !The trapeziodal rule with function Pointer input.
  function trapzPoint(funp,a,b,N) result(res)
      procedure(funPointer), pointer :: funp
      real :: a,b,h, res
      integer :: N, k

      h = (b - a) / N
      res = h * (funp(a) + funp(b))/2

      do k = 1,N
        res = res + h * funp(a + h*k)
      end do
  end function

  !Use quadpack to integrate.
  function quadIntegral(f, a, b) result(res)
    procedure (funPointer), intent(in), pointer :: f
    procedure (funPointer), pointer :: f2
    real :: a,b
    real :: res
    real :: ABSERR,EPSABS,EPSREL
    integer :: IER,NEVAL

    !not pretty but necessary, if f is used directly this leads to a
    !segmentation fault. (probably this comes form the legacy EXTERNAL command
    !that is used in the quadpack code.)
    f2 => f
  	!absolute and relative accuracy requested.
    EPSABS = 0.0E0
    EPSREL = 1.0E-3

    !call quadpacks universal integration routine for
    !bounded integrals.
    CALL QAGS(f2,A,B,EPSABS,EPSREL,res,ABSERR, &
              NEVAL,IER)!,LIMIT,LENW,LAST,IWORK,WORK)
  end function quadIntegral


  !Gaussian quadrature integration.
  function integral(alpha, beta, N, funValues, funp) result(sum)
    integer, intent(in) :: N !number of integration points.
    procedure(funPointer), pointer, optional :: funp
    real(kind = PRECISION), optional :: funValues(N)
    real(kind = PRECISION), dimension(N+1,2) :: ab
    real(Kind = PRECISION), dimension(N,2) :: xw
    real(kind = PRECISION) :: alpha, beta
    real(kind = PRECISION) :: sum
    integer :: i

    ab = rJacobi(N,alpha,beta)
    xw = weightGauss(N, ab)
    if (present(funValues)) then
      sum = 0
      do i= 1,N
        sum = sum + xw(i,2) * funValues(i)
      end do
    else if (present(funp)) then
      sum = 0
      do i = 1,N
        sum = sum + xw(i,2) * funp(xw(i,1))
      end do
    else
      print *, 'ERROR: Function values or a function-pointer have to be present.'
      sum = 0
    end if
  end function

  !Compute the orthogonal jacobi coefficients.
  function rJacobi(N,a,b) result(ab)
    integer, intent(in) :: N
    real(kind = PRECISION), optional :: a, b
    real(kind = PRECISION), dimension(N+1,2) :: ab
    real(kind = PRECISION) :: nu
    real(kind = PRECISION) :: mu
    real(kind = PRECISION), dimension(N+1) :: vA
    real(kind = PRECISION), dimension(N) :: nab
    real(kind = PRECISION), dimension(N-1) :: nab2
    real(kind = PRECISION),dimension(N) :: vB
    integer, dimension(N) :: sn
    integer, dimension(N-1) :: sn2
    integer :: i

    if(.not. present(a)) then
      a = 0
    end if
    if (.not. present(b)) then
      b = 0
    end if
    if ((N <= 0) .or. (a <= -1) .or. (b <= -1)) then
      print *, 'WeightModule: Error parameters out of range.'
    else
      nu = (b - a)/(a + b + 2)
      mu = 2**(a + b+1) * gamma(a + 1) * gamma(b + 1)/gamma(a+b+2)
      if (N .eq. 1) then
        ab(1,1) = nu
        ab(1,2) = mu
      else
        sn = (/ (i, i=1,N) /)
        nab = 2*sn + a + b;
        vA = 1
        vA(2:N+1) = (b**2 - a**2)*vA(2:N+1)  / (nab*(nab + 2))
        vA(1) = nu
        nab2 = nab(2:N)
        sn2 = sn(2:N)
        vB(1) = 4*(a+1)*(b+1)/((a+b+2)**2*(a+b+3));
        vB(2:N) = 4*(sn2+a)*(sn2+b)*sn2*(sn2+a+b)/((nab2**2)*(nab2+1)*(nab2-1));
        ab(:,1) = vA;
        ab(1,2) = mu;
        ab(2:N+1,2) = vB;
      end if
    end if
  end function

  !compute the gauss quadrature weights.
  function weightGauss(N,ab) result(xw)
    integer, intent(in) :: N
    real(kind = PRECISION), dimension(N+1,2) :: ab
    real(kind = PRECISION), dimension(N,N)     :: J
    real(kind = PRECISION), dimension(N) :: rD
    real(kind = PRECISION), dimension(N) :: iD
    real(Kind = PRECISION), dimension(N,N) :: V
    real(Kind = PRECISION), dimension(N)   :: vV
    real(Kind = PRECISION), dimension(N,N) :: W
    real(Kind = PRECISION), dimension(N,2) :: xw
    real(kind = PRECISION), dimension(7*N) :: work
    integer :: info
    integer :: sn

    if (size(ab) < N) then
      print *, 'GaussModule Error: Input array too short.'
    else
      J = 0
      do sn = 1,N
        J(sn,sn) = ab(sn,1)
      enddo
      do sn=2,N
        J(sn,sn-1)=sqrt(ab(sn,2))
        J(sn-1,sn)=J(sn,sn-1)
      enddo
      !Call the lapack function do replace the eig...
      !SUBROUTINE DGEEV( JOBVL, JOBVR, N, A, LDA, WR, WI, VL, LDVL, VR, LDVR, WORK, LWORK, INFO )
      call DGEEV(     'N',    'V', N, J, N,   rD, iD,  W, N,    V,   N,   work, 7*N, info)
      vV = V(1,:);
      call sort2(rD,vV)
      xw(:,1) = rD
      xw(:,2) = ab(1,2) * vV **2
    end if
  end function
end module
