module LagrangeModule
  use FunctionModule
  use UtilModule
  implicit none

contains

  !use the function data and the desired points to interpolate
  function lagrange(givenFunction) result(interpValues)
    type(FunctionType) :: givenfunction
    integer :: i
    integer :: j
    integer :: k
    real(kind = PRECISION), dimension(givenFunction%funSize)  :: l
    real(kind = PRECISION) :: tmp
    real(kind = PRECISION), dimension(givenFunction%desSize)  :: interpValues
    real(kind = PRECISION), dimension(givenFunction%funSize)  :: lambda

    lambda = 0
    tmp = 0
    !compute the lambdas.
    do i =1,givenFunction%funSize
      tmp = 1
      do j =1,givenFunction%funSize
        if (i .ne. j) then
          tmp =  tmp * (givenFunction%xValues(i) - givenFunction%xValues(j))
        end if
      end do
      lambda(i) = 1/tmp
    end do

    l = 0
    interpValues = 0
    !compute l.
    do k= 1,givenFunction%desSize
      do i = 1,givenFunction%funSize
        tmp = 0
        do j = 1,givenFunction%funSize
          !compute the sum of different mu's.
          tmp = tmp + lambda(j)/(givenFunction%desValues(k) - givenFunction%xValues(j))
        end do
      l(i) = (lambda(i)/(givenFunction%desValues(k) - givenFunction%xValues(i)))/tmp
      end do
    interpValues(k) =  dot_product(givenFunction%functionValues, l)
    end do
  end function

  !Generate lagrange polynomials up to a desired degree and return the result in
  !the matrix polyMat.
  subroutine genLagrangePols(deg,polyMat)
    integer, intent(in) :: deg
    real(kind = PRECISION), dimension (:,:), intent(inout) :: polyMat
    real(kind = PRECISION), dimension(deg)  :: x
    real(kind = PRECISION) :: denom
    real(kind = PRECISION), dimension(deg)  :: nom
    real(kind = PRECISION), dimension(deg+1)  :: res
    integer :: i,j

    call linspace(-1.0,1.0,x)
    !print *, x
    !print *, deg
    do i = 1,deg
      denom = 1
      nom = 0
      nom(deg) = 1
      do j = 1,deg

        if (i .ne. j) then
            denom = denom * 1/(x(i) - x(j))
            res = conv(nom, (/ 1.0, - x(j) /),deg,2)
            !print *, 'denom:', denom
            !print *, 'convResult', conv(nom, (/ 1.0, - x(j) /),deg,2)
            !print *, 'res:', res
            nom = res(2:deg+1)
        end if
      end do
      polyMat(i,:) = nom * denom
    end do
  end subroutine

end module
