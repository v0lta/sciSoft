module projectionModule
  use UtilModule
  use LagrangeModule
  use GramSchmidtModule

contains

  function projC(gAlpha,gBeta,wAlpha,wBeta,polyMat,f, M, method) result(c)
    integer, intent(in) :: method, M
    integer :: i,j
    real(kind = PRECISION) :: gAlpha, gBeta, wAlpha, wBeta
    real(kind = PRECISION), dimension(:), allocatable :: fVals
    real(kind = PRECISION), dimension(:,:) :: polyMat
    real(kind = PRECISION), dimension(M) :: c
    type(functionType) :: f

    allocate(fVals(f%funSize))
    c = 0
    do j = 1,M
      do i = 1,f%funSize
        fVals(i) = omega(f%xValues(i), wAlpha,wBeta) * polyMat(i,j) &
        * f%functionValues(i)
      end do
      if (method .eq. 1) then
        c(j) = integral(gAlpha, gBeta, f%funSize, fVals)
      else
        c(j) = trapz(f%funSize, f%xValues, fVals)
      end if
    end do
    deallocate(fVals)
  end function



end module
