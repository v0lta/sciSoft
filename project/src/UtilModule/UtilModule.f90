!A module which provides useful function.
!And global variables.
module UtilModule
implicit none
save
public

integer, parameter :: sp = 6
integer, parameter :: dp = 15
integer, parameter :: ep = 18
integer, parameter :: qp = 33

integer, parameter :: PRECISION = selected_int_kind(dp)
real(kind = PRECISION), parameter :: PI = 4.0*atan(1.0)

real(kind = PRECISION), public :: pAlpha, pBeta
real(kind = PRECISION), public, allocatable, dimension(:) :: pol1, pol2

contains

  subroutine printRows(mat,form,cmplxMatrix)
    character(10), optional :: form
    real(kind = PRECISION), dimension(:,:), intent(in), optional  :: mat
    complex(kind = PRECISION), dimension(:,:), intent(in), optional  :: cmplxMatrix
    integer :: i

    if (present(mat)) then
      if (present(form)) then
        do i = 1,size(mat,1)
          print *, (/ ('-',i=1,60) /), 'Row:', i
          print form, mat(i,:)
        end do
      else
        do i = 1,size(mat,1)
          print *, (/ ('-',i=1,60) /), 'Row:', i
          print *, mat(i,:)
        end do
      end if
      print *, (/ ('-',i=1,60) /), 'End Matrix'
    else if (present(cmplxMatrix)) then
      if (present(form)) then
        do i = 1,size(cmplxMatrix,1)
          print *, (/ ('-',i=1,60) /), 'Row:', i
          print form, cmplxMatrix(i,:)
        end do
      else
        do i = 1,size(cmplxMatrix,1)
          print *, (/ ('-',i=1,60) /), 'Row:', i
          print *, cmplxMatrix(i,:)
        end do
      end if
      print *, (/ ('-',i=1,60) /), 'End Matrix'
    else
      print *, ''//achar(27)//'[31m Error. print Row is missing input.'//achar(27)//'[0m'
    end if
  end subroutine

  subroutine eye(matrix, cmplxMatrix)
    integer :: rows, columns, stopDo, i
    real(kind = PRECISION), dimension(:,:), intent(inout), optional :: matrix
    complex(kind = PRECISION), dimension(:,:), intent(inout), optional :: cmplxMatrix

    if (present(matrix)) then
      rows = size(matrix,1)
      columns = size(matrix,2)

      stopDo = min(rows,columns)
      matrix = 0
      do i = 1,stopDo
        matrix(i,i) = 1
      end do
    else if (present(cmplxMatrix)) then
      rows = size(cmplxMatrix,1)
      columns = size(cmplxMatrix,2)

      stopDo = min(rows,columns)
      cmplxMatrix = 0
      do i = 1,stopDo
        cmplxMatrix(i,i) = cmplx(1.0,0)
      end do
    else
      print *, ''//achar(27)//'[31m Error. eye error. Missing input.'//achar(27)//'[0m'
    end if

  end subroutine eye

  !Flip updown of an array. (column based)
  subroutine flipud(matrix)
    integer :: rows, columns, i
    real(kind = PRECISION), dimension(:,:), intent(inout) :: matrix
    real(kind = PRECISION), dimension(:), allocatable :: tmpRow

    rows = size(matrix,1)
    columns = size(matrix,2)
    allocate(tmpRow(columns))

    do i = 1,(columns/2)
      tmpRow = matrix(i,:)
      matrix(i,:) = matrix(rows-i+1,:)
      matrix(rows-i+1,:) = tmpRow
    end do
    deallocate(tmpRow)
  end subroutine

  function flipVec(vector) result(flipped)
    complex(kind = PRECISION), dimension(:), intent(in) :: vector
    complex(kind = PRECISION), dimension(size(vector)) :: flipped
    integer :: i
    do  i = 1,size(vector)
      flipped(size(vector)-i+1) = vector(i)
    end do
  end function

  function flipTensor(tensor) result(flipped)
    complex(kind = PRECISION), intent(in), dimension(:,:,:) :: tensor
    complex(kind = PRECISION), dimension(size(tensor,1),size(tensor,2),size(tensor,3)) :: flipped
    integer :: i
    do i = 1,size(tensor,3)
      flipped(:,:,size(tensor,3)-i+1) = tensor(:,:,i)
    end do
  end function flipTensor

  !Specific repMat function for this useCase!
  function repMat(mat,n) result(bigMat)
    integer, intent(in) :: n
    complex(kind = PRECISION), dimension(1,1,n) :: mat
    complex(kind = PRECISION), dimension(2,2,n) :: bigMat
    integer :: i
    do i=1,n
      bigMat(:,:,i) = mat(1,1,i)
    end do
  end function repMat


  !Evaluate polynomials like matlab does (one one dimesion...)
  function polyval(poly,x) result(res)
    integer :: i,l
    real(kind = PRECISION) :: res,x
    real(kind = PRECISION), dimension(:), intent(in) :: poly

    l = size(poly)
    res = 0
    do i = 1,(l-1)
       res = res + poly(i) * x**(l-i)
    end do
    res = res + poly(l)
  end function polyval

  !Convolution for polinomial multiplication.
  function conv(p1,p2,s1,s2) result(res)
    integer, intent(in) :: s1,s2
    integer :: i
    real(kind = PRECISION), dimension(s1,s1+s2-1) :: conMat
    real(kind = PRECISION), dimension(:) :: p1,p2
    real(kind = PRECISION), dimension(s1+s2-1) :: res

    conMat = 0
    do i = 1,s1
      !Fill the rows of the convolution Matrix.
      conMat(i,i:(i+s2-1)) = p1(i) * p2
    end do

    do i = 1,(s1+s2-1)
      !Sum up the columns to obtain the desired result.
      res(i) = sum(conMat(:,i))
    end do
  end function

  ! A subroutine to generate a linear spacing in an array from the min value
  ! till the max value (both end included) (NOT TESTED)
  subroutine linspace(minArg,maxArg,array)
    real(kind=PRECISION),dimension(:)::array
    real(kind=PRECISION) :: minArg,maxArg,dx
    integer :: i,sizeOfArray

    sizeOfArray = size(array)
    dx = (maxArg-minArg)/(sizeOfArray-1)
    array = (/ ((i*dx+minArg),i=0,(sizeOfArray-1)) /)
  end subroutine linspace

  !Bubblesort for now to be replaced later with something faster.
  subroutine sort(array)
    real(kind = PRECISION), dimension(:), intent(inout) :: array
    integer :: n
    integer :: i
    logical :: swapped = .true.

    n = size(array)
    do while(swapped .eqv. .true.)
      swapped = .false.
      do i = 2,n
        if (array(i-1) > array(i)) then
          call swap(i-1,i,array)
          swapped = .true.
        end if
      end do
      n = n - 1
    end do
  end subroutine

  subroutine sort2(array,carryAlong)
    real(kind = PRECISION), dimension(:), intent(inout) :: array
    real(kind = PRECISION), dimension(:), intent(inout) :: carryAlong
    integer :: n
    integer :: i
    logical :: swapped = .true.

    n = size(array)
    do while(swapped .eqv. .true.)
      swapped = .false.
      do i = 2,n
        if (array(i-1) > array(i)) then
          call swap(i-1,i,array)
          call swap(i-1,i,carryAlong)
          swapped = .true.
        end if
      end do
      n = n - 1
    end do
  end subroutine

  !sort and apply the same changes to another array. (NOT TESTED)
  subroutine swap(pos1,pos2,array)
    integer :: pos1,pos2
    real(kind = PRECISION), dimension(:), intent(inout) :: array
    real(kind = PRECISION) :: tmp

    tmp = array(pos1)
    array(pos1) = array(pos2)
    array(pos2)= tmp
  end subroutine
end module UtilModule
