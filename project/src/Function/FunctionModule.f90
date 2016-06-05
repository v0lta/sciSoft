!A module to represent functions.
module FunctionModule
  use UtilModule
  implicit none
  private

  public FunctionType
  public sortFunction
  public destroyFunction

  interface FunctionType
    module procedure newFunction
  end interface

  type FunctionType
    real(kind = PRECISION), dimension(:), allocatable  :: desValues
    real(kind = PRECISION), dimension(:), allocatable  :: functionValues
    real(kind = PRECISION), dimension(:), allocatable  :: xValues
    integer :: funSize
    integer :: desSize
  end type

contains

  !This function creats a new function pseudoobject.
  function newFunction(inputSize,functionValues,xValues,inputSize2,desValues) result(thisFunction)
    integer, intent(in) :: inputSize
    integer, intent(in), optional :: inputSize2
    real(kind = PRECISION), dimension(:), intent(in)  :: functionValues
    real(kind = PRECISION), dimension(:), intent(in), optional  :: xValues
    real(kind = PRECISION), dimension(:), intent(in), optional  :: desValues
    type(FunctionType) :: thisFunction

    !dont forget deallocation!!!!
    allocate(thisFunction%functionValues(inputSize))
    thisFunction%functionValues = functionValues
    thisFunction%funSize = size(thisFunction%functionValues)

    if (present(xValues)) then
     allocate(thisFunction%xValues(inputSize))
     thisFunction%xValues = xValues
    else
      allocate(thisFunction%xValues(0))
    end if

    if (present(desValues)) then
      allocate(thisFunction%desValues(inputSize2))
      thisFunction%desValues = desValues
      thisFunction%desSize = size(thisFunction%desValues)
    else
      allocate(thisFunction%desValues(0))
    end if

  end function

  !sort x and function values according to x values.
  subroutine sortFunction(function)
    type(functionType), intent(inout) :: function
    integer :: n
    integer :: i
    logical :: swapped = .true.

    n = function%funSize
    do while(swapped .eqv. .true.)
      swapped = .false.
      do i = 2,n
        if (function%xValues(i-1) > function%xValues(i)) then
          call swap(i-1,i,function%xValues)
          call swap(i-1,i,function%functionValues)
          swapped = .true.
        end if
      end do
      n = n - 1
    end do

    !dont forget the desired values.
    if (size(function%desValues) .ne. 0) then
      call sort(function%desValues)
    end if

  end subroutine

  !Deallocate the function variables.
  subroutine destroyFunction(function)
    type(functionType), intent(inout) :: function
    if (allocated(function%desValues) .eqv. .true.) then
      deallocate(function%desValues)
    end if
    if (allocated(function%functionValues) .eqv. .true.) then
      deallocate(function%functionValues)
    end if
    if (allocated(function%xValues) .eqv. .true.) then
      deallocate(function%xValues)
    end if
  end subroutine



end module FunctionModule
