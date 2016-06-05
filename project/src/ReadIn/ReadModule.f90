!A module to read in function files.
!Amd retirms a function object
module ReadModule
  use FunctionModule
  use UtilModule
  implicit none
  private

  public readFunctionData

contains

  !Read in input with function data and return a new function pseudoobject.
  function readFunctionData(CHROMOSOME_UNIT) result(returnFunction)
    integer, optional, intent(in) :: CHROMOSOME_UNIT
    type(FunctionType)  :: returnFunction
    real(kind=PRECISION), dimension(:),allocatable :: desValues, functionValues, xValues
    integer :: inputSize = 0
    integer :: inputSize2 = 0
    integer :: iostat = 0

    if( present(CHROMOSOME_UNIT) ) then
      !Read from file.
      read (CHROMOSOME_UNIT, *, iostat=iostat) inputSize
      allocate(functionValues(inputSize))
      allocate(xValues(inputSize))

      functionValues = 0
      xValues = 0
      read (CHROMOSOME_UNIT, *, iostat=iostat) functionValues
      read (CHROMOSOME_UNIT, *, iostat=iostat) xValues

      !Check if there is more data to read.
      if (iostat .ne. 0) then
        returnFunction = FunctionType(inputSize, functionValues)
      else
        read (CHROMOSOME_UNIT, *, iostat=iostat) inputSize2
        allocate(desValues(inputSize2))
        desValues = 0
        read (CHROMOSOME_UNIT, *, iostat=iostat) desValues
        returnFunction = FunctionType(inputSize,functionValues,xValues,inputSize2,desValues)
      end if
    else
      !Read from standart input
      read *, inputSize
      allocate(functionValues(inputSize))
      allocate(xValues(inputSize))

      functionValues = 0
      xValues = 0
      read (*,*, iostat=iostat) functionValues
      read (*,*, iostat=iostat) xValues

      !Check if there is more data to read.
      if (iostat .ne. 0) then
        returnFunction = FunctionType(inputSize, functionValues)
      else
        read (*,*, iostat=iostat) inputSize2
        allocate(desValues(inputSize2))
        desValues = 0
        read (*,*, iostat=iostat) desValues
        returnFunction = FunctionType(inputSize,functionValues,xValues,inputSize2,desValues)
      end if
    end if

  end function readFunctionData
end module ReadModule
