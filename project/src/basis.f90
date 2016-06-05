program main
  use ReadModule
  use IntegrationModule
  use UtilModule
  use FunctionModule
  use LagrangeModule
  use GramSchmidtModule
  use ProjectionModule
  use ParseModule
  use FavardModule
  use AsyModule
  implicit none

  integer :: i, intPoints, M, intMethod, basisMethod
  logical  :: hlp, tst, lag, tme, int, intGaus, quad, pNoS, ws, degS, orth, proj
  real(kind = PRECISION)  :: wAlpha, wBeta, gAlpha, gBeta, sum
  integer, parameter :: MAX_PATH_LEN = 1024, CHROMOSOME_UNIT = 10
  !part 2
  logical :: evOrth, isGM, isFav, isAsy
  integer :: T, h
  real(kind = PRECISION) :: x

  !Set defaults.
  gAlpha = 0
  gBeta = 0
  wAlpha = 0
  wBeta = 0
  M = 6
  intPoints = 500
  intMethod = 0
  basisMethod = 0

  call readCmdLine(hlp, tst, lag, tme, int, intGaus, quad, pNoS, ws, degS, orth, &
                  proj, wAlpha, wBeta, gAlpha, gBeta, intPoints, M, evOrth, h, x, &
                  isGM, isFav, isAsy, T)

  !----------------------------Parsing done. Begin execution.---------------------
  !figure out which methods to use.
  if (isFav .eqv. .true.) then
    basisMethod = 1
  else if (isAsy .eqv. .true.) then
    basisMethod = 2
  else
    !Default method Gram-Schmidt with lagrange polys.
    basisMethod = 0
  end if
  if (intGaus .eqv. .true.) then
    intMethod = 1
  else if (quad .eqv. .true.) then
    intMethod = 2
  else
    !Default integration use the trapezoidal rule.
    intMethod = 0
  end if

  !are we going to test things...
  if (hlp .eqv. .true.) then
    call printHelp()
   else if (tst .eqv. .true.) then
     call runTests()
   !Check if integration is desired.
   else if (int .eqv. .true.) then
     if (intGaus .eqv. .true.) then
       call executeIntegration(gAlpha,gBeta,tme,sum,1)
     else
       call executeIntegration(0.0,0.0,tme,sum,0)
     end if
   !Lagrange interpolation?
   else if (lag .eqv. .true.) then
      call plainLagrangeInterpolation(tme)
   !Take care of the orthogonal things. Figure out the integration method to use.
   else if (orth .eqv. .true.) then
    call orthogonalInterpolation(M,intPoints,gAlpha,gBeta,wAlpha,wBeta,intMethod, &
                                basisMethod,T, proj,tme)
    !possibly the polynomial is to be evaluated.
    else if (evOrth .eqv. .true.) then
      call evalOrth(basisMethod, intMethod, M, T, galpha, gbeta, walpha, wbeta,&
                            intPoints, h,x)
    !The default operation, when no input Arguments could be parsed.
    else if ((hlp .eqv. .false.) .and. (tst .eqv. .false.) .and. (int .eqv. .false.) &
          .and. (orth .eqv. .false.)) then
      call orthogonalInterpolation(M,intPoints,gAlpha,gBeta,wAlpha,wBeta,1, &
                                      basisMethod,T, proj,tme)
  end if

contains

  subroutine evalOrth(basisMethod, intMethod, M, T, galpha, gbeta, walpha, &
                      wbeta, intPoints, h,x)
    integer, intent(in) :: basisMethod, intMethod, M,T, intPoints, h
    real(kind = PRECISION), dimension(:,:), allocatable :: polyMat2
    real(kind = PRECISION), intent(in) :: galpha, gbeta, walpha, wbeta, x

    !Get the basis.
    call getBasis(basisMethod, intMethod, polyMat2, M,T, galpha, gbeta, walpha, wbeta, &
                          intPoints)
    !Evaluate it.
    !call printRows(polyMat2)
    print *, polyval(polyMat2(:,h), x)
    deallocate(polyMat2)
  end subroutine evalOrth


  subroutine getBasis(basisMethod, intMethod, polyMat, M,T, galpha, gbeta, walpha, wbeta,&
                        intPoints)
    integer, intent(in) :: basisMethod, M,T, intPoints, intMethod
    real(kind = PRECISION), dimension(:,:), allocatable :: polyMat
    real(kind = PRECISION), intent(in) :: galpha, gbeta, walpha, wbeta
    !Asy variables
    integer :: j
    real(kind = PRECISION) :: Dinf
    complex(kind = PRECISION) :: ipiVal
    complex(kind = PRECISION), dimension(4) :: asyC, asyD
    complex(kind = PRECISION), dimension(2,2,T,ceiling((T-1.0)/2.0)) :: Uleft
    complex(kind = PRECISION), dimension(2,2,T,ceiling((T-1.0)/2.0)) :: Uright
    real(kind = PRECISION), dimension(M) :: z

    !Get the basis.
    if (basisMethod .eq. 1) then
      !use the recursive Method
      allocate(polyMat(M,M))
      polyMat = 1
      call favard(M, walpha, wbeta, polyMat)
    else if (basisMethod .eq. 2) then
      !use asymtotic expansions.
      allocate(polyMat(T,M))
      polyMat = 1
      Dinf = 2**(-wAlpha/2.0-wBeta/2.0);
      Uright = 0; Uleft = 0;
      asyD = 0; asyC = 0;
      call compUExplicit(wAlpha, wBeta, Dinf, asyC , asyD, Uright, Uleft)

      !pick a point somewhere in the lens.
      call linspace(-0.5,0.5,z)
      Dinf = 2**(-wAlpha/2.0-wBeta/2.0);
      do i = 1,M
            do j = 1,T
            call asyLens(i**2,j,cmplx(z(i),0.0),wAlpha,wBeta,Dinf,Uright,Uleft,ipiVal)
            !negecting the complex part is ok. It should be very small or zero!
            polyMat(j,i) = real(ipiVal)
            end do
      end do
      !polyMat = transpose(polyMat)
      !call printRows(polyMat)
    else
      !Default use Gram-Schmidt orthogonalized Lagrange polynomials.
      allocate(polyMat(M,M))
      call genLagrangePols(M,polyMat)
      !orthogonalize the basis.
      call gramSchmidt(intPoints,gAlpha,gBeta,wAlpha,wBeta,polyMat,intMethod)
    end if
  end subroutine getBasis

  subroutine orthogonalInterpolation(M,intPoints,gAlpha,gBeta,wAlpha,wBeta,intMethod, &
                                  basisMethod, T,isProj,time, &
                                  CHROMOSOME_UNIT, resultOut)
    integer, intent(in) :: M, T, intPoints, intMethod, basisMethod
    integer, intent(in), optional :: CHROMOSOME_UNIT
    real(kind = PRECISION), optional, intent(out), dimension(:) :: resultOut
    logical, intent(in) :: isProj, time
    real(kind = PRECISION) :: gAlpha, gBeta, wAlpha, wBeta
    real(kind = PRECISION), dimension(:,:), allocatable :: polyMat
    type(FunctionType)  :: f
    real(kind = PRECISION) :: transA,transB
    real(kind = PRECISION),dimension(M) :: c
    real(kind = PRECISION), dimension(:), allocatable :: result
    real(kind = PRECISION), dimension(:,:),allocatable :: evalMat
    real :: start,finish

    if (time .eqv. .true.) then
      call cpu_time(start)
    end if

    !Read in data, sort it and allocate variables accordingly.
    if (present(CHROMOSOME_UNIT)) then
      f = readFunctionData(CHROMOSOME_UNIT)
    else
      f = readFunctionData()
    end if

    call sortFunction(f)
    allocate(evalMat(f%funSize,M))
    allocate(result(f%desSize))

    !Get the basis.
    call getBasis(basisMethod, intMethod, polyMat, M,T, galpha, gbeta, walpha, wbeta, &
                          intPoints)

    !Transform x to -1,1
    call transformX(f%xValues,f%desValues, transA, transB)
    !Evaluate the polynomials in x
    call evalMatrix(polyMat, f%xValues, M, evalMat)
    !compute a orthogonal function representation.
    if (isProj .eqv. .false.) then
      c = compC(evalMat,f,f%funSize, M)
    else
      c = projC(gAlpha,gBeta,wAlpha,wBeta,evalMat,f,M,intMethod)
    end if
    !Ajust size according to the desired values.
    deallocate(evalMat)
    allocate(evalMat(f%desSize,M))
    !evaluate the polynomials at the desired points
    call evalMatrix(polyMat, f%desValues, M, evalMat)
    !compute the result using the representation found earlier.
    result = matmul(evalMat,c)

    !check if subroutine output is desired or wether the result should go
    !onto the command line.
    if(present(resultOut)) then
      resultOut = result
    else
      print '(f16.4)', result
    end if

    if (time .eqv. .true.) then
      call cpu_time(finish)
      print *, 'computation time:' , finish-start
    end if
    !done...deallocate everything
    call destroyFunction(f)
    deallocate(result)
    deallocate(evalMat)
    deallocate(polyMat)
  end subroutine orthogonalInterpolation

  !Call the non-orthogonal Lagrange interpolation routines.
  subroutine plainLagrangeInterpolation(tme)
    type(FunctionType) :: f
    logical, intent(in) :: tme
    real :: start,finish = 0
    !Start the clock.
    call cpu_time(start)
    f = readFunctionData()
    call sortFunction(f)
    !Du the computation and put the result on the standart output.
    print '(es24.16)', lagrange(f)
    call cpu_time(finish)
    if (tme .eqv. .true.) then
      !Print the time it took to do the computations.
      print *, 'computation time:' , finish-start
    end if
    call destroyFunction(f)
  end subroutine plainLagrangeInterpolation

  !Call the integration related subroutines.
  subroutine executeIntegration(gAlpha,gBeta,time,sum,method,CHROMOSOME_UNIT)
    logical, intent(in) :: time
    real(kind = PRECISION), intent(in) :: gAlpha,gBeta
    integer, intent(in) :: method
    integer, intent(in), optional :: CHROMOSOME_UNIT
    real :: start,finish = 0
    real(kind = PRECISION), intent(out) :: sum
    type(FunctionType) :: f
    real(kind = PRECISION), allocatable, dimension(:) :: x

    if (present(CHROMOSOME_UNIT)) then
      !Read from file.
      f = readFunctionData(CHROMOSOME_UNIT)
    else
      !Read from standart input.
      f = readFunctionData()
      call sortFunction(f)
    end if
    !start the timer.
    call cpu_time(start)
    !See which integration method is desired.
    if (method .eq. 1) then
      sum = integral(gAlpha,gBeta,f%funSize, f%functionValues )
    else
      allocate(x(f%funSize))
      sum = trapz(f%funSize, f%xValues, f%functionValues)
      deallocate(x)
    end if
    call cpu_time(finish)
    if (time .eqv. .true.) then
      print *, 'computation time:' , finish-start
    end if
    print '(f16.12)', sum
    call destroyFunction(f)
  end subroutine executeIntegration

!------------------------------tests from part 1 --------------------------------
  !This subroutine goes trought the predefined tests.
  subroutine runTests()
    call integrationTest()
    call expInterpolationTest()
    call ref3Test()
    call eliaTest()
    call abiTest()
    print *, (/ ('-',i=1,60)  /), 'done testing.'
  end subroutine runTests

  !Predict the Annheuser-Bush stock values.
  subroutine abiTest()
    !File i/o variables.
    character(MAX_PATH_LEN) :: filename
    integer :: open_status, close_status,i
    real(kind = PRECISION), dimension(192) :: resultOut
    type(FunctionType) :: solFunction


    print *, (/ ('-',i=1,60)  /), 'Abi test.'

    !open the proper data file.
    filename = './tests/ABI/abi.in'
    open (unit=CHROMOSOME_UNIT, file=filename, status='old', &
          iostat=open_status, action='read', position='rewind')

    !Somewhat arbitraty weights I know...but they work.
    walpha = 0.1
    wbeta = 1.5
    galpha = 0
    gbeta = -0.4

    call orthogonalInterpolation(4,intPoints,gAlpha,gBeta,wAlpha,wBeta, &
                          0,0,0,.false.,.true., CHROMOSOME_UNIT,resultOut)

    !close the current file.
    close(CHROMOSOME_UNIT, iostat=close_status)

    !Read the solution.
    filename = './tests/elia/matlab/matlabAbi.test'
    open (unit=CHROMOSOME_UNIT, file=filename, status='old', &
          iostat=open_status, action='read', position='rewind')

    solFunction = readFunctionData(CHROMOSOME_UNIT)
    call sortFunction(solFunction)

    print *, 'Two norm distance to matlab approx: ',  &
          norm2(abs(resultOut-solFunction%functionValues))
    close(CHROMOSOME_UNIT, iostat=close_status)
    !deallocate the function variables.
    call destroyFunction(solFunction)
  end subroutine abiTest

  !Perform the eila power test.
  subroutine eliaTest()
    !File i/o variables.
    character(MAX_PATH_LEN) :: filename
    integer :: open_status, close_status,i
    real(kind = PRECISION), dimension(192) :: resultOut
    type(FunctionType) :: solFunction

    print *, (/ ('-',i=1,60)  /), 'Elia test.'

    !open the first reference file.
    filename = './tests/elia/elia1115.in'
    open (unit=CHROMOSOME_UNIT, file=filename, status='old', &
          iostat=open_status, action='read', position='rewind')

    !Somewhat arbitraty weights I know...but they work.
    walpha = 0.1
    wbeta = 1.5
    galpha = 0
    gbeta = -0.4

    call orthogonalInterpolation(4,intPoints,gAlpha,gBeta,wAlpha,wBeta, &
                          0,0,0,.false.,.true., CHROMOSOME_UNIT,resultOut)
    close(CHROMOSOME_UNIT, iostat=close_status)

    !Read the solution.
    filename = './tests/elia/matlab/matlabElia.test'
    open (unit=CHROMOSOME_UNIT, file=filename, status='old', &
          iostat=open_status, action='read', position='rewind')

    solFunction = readFunctionData(CHROMOSOME_UNIT)
    call sortFunction(solFunction)

    print *, 'Two norm distance to matlab approx: ',  &
          norm2(abs(resultOut-solFunction%functionValues))
    print *, ''//achar(27)//'[37m This is not as bad as it looks. Keep in mind&
     &we are talking about 192 values with magnitude 1.0e+06 here. Please refer&
     &to the report for a plot '//achar(27)//'[0m.'
    close(CHROMOSOME_UNIT, iostat=close_status)
    call destroyFunction(solFunction)
  end subroutine eliaTest

  subroutine integrationTest()
    !File i/o variables.
    character(MAX_PATH_LEN) :: filename
    integer :: open_status, close_status,i
    !ref1.in variables
    real(kind = PRECISION) :: sol = 4.993379538113
    real(kind = PRECISION) :: sum

    !---Do the integration test-----------------------------------------------------
        print *, (/ ('-',i=1,60)  /), 'first test. ref1.in'
        !open the first reference file.
        filename = './tests/ref1.in'
        open (unit=CHROMOSOME_UNIT, file=filename, status='old', &
              iostat=open_status, action='read', position='rewind')
        !Do integration on it interpolation
        call executeIntegration(0.4,-0.13,.true.,sum,1,CHROMOSOME_UNIT)
        write (*,'("  ",f15.12)') sol
        if (abs(sum -sol) < 0.00000000001) then
          print *, ''//achar(27)//'[92m ref1.in passed'//achar(27)//'[0m.'
        else
            print *, 'ref1.in failed!.'
        end if
        close(CHROMOSOME_UNIT, iostat=close_status)
  end subroutine integrationTest

  subroutine expInterpolationTest()
    !Exp interpolation variables
    real(kind = PRECISION) :: sol
    real(kind = PRECISION), dimension(3) :: desValues
    real(kind = PRECISION), dimension(26) :: functionValues
    real(kind = PRECISION), dimension(26) :: xValues
    real(kind = PRECISION), dimension(1) :: lagApprox
    real :: start, finish
    !interpol variables
    integer, parameter :: M=21
    real(kind = PRECISION), dimension(M,M) :: polyMat
    real(kind = PRECISION) :: transA,transB
    real(kind = PRECISION), dimension(M) :: c
    real(kind = PRECISION), dimension(:), allocatable :: result
    real(kind = PRECISION), dimension(:,:),allocatable :: evalMat
    type(FunctionType) :: testFunction

    !---interpolate function values of the exponential function. (lagrange)---------
        print *, (/ ('-',i=1,50)  /), 'second test Exp() lagrange'
        call linspace(-2*PI,2*PI,xValues)
        functionValues = exp(xValues)
        desValues = (/ -PI/16, 0.0 , PI/16 /);
        testFunction = FunctionType(size(functionValues),functionValues,xValues, &
                          size(desValues),desValues)
        call cpu_time(start)
        lagApprox = lagrange(testFunction)
        call cpu_time(finish)
        print *,  'computation time:', finish-start
        print *, lagrange(testFunction)
        print *, exp(desValues)
        sol = exp(desValues(1))
        if (abs(lagApprox(1) -sol) < 0.00000000001) then
          print *, ''//achar(27)//'[92m exp(x)-test lagrange passed'//achar(27)//'[0m.'
        else
            print *, 'exp(x)-test lagrange failed!.'
        end if

        print *, (/ ('- - ',i=1,11)  /), 'second test Exp() interpol trapz'
        !interpolate a function value of the exponential function. (trapz-interpol)
        call cpu_time(start)
        allocate(evalMat(testFunction%funSize,M))
        allocate(result(testFunction%desSize))
        !call genLagrangePols(M,polyMat)
        call eye(polyMat)
        call flipud(polyMat)
        !Gram-schmidt method = 0 => trapz is used.
        call gramSchmidt(600,0.0,0.0,0.0,0.0,polyMat,0)
        call transformX(testFunction%xValues,testFunction%desValues, transA, transB)

        call evalMatrix(polyMat, testFunction%xValues, M, evalMat)
        c = compC(evalMat,testFunction,testFunction%funSize, M)
        !Ajust size according to the desired values.
        deallocate(evalMat)
        allocate(evalMat(testFunction%desSize,M))
        !evaluate the polynomials at the desired points
        call evalMatrix(polyMat, testFunction%desValues, M, evalMat)
        !compute the result using the representation found earlier.
        result = matmul(evalMat,c)


        call cpu_time(finish)
        print *,  'computation time:', finish-start
        print *, result
        print *, 'two norm distance to exact solution:', norm2(abs(exp(desValues)-result))
        deallocate(result)
        deallocate(evalMat)
        call destroyFunction(testFunction)
  end subroutine expInterpolationTest

  subroutine ref3Test()
    character(MAX_PATH_LEN) :: filename
    integer :: open_status, close_status,i
    real(kind = PRECISION), dimension(:), allocatable :: result
    type(FunctionType) :: solFunction

    print *, (/ ('-',i=1,60)  /), 'Third test. ref3.in'
    !Read the solution.
    filename = './tests/ref3.test'
    open (unit=CHROMOSOME_UNIT, file=filename, status='old', &
          iostat=open_status, action='read', position='rewind')

    solFunction = readFunctionData(CHROMOSOME_UNIT)
    call sortFunction(solFunction)

    !print *, solFunction%functionValues
    close(CHROMOSOME_UNIT, iostat=close_status)

    !open the third reference file.
    filename = './tests/ref3.in'
    open (unit=CHROMOSOME_UNIT, file=filename, status='old', &
          iostat=open_status, action='read', position='rewind')

    !compute the interpolation values.
    allocate(result(30))

    call orthogonalInterpolation(6,intPoints,0.4,-0.13,0.4,-0.13, &
                                  0,0,0, .false., .true., CHROMOSOME_UNIT,result)

    close(CHROMOSOME_UNIT, iostat=close_status)
    print *,'Interpolation-error two norm: ', norm2(abs(result-solFunction%functionValues))
    deallocate(result)
    call destroyFunction(solFunction)
  end subroutine ref3Test

end program
