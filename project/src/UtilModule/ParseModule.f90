module ParseModule
  implicit none


contains

!test, lagrange, time ,integration, gauss, pointNoSet, weightSet,
!degreeSet, doOrthogonal, isProj
subroutine readCmdLine(hlp, tst, lag, tme, int, intGaus, quad, pNoS, ws, degS, orth, &
                proj, wAlpha, wBeta, gAlpha, gBeta, N, M, evOrth, h, x, &
                isGM, isFav, isAsy, T)
  integer :: loop
  logical, intent(inout)  :: hlp, tst, lag, tme, int, intGaus, quad, pNoS, ws, &
                             degS, orth, proj, evOrth, isGM, isFav, isAsy
  real, intent(inout)  :: wAlpha, wBeta, gAlpha, gBeta, x
  integer, intent(inout) :: N, M, T , h
  character(len=32) :: arg
  !Set debug to false in order to supress additional output.
  logical, parameter :: debug = .false.
  !Set all the control variables to .false.
  call setToFalse(hlp, tst, lag, tme, int, intGaus, quad, pNoS, ws, degS, orth, &
                  proj, evOrth, isGM, isFav, isAsy)

  !Go trough all the command line arguments supplied and set control variables
  !accordingly.
  loop = 0
  do
    call get_command_argument(loop, arg)
    if (debug .eqv. .true.) then
      print *, "current arg: ", arg
    end if
    if ( (len_trim(arg) == 0) .and.  (loop .eq. 1)) then
        print *, "No valid command could be parsed."
        exit
    else if (len_trim(arg) == 0) then
       if (debug .eqv. .true.) then
        print *, "exiting at:" ,loop
       end if
       exit
    else
        if (arg .eq. "--help") then
          hlp = .true.
        else if (arg .eq. "--test") then
          tst = .true.
        else if (arg .eq. "--lagrange") then
          lag = .true.
        else if (arg .eq. "--int") then
          int = .true.
        else if (arg .eq. "--intMeth") then
          call get_command_argument(loop+1, arg)
          if (arg .eq. "gauss") then
            intGaus = .true.
            call get_command_argument(loop+2, arg)
            read (arg,'(f8.0)') gAlpha
            !print *, alpha
            call get_command_argument(loop+3, arg)
            read (arg,'(f8.0)') gBeta
            !print *, beta
            !read (arg,'(f8.0)') beta
          else if (arg .eq. "quadpack") then
            quad = .true.
          end if
        else if (arg .eq. "-t") then
          !print *, 'found tme'
          tme = .true.
        else if (arg .eq. "-i") then
          pNoS = .true.
          !Assumig the next argument is the point number:
          call get_command_argument(loop+1, arg)
          read (arg,'(i5)') N
          !no need to look at this argument again:
        else if(arg .eq. "--algo") then
          orth = .true.
          !print *, orth
          call get_command_argument(loop+1, arg)
          if (arg .eq. "proj") then
            proj = .true.
          else
            proj = .false.
          end if
        else if(arg .eq. "-w") then
          ws = .true.
          call get_command_argument(loop+1, arg)
          read (arg,'(f8.0)') wAlpha
          call get_command_argument(loop+2, arg)
          read (arg,'(f8.0)') wBeta
        else if (arg .eq. "-n") then
          call get_command_argument(loop+1, arg)
          read(arg,'(i4)') M
        else if (arg .eq. "--evalOrth") then
          evOrth = .true.
          !Read in h
          call get_command_argument(loop+1, arg)
          read(arg,'(i4)') h
          !read in x.
          call get_command_argument(loop+2, arg)
          read (arg,'(f8.0)') x
        else if (arg .eq. "--orthMeth") then
          call get_command_argument(loop+1, arg)
          if (arg .eq. "gs") then
            !do gram schmidt.
            isGM = .true.
          else if(arg .eq. "recur") then
            !do recursion.
            isFav = .true.
          else if (arg .eq. "asy") then
            !do asy.
            isAsy = .true.
            call get_command_argument(loop+3, arg)
            read (arg,'(i5)') T
          end if

       end if


    end if
    !write (*,*)trim(arg)
    loop = loop +1
  end do

  if (debug .eqv. .true.) then
  print *, "hlp, tst, lag, tme, int, intGaus, pNoS, ws, degS, orth, proj -- part1"
  print *, hlp, tst, lag, tme, int, intGaus, pNoS, ws, degS, orth, proj
  print *, "gAlpha, gBeta, wAlpha, wBeta, M, N:"
  print *, gAlpha, gBeta, wAlpha, wBeta, M, N
  print *, "evOrth, isGM, isFav, isAsy, quad -- part2"
  print *, evOrth, isGM, isFav, isAsy, quad
  print *, "h,x,T:"
  print *, h,x,T
  print *, '-----------------------parsing done---------------------------------'
  end if
end subroutine ReadCmdLine


  subroutine printHelp()
    print *, "The basis executable my be used with the following arguments:"
    print *, './basis'//achar(27)//'[1m --help  '//achar(27)//'[0m', "prints this ", &
               "help summary."
    print *, './basis'//achar(27)//'[1m --test  '//achar(27)//'[0m', "Runs a set of ", &
              "predefined performance tests."
    print *, './basis'//achar(27)//'[1m --lagrange  '//achar(27)//'[0m', "Approximation ",&
              "with plain Lagrange polynomials."
    print *,  '       '//achar(27)//'[1m -t  '//achar(27)//'[0m times this computation.'
    print *, './basis'//achar(27)//'[1m --int  '//achar(27)//'[0m', "Integration ",&
              "may be followd by  --intMeth"
    print *, '       '//achar(27)//'[1m --intMeth  '//achar(27)//'[0m', 'Integration method.', &
                      ' With'//achar(27)//'[1m gauss{alpha} {beta} '//achar(27)//'[0m', 'or', &
                      ''//achar(27)//'[1m trapz  '//achar(27)//'[0m.'
    print *,  '       '//achar(27)//'[1m -t  '//achar(27)//'[0m times this computation.'
    print *, './basis'//achar(27)//'[1m --algo  '//achar(27)//'[0m', "Orthogonal", &
                      "decomposition with", &
                      ' '//achar(27)//'[1m proj  '//achar(27)//'[0m or', &
                      ' '//achar(27)//'[1m interpol  '//achar(27)//'[0m.'
    print *, '       '//achar(27)//'[1m --intMeth  '//achar(27)//'[0m', 'Integration method.', &
                     ' With'//achar(27)//'[1m gauss{alpha} {beta} '//achar(27)//'[0m', ',', &
                     ''//achar(27)//'[1m trapz  '//achar(27)//'[0m or'
    print *, '                  '//achar(27)//'[1m quadpack  '//achar(27)//'[0m' 
    print *, '       '//achar(27)//'[1m -w {alpha} {beta} '//achar(27)//'[0m sets ', &
                      'the jacobi function weights.'
    print *,  '       '//achar(27)//'[1m -n {degree} '//achar(27)//'[0m sets the ', &
                      'degree of the basis.'
    print *,  '       '//achar(27)//'[1m -i {#points} '//achar(27)//'[0m sets the ', &
                      'number of integraton points.'
    print *,  '       '//achar(27)//'[1m -t  '//achar(27)//'[0m times this computation.'
    print *,  './basis'//achar(27)//'[1m --evalOrth{h}{x} '//achar(27)//'[0m get value ', &
                      'of degree h orthogonal polynomial @ x.'
    print *,  './basis'//achar(27)//'[1m --orthMeth{meth}({T})'//achar(27)//'[0m ', &
              'basis computation method:', &
              ''//achar(27)//'[1m gs, recur, asy{reg}{T}'//achar(27)//'[0m.'
  end subroutine printHelp

  !Make shure all logicals are .false.!
  subroutine setToFalse(hlp, tst, lag, tme, int, intMeth, quad, pNoS, &
                        ws, orth, deg, proj, evOrth, isGM, isFav, isAsy)
    logical, intent(inout) :: hlp,  tst, lag, tme, int, intMeth, quad, pNoS, &
                        ws, orth, deg, proj, evOrth, isGM, isFav, isAsy
    hlp = .false.; tst = .false. ; lag = .false. ; tme = .false. ; int = .false.;
    intMeth = .false.; quad = .false.; pNoS = .false. ; ws = .false. ; deg = .false. ;
    proj = .false. ; orth = .false.; evOrth = .false.; isGM = .false.;
    isFav = .false.; isAsy = .false.;
  end subroutine setToFalse

end module ParseModule
