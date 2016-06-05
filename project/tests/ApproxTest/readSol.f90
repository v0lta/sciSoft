!run intepolations of the ref3.in file and compare the results to ref3.text.
program approxTestRef3
  use LagrangeModule
  use UtilModule
  use FunctionModule
  use ReadModule
  implicit none

  type(FunctionType)  :: solFunction

  solFunction = readFunctionData()
  call sortFunction(solFunction)

  print '(es24.16)', solFunction%xValues





end program
