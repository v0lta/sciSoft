!run intepolations of the ref3.in file and compare the results to ref3.text.
program approxTestRef3
  use LagrangeModule
  use UtilModule
  use FunctionModule
  use ReadModule

  type(FunctionType)  :: testFunction
  type(FunctionType)  :: solFunction
  testFunction = readFunctionData()
  call sortFunction(testFunction)

  print '(es24.16)', lagrange(testFunction)

end program
