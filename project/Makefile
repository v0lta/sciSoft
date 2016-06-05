#---------------------------gfortran------------------------------------------------------
FC=gfortran
flags= -Og -fimplicit-none -Warray-bounds -Wconversion -pedantic -std=f2008 -freal-4-real-8
#flags= -fimplicit-none  -std=f2008 -freal-4-real-8 -O3
#------------everything works fine--------------------------------------------------------

#---------------------------nagfor--------------------------------------------------------
#FC=nagfor
#flags= -r8 -kind=byte -colour=error:red -gline
#flags= -r8 -kind=byte -ieee=full -O3
#------------------works------------------------------------------------------------------

#---------------------------ifort---------------------------------------------------------
#FC = ifort
#flags = -r8
#compiles but does not work. Ifort does not understand the public variable construct I used
#to implement functions for quadpack probably an f2003 problem.

#g95 does not compile due to missing f2003 features.

#----------------------------mpif90-------------------------------------------------------
#FC=mpif90
#flags= -freal-4-real-8
#-------ok--------------------------------------------------------------------------


lpack=-L/usr/lib/ -llapack  -lblas

basis: ReadModule.o UtilModule.o FunctionModule.o IntegrationModule.o \
       quadpack.o LagrangeModule.o GramSchmidtModule.o ProjectionModule.o \
       ParseModule.o FavardModule.o AsyModule.o
	$(FC) $(flags) -c ./src/basis.f90
	$(FC) $(flags) -o basis basis.o ReadModule.o UtilModule.o FunctionModule.o \
			IntegrationModule.o LagrangeModule.o GramSchmidtModule.o AsyModule.o \
			ProjectionModule.o ParseModule.o quadpack.o FavardModule.o $(lpack)

testLens:  AsyModule.o UtilModule.o
	$(FC) $(flags) -c ./tests/AsyTest/lensTest.f90
	$(FC) $(flags) -o lensTest.out lensTest.o UtilModule.o AsyModule.o

testAsyExpl: AsyModule.o UtilModule.o IntegrationModule.o
	$(FC) $(flags) -c ./tests/AsyTest/asyTestExplicit.f90
	$(FC) $(flags) -o asyTestExpl.out asyTestExplicit.o UtilModule.o AsyModule.o \
										IntegrationModule.o quadpack.o $(lpack)

testAsy: AsyModule.o UtilModule.o
	$(FC) $(flags) -c ./tests/AsyTest/asyTest.f90
	$(FC) $(flags) -o asyTest.out asyTest.o UtilModule.o AsyModule.o

testShape: AsyModule.o UtilModule.o
	$(FC) $(flags) -c ./tests/AsyTest/testShape.f90
	$(FC) $(flags) -o shapeTest.out testShape.o UtilModule.o AsyModule.o

testQuadInt: quadpack.o IntegrationModule.o UtilModule.o
	$(FC) $(flags) -c ./tests/IntegratorTest/quadIntegralTest.f90
	$(FC) $(flags) -o testIntQuad.out quadpack.o quadIntegralTest.o IntegrationModule.o \
		          UtilModule.o $(lpack)

testQuad: quadpack.o
	$(FC) $(flags) -c ./tests/ExtTest/quadTest.f90
	$(FC) $(flags) -o testQuad.out quadpack.o quadTest.o

testFavard: FavardModule.o UtilModule.o IntegrationModule.o quadpack.o GramSchmidtModule.o
	$(FC) $(flags) -c ./tests/AsyTest/testFavard.f90
	$(FC) $(flags) -o testFavard.out testFavard.o UtilModule.o FavardModule.o \
			  IntegrationModule.o quadpack.o GramSchmidtModule.o $(lpack)

testElia: IntegrationModule.o UtilModule.o FunctionModule.o ReadModule.o \
	 				GramSchmidtModule.o LagrangeModule.o
	$(FC) $(flags) -c ./tests/elia/eliaTest.f90
	$(FC) $(flags) -o eliaTest.out eliaTest.o UtilModule.o GramSchmidtModule.o \
			  FunctionModule.o IntegrationModule.o ReadModule.o LagrangeModule.o $(lpack)

testEvalMat: GramSchmidtModule.o UtilModule.o IntegrationModule.o ProjectionModule.o
	$(FC) $(flags) -c ./tests/ApproxTest/evalMatTest.f90
	$(FC) $(flags) -o evalMatTest.out evalMatTest.o GramSchmidtModule.o UtilModule.o \
			  IntegrationModule.o ProjectionModule.o $(lpack)

testParse: ParseModule.o
	$(FC) $(flags) -c ./tests/UtilTest/parseTest.f90
	$(FC) $(flags) -o basis parseTest.o ParseModule.o

testTrapz : IntegrationModule.o  UtilModule.o FunctionModule.o
	$(FC) $(flags) -c ./tests/IntegratorTest/trapzTest.f90
	$(FC) $(flags) -o trapzTest.out trapzTest.o IntegrationModule.o  UtilModule.o \
		          FunctionModule.o $(lpack)

testProj: IntegrationModule.o UtilModule.o FunctionModule.o LagrangeModule.o\
	  ReadModule.o GramSchmidtModule.o

testLagApprox: IntegrationModule.o UtilModule.o FunctionModule.o  quadpack.o\
		ReadModule.o GramSchmidtModule.o LagrangeModule.o
	$(FC) $(flags) -c ./tests/ApproxTest/lagApproxTest.f90
	$(FC) $(flags) -o lagApprox.out lagApproxTest.o UtilModule.o GramSchmidtModule.o \
	 		  FunctionModule.o IntegrationModule.o ReadModule.o LagrangeModule.o \
			  quadpack.o $(lpack)

testLagGen: UtilModule.o FunctionModule.o LagrangeModule.o
	$(FC) $(flags) -c ./tests/ApproxTest/lagrangeTest.f90
	$(FC) $(flags) -o testLagGen.out lagrangeTest.o UtilModule.o FunctionModule.o LagrangeModule.o

testMatEval: UtilModule.o GramSchmidtModule.o IntegrationModule.o
	$(FC) $(flags) -c ./tests/ApproxTest/testEvalMat.f90
	$(FC) $(flags) -o testEval.out testEvalMat.o UtilModule.o GramSchmidtModule.o \
			  IntegrationModule.o $(lpack)

testOrthApprox: IntegrationModule.o UtilModule.o FunctionModule.o \
								 ReadModule.o GramSchmidtModule.o quadpack.o
	$(FC) $(flags) -c ./tests/ApproxTest/monApproxTest.f90
	$(FC) $(flags) -o testApprox.out monApproxTest.o UtilModule.o GramSchmidtModule.o \
	 		  FunctionModule.o IntegrationModule.o ReadModule.o quadpack.o \
				$(lpack)

testIntOmegaPoly: GramSchmidtModule.o IntegrationModule.o UtilModule.o FunctionModule.o
	$(FC) $(flags) -c ./tests/IntegratorTest/testIntOmegaPoly.f90
	$(FC) $(flags) -o testOmegaPoly.out testIntOmegaPoly.o UtilModule.o GramSchmidtModule.o \
	 		  FunctionModule.o IntegrationModule.o $(lpack)

testOrthogonality: IntegrationModule.o UtilModule.o FunctionModule.o \
		   ReadModule.o GramSchmidtModule.o LagrangeModule.o quadpack.o
	$(FC) $(flags) -c ./tests/ApproxTest/checkOrthogonality.f90
	$(FC) $(flags) -o testOrth.out checkOrthogonality.o UtilModule.o GramSchmidtModule.o \
	 		  FunctionModule.o IntegrationModule.o ReadModule.o  \
			  quadpack.o LagrangeModule.o $(lpack)

testConv: UtilModule.o
	$(FC) $(flags) -c ./tests/UtilTest/convTest.f90
	$(FC) $(flags) -o testConv.out convTest.o UtilModule.o

testPoly: UtilModule.o
	$(FC) $(flags) -c ./tests/UtilTest/polyvalTest.f90
	$(FC) $(flags) -o testPolyVal.out polyvalTest.o UtilModule.o


testApprox: LagrangeModule.o UtilModule.o FunctionModule.o ReadModule.o quadpack.o
	$(FC) $(flags) -c ./tests/ApproxTest/approxTestRef3.f90
	$(FC) $(flags) -c ./tests/ApproxTest/readSol.f90
	$(FC) $(flags) -o testApproxRef3.out approxTestRef3.o UtilModule.o FunctionModule.o \
		          LagrangeModule.o ReadModule.o
	$(FC) $(flags) -o readSol.out readSol.o UtilModule.o FunctionModule.o \
	         	  LagrangeModule.o ReadModule.o quadpack.o

testApproxExp: LagrangeModule.o UtilModule.o FunctionModule.o
	$(FC) $(flags) -c ./tests/ApproxTest/approxTestExp.f90
	$(FC) $(flags) -o testApproxExp.out approxTestExp.o UtilModule.o FunctionModule.o \
			  LagrangeModule.o

testUtil: UtilModule.o
	$(FC) $(flags) -c ./tests/UtilTest/utilTest.f90
	$(FC) $(flags) -o testUtil.out utilTest.o UtilModule.o

testGauss2: IntegrationModule.o UtilModule.o FunctionModule.o
		$(FC) $(flags) -c ./tests/IntegratorTest/gaussTest2.f90 $(lpack)
		$(FC) $(flags) -o testGauss.out gaussTest2.o IntegrationModule.o \
		 		UtilModule.o quadpack.o $(lpack)

testGauss: IntegrationModule.o UtilModule.o FunctionModule.o ReadModule.o
	$(FC) $(flags) -c ./src/Integrator/IntegrationModule.f90 IntegrationModule.o  \
			FunctionModule.o ReadModule.o $(lpack)
	$(FC) $(flags) -c ./tests/IntegratorTest/gaussTest.f90 $(lpack)
	$(FC) $(flags) -o testGauss.out gaussTest.o IntegrationModule.o FunctionModule.o \
	 		UtilModule.o ReadModule.o quadpack.o $(lpack)

testWeight: IntegrationModule.o
	$(FC) $(flags) -c ./tests/IntegratorTest/weightTest.f90 $(lpack)
	$(FC) $(flags) -o testWeight.out weightTest.o IntegrationModule.o UtilModule.o \
			  quadpack.o $(lpack)


testRead: ReadModule.o
	$(FC) $(flags) -c ./tests/ReadTest/ReadTest.f90
	$(FC) $(flags) -o testRead.out FunctionModule.o UtilModule.o ReadTest.o ReadModule.o

testFunction: FunctionModule.o
	$(FC) $(flags) -c ./tests/FunctionTest/functionTest.f90
	$(FC) $(flags) -o function FunctionModule.o UtilModule.o functionTest.o

quadpack.o:
	$(FC) $(flags) -c ./src/External/quadpack.f90


FavardModule.o: FunctionModule.o IntegrationModule.o UtilModule.o
		$(FC) $(flags) -c ./src/Approx/FavardModule.f90

AsyModule.o: UtilModule.o
	$(FC) $(flags) -c ./src/Approx/AsyModule.f90

ParseModule.o:
	$(FC) $(flags) -c ./src/UtilModule/ParseModule.f90

ProjectionModule.o: UtilModule.o LagrangeModule.o GramSchmidtModule.o
	$(FC) $(flags) -c ./src/Approx/ProjectionModule.f90

GramSchmidtModule.o: UtilModule.o FunctionModule.o IntegrationModule.o
	$(FC) $(flags) -c ./src/Approx/GramSchmidtModule.f90

ReadModule.o: FunctionModule.o
	$(FC) $(flags) -c ./src/ReadIn/ReadModule.f90

LagrangeModule.o: UtilModule.o FunctionModule.o
	$(FC) $(flags) -c ./src/Approx/LagrangeModule.f90

IntegrationModule.o: UtilModule.o FunctionModule.o quadpack.o
	$(FC) $(flags) -c ./src/Integrator/IntegrationModule.f90 $(lpack)

FunctionModule.o: UtilModule.o
	$(FC) $(flags) -c ./src/Function/FunctionModule.f90

UtilModule.o:
	$(FC) $(flags) -c ./src/UtilModule/UtilModule.f90

clean:
	rm *.mod *.o *.out
