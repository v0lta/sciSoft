#!/bin/bash

#make the files
make testWeight
make testGauss
make testRead
make testUtil
make testApproxExp
make testPoly
make testConv
make testOrthApprox
make testOrthogonality
make testQuad
make testQuadInt
make testLens

./tests/ExtTest/quadTest_run.sh
./tests/ReadTest/ReadTest_run.sh
./tests/IntegratorTest/IntTest_run.sh
./tests/UtilTest/utilTest_run.sh
./tests/UtilTest/convTest_run.sh
./tests/ApproxTest/ApproxTest_run.sh
./tests/ApproxTest/orthTest_run.sh
./tests/AsyTest/lensTest_run.sh


#w're done here... cleanup.
make clean
