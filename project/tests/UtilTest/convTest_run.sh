#!/bin/bash


#some color definitions
red='\033[0;31m'
green='\033[0;32m'
NC='\033[0m' # No Color

octave=$(echo -e $(./tests/UtilTest/conv/octTest.sh))
fortran=$(echo -e $(echo -e $(cat tmpP1.dat) "\n" $(cat tmpP2.dat)| ./testConv.out))
#diff -c  <(echo "$octave" ) <(echo "$fortran")

if [ "$octave" == "$fortran" ]; then
  echo -e "${green}UtilModule::conv: passed!${NC}"
else
  echo -e "${red}UtilModule::conv: failed!${NC}"
fi

rm *.dat
