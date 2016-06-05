#!/bin/bash
#some color definitions
red='\033[0;31m'
green='\033[0;32m'
NC='\033[0m' # No Color

fortran=$(echo -e $(./testApprox.out < tests/ref3.in))
sol=$(echo -e $(cat tests/ApproxTest/matlab/fortTest.res))

diff -c  <(echo "$sol" ) <(echo "$fortran")

if [ "$fortran" == "$sol" ]; then
  echo -e "${green}GramSchmidtModule:monBasis: passed!${NC}"
else
  echo -e "${red}GramSchmidtModule:monBasis: failed!${NC}"
fi
