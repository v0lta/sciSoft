#!/bin/bash

#some color definitions
red='\033[0;31m'
green='\033[0;32m'
NC='\033[0m' # No Color

fortran=$(./testQuad.out)
sol=$(echo -e '0.66666667' )
diff -c  <(echo "$sol" ) <(echo "$fortran")

if [ "$sol" == "$fortran" ]; then
  echo -e "${green}quadpack::QAGS: passed!${NC}"
else
  echo -e "${red}quadpack::QAGS: passed!${NC}"
fi
