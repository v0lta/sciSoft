#!/bin/bash


#some color definitions
red='\033[0;31m'
green='\033[0;32m'
NC='\033[0m' # No Color

fort=$(echo -e $(./testApproxRef3.out < ./tests/ref3.in))
octv=$(echo -e $(./readSol.out < ./tests/ref3.test))

if [ "$fort" == "$octv" ]; then
  echo -e "${green}ApproxModule::Lagrange: passed!${NC}"
else
  echo -e "${red}ApproxModule::Lagrange: failed!${NC}"
fi
