#!/bin/bash


#some color definitions
red='\033[0;31m'
green='\033[0;32m'
NC='\033[0m' # No Color

fort=$(echo -e $(./testApproxExp.out))
octv=$(echo -e $(octave --silent ./tests/ApproxTest/matlab/plainLagrange.m))


if [ "$fort" == "$octv" ]; then
  echo -e "${green}ApproxModule::Lagrange: passed!${NC}"
else
  echo -e "${red}UtilModule::Lagrange: failed!${NC}"
fi
