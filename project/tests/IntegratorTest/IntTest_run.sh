#!/bin/bash

#some color definitions
red='\033[0;31m'
green='\033[0;32m'
NC='\033[0m' # No Color

#commence testing.....
computed=$(echo -e $(./testWeight.out)| awk '{print tolower($0)}')
given=$(echo -e $(octave --quiet tests/IntegratorTest/scriptWeight.m))

if [ "$computed" == "$given" ]; then
  echo -e "${green}IntegrationModule::rJacobi: passed!${NC}"
else
  echo -e "${red}IntegrationModule::rJacobi: failed!${NC}"
fi

intResult=$(echo -e $(./testGauss.out < tests/ref1.in))

if [ "$intResult" == "4.993379538113" ]; then
  echo -e "${green}IntegrationModule::gauss: passed!${NC}"
else
  echo -e "${red}IntegrationModule::gauss: failed!${NC}"
fi

intResult=$(./testIntQuad.out)

if [ "$intResult" == "2.00000000" ]; then
  echo -e "${green}IntegrationModule::quadIntegral: passed!${NC}"
else
  echo -e "${red}IntegrationModule::quadIntegral: failed!${NC}"
fi
