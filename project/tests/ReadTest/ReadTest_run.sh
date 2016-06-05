#!/bin/bash

#some color definitions
red='\033[0;31m'
green='\033[0;32m'
NC='\033[0m' # No Color

#commence testing.....
computed=$(echo -e $(./testRead.out < tests/ref1.in)| awk '{print tolower($0)}')
given=$(echo -e $(cat tests/ref1.in))


if [ "$computed" == "$given" ]; then
  echo -e "${green}ReadModule: ref1.in, passed!${NC}"
else
  echo -e "${red}ReadModule: ref1.in, failed!${NC}"
fi


computed=$(echo -e $(./testRead.out < tests/ref3.in)| awk '{print tolower($0)}')
given=$(echo -e $(cat tests/ref3.in))

if [ "$computed" == "$given" ]; then
  echo -e "${green}ReadModule::readFunctionData: ref3.in,passed!${NC}"
else
  echo -e "${red}ReadModule::readFunctionData: ref3.in,failed!${NC}"
fi
