#!/bin/bash


#some color definitions
red='\033[0;31m'
green='\033[0;32m'
NC='\033[0m' # No Color

touch tests/tmp.out
octave --quiet tests/UtilTest/generateRand.m > tests/tmp.out

sortedBash=$(echo -e $(sort -n tests/tmp.out))

sorted=$(echo -e $(./testUtil.out < tests/tmp.out))


if [ "$sortedBash" == "$sorted" ]; then
  echo -e "${green}UtilModule::Sort: passed!${NC}"
else
  echo -e "${red}UtilModule::Sort: failed!${NC}"
fi


polyRes=$(echo -e $(./testPolyVal.out))
if [ "$polyRes" == "24.750000000000000" ]; then
  echo -e "${green}UtilModule::polyVal: passed!${NC}"
else
  echo -e "${red}UtilModule::polyVal: failed!${NC}"
fi



rm tests/tmp.out
