


red='\033[0;31m'
green='\033[0;32m'
NC='\033[0m' # No Color

fortran=$(echo -e $(./testOrth.out))

if [ "$fortran" == "T" ]; then
  echo -e "${green}gramSchmidt::GramSchmidt: ok.${NC}"
else
  echo -e "${red}gramSchmidt::GramSchmidt: failed.${NC}"
fi


#make clean
#rm monOrth.dat
