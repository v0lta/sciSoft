
#some color definitions
red='\033[0;31m'
green='\033[0;32m'
NC='\033[0m' # No Color

echo -e "sorry for the delay matlab is starting up..."
cd ./tests/AsyTest/matlab/
matlab -nosplash -nodesktop -r combine2 > dat.out
cd ../../..



matlab=$(echo -e $(sed '$d;1d;/[a-zA-Z]/d;s/i//g;s/+//g;s/- /-/g' ./tests/AsyTest/matlab/dat.out))
fortran=$(echo -e $(./lensTest.out))

diff <(echo "$matlab" ) <(echo "$fortran")

if [ "$matlab" == "$fortran" ]; then
  echo -e "${green}AsyModule::lens: passed!${NC}"
else
  echo -e "${red}AsyModule::lens: failed!${NC}"
fi


#make clean
