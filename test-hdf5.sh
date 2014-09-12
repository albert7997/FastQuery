#!/bin/sh
rm -f examples/test.h5 examples/test-index.h5 
#(cd examples && make -j 2 nompi)
./examples/read -f examples/h5uc-data.h5 -n px -p TimeStep1 -v $1 -g $2
./examples/read -f examples/h5uc-data.h5 -n "px[0:10]" -p TimeStep1 -v $1 -g $2
./examples/writeData -f examples/test.h5 -i examples/test.dat -n x -p "/time1" -d "4:3" -t double -c -v $1 -g $2
./examples/writeData -f examples/test.h5 -i examples/test.dat -n x -p "/" -d "4:3" -t int -c -v $1 -g $2
./examples/writeAttr -f examples/test.h5 -i examples/test.dat -n x -p "/time1" -a "attr" -l 2 -t int -c -v $1 -g $2
./examples/buildIndex -f examples/test.h5 -i examples/test-index.h5 -r -l 2 -v $1 -g $2
./examples/queryIndex -f examples/test.h5 -i examples/test-index.h5 -q "x>3 && x<8" -n x -p time1 -l 2 -v $1 -g $2
./examples/queryIndex -f examples/test.h5 -i examples/test-index.h5 -q "x>3 && x<8" -n x -p time1 -l 2 -b -v $1 -g $2
./examples/histogram -x 1 -y 2 -e 9 -s 2 -f examples/test.h5 -i examples/test-index.h5 -q "x>1 && x<10" -n x -p time1 -l 2 -v $1 -g $2
rm examples/test.h5 examples/test-index.h5
