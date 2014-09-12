#!/bin/sh
rm -f examples/test.nc examples/test-index.nc
./examples/read -f examples/tos_O1_2001-2002.nc -n "tos" -m NETCDF -v $1 -g $2
./examples/read -f examples/tos_O1_2001-2002.nc -n "tos[0:10,1,1]" -m NETCDF -v $1 -g $2
./examples/writeData -f examples/test.nc -i examples/test.dat -n x -p "/time1" -d "4:3" -t double -c -m NETCDF -v $1 -g $2
./examples/writeData -f examples/test.nc -i examples/test.dat -n x -p "/" -d "4:3" -t int -c -m NETCDF -v $1 -g $2
./examples/writeAttr -f examples/test.nc -i examples/test.dat -n x -p "/time1" -a "attr" -l 2 -t int -c -m NETCDF -v $1 -g $2
./examples/buildIndex -f examples/test.nc -i examples/test-index.nc -r -l 2 -m NETCDF -v $1 -g $2
./examples/queryIndex -f examples/test.nc -i examples/test-index.nc -q "x>3 && x<8" -n x -p time1 -l 2 -m NETCDF -v $1 -g $2
./examples/queryIndex -f examples/test.nc -i examples/test-index.nc -q "x>3 && x<8" -n x -p time1 -l 2 -m NETCDF -b -v $1 -g $2
rm examples/test.nc examples/test-index.nc
