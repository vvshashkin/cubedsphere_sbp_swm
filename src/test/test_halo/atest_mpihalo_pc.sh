#!/bin/bash

TESTNAME="equiangular cubed sphere halo test (mpi)"
EXE=./$1/TEST_HALO_MAIN
echo $TESTNAME $EXE
mpirun -n 48 $EXE
