#!/bin/bash

TESTNAME="exchange test"
EXE=./$1/TEST_EXCH_MAIN
echo $TESTNAME $EXE
mpirun -n 48 $EXE
