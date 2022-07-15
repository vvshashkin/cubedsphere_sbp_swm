#!/bin/bash

TESTNAME="namelist read&broadcast test"
EXE=./$1/TEST_NAMELIST
echo $TESTNAME $EXE
mpirun -n 12 $EXE
