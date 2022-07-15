#!/bin/bash

TESTNAME="Latlon output test"

EXE=./$1/TEST_LATLON_OUTPUT
echo $TESTNAME $EXE

mpirun -n 12 $EXE
