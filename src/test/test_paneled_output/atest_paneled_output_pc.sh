#!/bin/bash

TESTNAME="paneled output test"
EXE=./$1/TEST_PANELED_OUTPUT
echo $TESTNAME "(" $EXE")"
mpirun -n 24 $EXE
