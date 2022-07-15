#!/bin/bash
#first arg is path to test executables 

TESTNAME="equiangular cubed sphere metric test"
EXE=.$1/TEST_METRIC_MAIN
echo $TESTNAME $EXE

$EXE
