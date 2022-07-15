#!/bin/bash
#first arg is path to test executables

TESTNAME="time stepping schemes test"
EXE=.$1/TEST_TS
echo $TESTNAME $EXE

$EXE
