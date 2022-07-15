#!/bin/bash

TESTNAME="Horizontal regrid test (regrid to latlon)"
EXE=./$1/TEST_REGRID
echo $TESTNAME $EXE

PREV_RESULTS=./$1"/src/test/test_regrid/reference_results.txt"

$EXE > results.txt && \
diff $PREV_RESULTS  results.txt > /dev/null && \
echo "serial $TESTNAME passed" || echo "parallel $TESTNAME failed" && \
diff -C1 $PREV_RESULTS results.txt
