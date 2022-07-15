#!/bin/bash

TESTNAME="Horizontal differential operators test"
EXE=./$1/TEST_DIFFOPS_ALL
echo $TESTNAME $EXE

PREV_RESULTS=./$1"/src/test/test_diff_ops/reference_results.txt"

$EXE > results.txt && \
diff $PREV_RESULTS  results.txt > /dev/null && \
echo "serial $TESTNAME passed" || echo "parallel $TESTNAME failed" && \
diff -C1 $PREV_RESULTS results.txt

mpirun -n 48 $EXE > results_mpi.txt && \
diff $PREV_RESULTS results_mpi.txt > /dev/null && \
echo "parallel $TESTNAME passed" || echo "parallel $TESTNAME failed" && \
diff -C1 $PREV_RESULTS results_mpi.txt
