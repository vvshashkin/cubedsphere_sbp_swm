#!/bin/bash

TESTNAME="cmd line args test"
EXE=./$1/TEST_CMD_LINE
echo $TESTNAME $EXE
args=($EXE "Pozdravlyau" "Sharik," "ti" "balbes!")
i=0
for arg in `mpirun -n 2 $EXE ${args[@]:1}`
do
    [ $arg  == ${args[i]} ] || echo "failed"
    [ $arg  == ${args[i]} ] || exit
    let i=i+1
done

echo "passed"
