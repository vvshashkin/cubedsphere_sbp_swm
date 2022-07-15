#!/bin/bash

test_exe_path="."
test_res_dir="testres"
test_res_short="test_short.txt"
test_res_long="test_long.txt"
[ $1 ] && test_res_dir=$1

srcdir=src
genmake=genMakefile.sh
makecmd="make -k all" # -k to compile & test all that is compileable at the moment
                      # despite of the some test's compilation is failed

pref="atest"
suff="pc"

cleanexit(){
    msgall $@
    rm `ls | grep -Ev "Makefile|make.log|$test_res_short|$test_res_long"` -rf
    #rm $srcdir obj mod -rf
    exit
}
msgsh_n_io(){
    echo $@
    echo $@ >> $test_res_short
}
msglong(){
    echo $@ >> $test_res_long
}
msgall(){
    msgsh_n_io $@
    msglong $@
}

[ -e $test_res_dir ] && echo $test_res_dir "already exists, please remove it or use another directory for test results" && exit
mkdir $test_res_dir
cd $test_res_dir
if [ $? -ne 0 ]
then
    echo "cannot cd to " $test_res_dir "exit!" ; exit
fi

echo > $test_res_short
echo > $test_res_long
msgall `date`

cp -r ../$srcdir $srcdir


../$genmake || cleanexit "cannot generate Makefile. Exit!"

echo "making..."
$makecmd &>make.log || msgall "Warning: make errors, see make.log"
echo "make complete!"

for tst in `ls ../$srcdir/test/*/$pref*$suff.sh`
do
    $tst &> tst.tmp
    stat=$?

    msglong ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
    msglong "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"
    cat tst.tmp >> $test_res_long

    tname=`head -n1 tst.tmp`
    tres="undefined"
    grep "pass" tst.tmp &>/dev/null && tres="passed"
    grep "fail" tst.tmp &>/dev/null && tres="failed"
    [ $stat -ne 0 ]         && tres="failed"
    msgsh_n_io $tname ":::" $tres
done

msgall ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
msgall "tests completed, please see results in " $test_res_dir"/"$test_res_short", "$test_res_dir"/"$test_res_long

cleanexit
