#!/bin/bash

# $1 = TransferApp
# $2 = CompareApp
# $3 = Input file
# $4 = Test directory
# $5 = Comparison tolerance

echo "Transfer app: $1"
echo "Compare app: $2"
echo "Input file: $3"
echo "Test directory: $4"
echo "Tolerance: $5"

ps=`xml_grep 'patchfile' $4/$3 --text_only`
pfile1=`echo $ps | awk -F ' ' '{print $1}'`
pfile2=`echo $ps | awk -F ' ' '{print $2}'`
rs=`xml_grep 'resultfile' $4/$3 --text_only`
rfile1=`echo $rs | awk -F ' ' '{print $1}'`
rfile2=`echo $rs | awk -F ' ' '{print $2}'`
bs=`xml_grep 'boundaryfile' $4/$3 --text_only`
bfile1=`echo $bs | awk -F ' ' '{print $1}'`
bfile2=`echo $bs | awk -F ' ' '{print $2}'`

tmpdir=`mktemp -d -t simraXXXXXX`
cp $4/$3 $tmpdir
cp $4/$pfile1 $tmpdir
cp $4/$pfile2 $tmpdir
cp $4/$rfile1 $tmpdir
cp $4/$bfile1 $tmpdir

pushd $tmpdir
$1 $3
test $? -eq 0 || exit 1
$2 $rfile1 $rfile2 $bfile1 $bfile2 $pfile2 -tol $5
test $? -eq 0 || exit 1
popd
rm -Rf $tmpdir
