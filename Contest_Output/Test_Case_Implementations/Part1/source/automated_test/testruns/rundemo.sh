#!/bin/sh

preptest () {
	testn=$1
	cp $testn.wam fnames.wam
	filelist="$testn.out $testn.fpt $testn.idf $testn.pnl $testn.p2f $testn.5pb $testn.5v? $testn.?? error*.$testn $testn'_low.gdf' $testn'_pan.gdf' $testn'_pan.dat' $testn'_pat.dat' $testn.crp wamitlog.$testn"
	for file in $filelist
	do
		if [ -f $file ] ; then
		       	rm $file
		fi
	done
}

runwamit () {
	testn=$1
	wamit_demo
}

postwamit () {
	testn=$1
	mv errorp.log errorp.$testn
	mv errorf.log errorf.$testn
	mv wamitlog.txt wamitlog.$testn
}

testlist='test01 test01a test02 test03 test04 test05 test05a test06 test07 test08 test09 test09a test11 test11a test11b test12 test13 test13a test14 test14a test15 test16 test16a test17 test17a test17b test17c test18 test19 test20 test21 test22 test23 test24 test25'

for testcase in $testlist
do
	echo "preparing to run $testcase"
	preptest $testcase
	runwamit $testcase 
	postwamit $testcase
done
