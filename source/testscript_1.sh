#! /bin/bash
# [Topcoder]
# This script consists of Test-case from 7 to 11 
exec >> automated_test.log

# I don't want to stop script when shell script returns zero
set -e
# http://stackoverflow.com/questions/59895/can-a-bash-script-tell-what-directory-its-stored-in?page=1&tab=votes#tab-top
SOURCE="${BASH_SOURCE[0]}"
while [ -h "$SOURCE" ]; do # resolve $SOURCE until the file is no longer a symlink
  DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"
  SOURCE="$(readlink "$SOURCE")"
  [[ $SOURCE != /* ]] && SOURCE="$DIR/$SOURCE" # if $SOURCE was a relative symlink, we need to resolve it relative to the path where the symlink file was located
done
DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"

INSTALL_PATH="$DIR"
ROOT="$DIR"

#cd "$INSTALL_PATH"/openwarpgui
cd "$INSTALL_PATH"

chmod +x "$INSTALL_PATH"/openwarpgui/bundled/mesh-generator/bin/nglib_mesh

chmod +x "$INSTALL_PATH"/openwarpgui/bundled/paraview_linux/bin/paraview

export LD_LIBRARY_PATH="$ROOT/openwarpgui/bundled/simulation/libs"
export LDFLAGS="-L$ROOT/openwarpgui/bundled/simulation/libs"
export PYTOHNPATH="$INSTALL_PATH"/openwarpgui

mkdir -p ~/OpenWarpFiles/temp
mkdir -p ~/OpenWarpFiles/user_data

echo "Starting Test Script 1 "


# Automated Test Script and Deployment guide #-------------------------------------------------------------------------------------
# Test1. Created GDFMesh.py to convert .GDF files to .dat files. (for runing GDFMesh, matlab engine is required, see deployment guide).
#	 -f is File Name
# 	 -t is for 
# 	 -g is Gravity Center
#	 -p is for
   echo "-------------------------starting Test1--------------------------------"	 
#  python Test1/GDFMesh.py -f 'testruns/test01.gdf' -t 0 -g [0,0,0] -p 0
#  Currently, there is an compilation error with GDFMesh.m, so it will throw errors 

#  RunNemoh_v5.m 
#  python Test1/RunNemoh_v5.py 

#  bemio

# These will create a new file named filename_bemio_output in the test run folder   
  echo "Using Bemio"
python automated_test/Test1/GDFtoDAT.py 'automated_test/testruns/test01.gdf'
python automated_test/Test1/GDFtoDAT.py 'automated_test/testruns/test02.gdf'
python automated_test/Test1/GDFtoDAT.py 'automated_test/testruns/test03.gdf'
python automated_test/Test1/GDFtoDAT.py 'automated_test/testruns/test04.gdf'
python automated_test/Test1/GDFtoDAT.py 'automated_test/testruns/test05.gdf'
python automated_test/Test1/GDFtoDAT.py 'automated_test/testruns/test05a.gdf'
python automated_test/Test1/GDFtoDAT.py 'automated_test/testruns/test06.gdf'
python automated_test/Test1/GDFtoDAT.py 'automated_test/testruns/test07.gdf'
python automated_test/Test1/GDFtoDAT.py 'automated_test/testruns/test09.gdf'
python automated_test/Test1/GDFtoDAT.py 'automated_test/testruns/test09a.gdf'
python automated_test/Test1/GDFtoDAT.py 'automated_test/testruns/test11.gdf'
python automated_test/Test1/GDFtoDAT.py 'automated_test/testruns/test11b.gdf'
python automated_test/Test1/GDFtoDAT.py 'automated_test/testruns/test11c.gdf'
python automated_test/Test1/GDFtoDAT.py 'automated_test/testruns/test13ac.gdf'
python automated_test/Test1/GDFtoDAT.py 'automated_test/testruns/test14.gdf'
python automated_test/Test1/GDFtoDAT.py 'automated_test/testruns/test15.gdf'
python automated_test/Test1/GDFtoDAT.py 'automated_test/testruns/test17.gdf'
python automated_test/Test1/GDFtoDAT.py 'automated_test/testruns/test17c.gdf'
python automated_test/Test1/GDFtoDAT.py 'automated_test/testruns/test18.gdf'
python automated_test/Test1/GDFtoDAT.py 'automated_test/testruns/test19.gdf'
python automated_test/Test1/GDFtoDAT.py 'automated_test/testruns/test20.gdf'
python automated_test/Test1/GDFtoDAT.py 'automated_test/testruns/test21.gdf'
python automated_test/Test1/GDFtoDAT.py 'automated_test/testruns/test22.gdf'
python automated_test/Test1/GDFtoDAT.py 'automated_test/testruns/test22a.gdf'

echo "--------------------------finishing Test1 ---------------------------------"
#--------------------------------------------------------------------------------------------------------------------
echo "--------------------------starting Test2 -------------------------------"

	# running Openwarpgui in background, so that I can use openwarp_cli
	# It will open a browser.
	# source run.sh &
	
	# wait for 3 seconds 	
	# sleep 5
	
	chmod u+x automated_test/testruns/openwarp/
	#Provide filename, not extension '.out'
	# Create .H5 files from all Wamit's .out file 
	python automated_test/Test2/OUTtoHDF5.py 'automated_test/testruns/out/test01'
	python automated_test/Test2/OUTtoHDF5.py 'automated_test/testruns/out/test01a'
	python automated_test/Test2/OUTtoHDF5.py 'automated_test/testruns/out/test02'
	python automated_test/Test2/OUTtoHDF5.py 'automated_test/testruns/out/test03'
	python automated_test/Test2/OUTtoHDF5.py 'automated_test/testruns/out/test04'
	python automated_test/Test2/OUTtoHDF5.py 'automated_test/testruns/out/test05'
	python automated_test/Test2/OUTtoHDF5.py 'automated_test/testruns/out/test05a'
	python automated_test/Test2/OUTtoHDF5.py 'automated_test/testruns/out/test06'
	python automated_test/Test2/OUTtoHDF5.py 'automated_test/testruns/out/test07'
	python automated_test/Test2/OUTtoHDF5.py 'automated_test/testruns/out/test09'
	python automated_test/Test2/OUTtoHDF5.py 'automated_test/testruns/out/test09a'
	python automated_test/Test2/OUTtoHDF5.py 'automated_test/testruns/out/test11'
	python automated_test/Test2/OUTtoHDF5.py 'automated_test/testruns/out/test11b'
	python automated_test/Test2/OUTtoHDF5.py 'automated_test/testruns/out/test11c'
	python automated_test/Test2/OUTtoHDF5.py 'automated_test/testruns/out/test12'
	python automated_test/Test2/OUTtoHDF5.py 'automated_test/testruns/out/test13'
	python automated_test/Test2/OUTtoHDF5.py 'automated_test/testruns/out/test13a'
	python automated_test/Test2/OUTtoHDF5.py 'automated_test/testruns/out/test14'
	python automated_test/Test2/OUTtoHDF5.py 'automated_test/testruns/out/test14a'
	python automated_test/Test2/OUTtoHDF5.py 'automated_test/testruns/out/test15'
	python automated_test/Test2/OUTtoHDF5.py 'automated_test/testruns/out/test17'
	python automated_test/Test2/OUTtoHDF5.py 'automated_test/testruns/out/test17c'
	python automated_test/Test2/OUTtoHDF5.py 'automated_test/testruns/out/test18'
	python automated_test/Test2/OUTtoHDF5.py 'automated_test/testruns/out/test19'
	python automated_test/Test2/OUTtoHDF5.py 'automated_test/testruns/out/test20'
	python automated_test/Test2/OUTtoHDF5.py 'automated_test/testruns/out/test21'
	python automated_test/Test2/OUTtoHDF5.py 'automated_test/testruns/out/test22'
	python automated_test/Test2/OUTtoHDF5.py 'automated_test/testruns/out/test22a'
	
	
	# Create HDF5 file from recently created .dat files inside /testruns/openwarp/openwarpH5 folder.
	# Same Configuration, different geometry ! 	
	python openwarpgui/openwarp_cli.py automated_test/Test2/test2.json
	
	# Compare hdf5 file using H5DIFF tool:
	h5diff -vc automated_test/Test2/cli_results/test01/simulation/db.hdf5 automated_test/testruns/openwarp/test01.h5 >automated_test/Test2/test01.log
	
	h5diff -vc automated_test/Test2/cli_results/test01a/simulation/db.hdf5 automated_test/testruns/openwarp/test01a.h5 >automated_test/Test2/test01a.log
	
	h5diff -vc automated_test/Test2/cli_results/test02/simulation/db.hdf5 automated_test/testruns/openwarp/test02.h5 >automated_test/Test2/test02.log
	
	h5diff -vc automated_test/Test2/cli_results/test03/simulation/db.hdf5 automated_test/testruns/openwarp/test03.h5 >automated_test/Test2/test03.log
	
	h5diff -vc automated_test/Test2/cli_results/test04/simulation/db.hdf5 automated_test/testruns/openwarp/test04.h5 >automated_test/Test2/test04.log
	
	h5diff -vc automated_test/Test2/cli_results/test05/simulation/db.hdf5 automated_test/testruns/openwarp/test05.h5 >automated_test/Test2/test05.log
	
	h5diff -vc automated_test/Test2/cli_results/test05a/simulation/db.hdf5 automated_test/testruns/openwarp/test05a.h5 >automated_test/Test2/test05a.log
	
	h5diff -vc automated_test/Test2/cli_results/test06/simulation/db.hdf5 automated_test/testruns/openwarp/test06.h5 >automated_test/Test2/test06.log

	h5diff -vc automated_test/Test2/cli_results/test07/simulation/db.hdf5 automated_test/testruns/openwarp/test07.h5 >automated_test/Test2/test07.log
	
	h5diff -vc automated_test/Test2/cli_results/test09/simulation/db.hdf5 automated_test/testruns/openwarp/test09.h5 >automated_test/Test2/test09.log

	h5diff -vc automated_test/Test2/cli_results/test09a/simulation/db.hdf5 automated_test/testruns/openwarp/test09a.h5 >automated_test/Test2/test09a.log

	h5diff -vc automated_test/Test2/cli_results/test11/simulation/db.hdf5 automated_test/testruns/openwarp/test11.h5 >automated_test/Test2/test11.log
	
	h5diff -vc automated_test/Test2/cli_results/test012/simulation/db.hdf5 automated_test/testruns/openwarp/test12.h5 >automated_test/Test2/test12.log

	h5diff -vc automated_test/Test2/cli_results/test11b/simulation/db.hdf5 automated_test/testruns/openwarp/test11b.h5 >automated_test/Test2/test11b.log

	h5diff -vc automated_test/Test2/cli_results/test11c/simulation/db.hdf5 automated_test/testruns/openwarp/test11c.h5 >automated_test/Test2/test11c.log
	
	h5diff -vc automated_test/Test2/cli_results/test12/simulation/db.hdf5 automated_test/testruns/openwarp/test12.h5 >automated_test/Test2/test12.log

	h5diff -vc automated_test/Test2/cli_results/test13/simulation/db.hdf5 automated_test/testruns/openwarp/test13.h5 >automated_test/Test2/test13.log

	h5diff -vc automated_test/Test2/cli_results/test13a/simulation/db.hdf5 automated_test/testruns/openwarp/test13a.h5 >automated_test/Test2/test13a.log

	h5diff -vc automated_test/Test2/cli_results/test14/simulation/db.hdf5 automated_test/testruns/openwarp/test14.h5 >automated_test/Test2/test14.log

	h5diff -vc automated_test/Test2/cli_results/test14a/simulation/db.hdf5 automated_test/testruns/openwarp/test14a.h5 >automated_test/Test2/test14a.log

	h5diff -vc automated_test/Test2/cli_results/test15/simulation/db.hdf5 automated_test/testruns/openwarp/test15.h5 >automated_test/Test2/test15.log
	
	h5diff -vc automated_test/Test2/cli_results/test17/simulation/db.hdf5 automated_test/testruns/openwarp/test17.h5 >automated_test/Test2/test17.log
	
	h5diff -vc automated_test/Test2/cli_results/test17c/simulation/db.hdf5 automated_test/testruns/openwarp/test17c.h5 >automated_test/Test2/test17c.log

	h5diff -vc automated_test/Test2/cli_results/test18/simulation/db.hdf5 automated_test/testruns/openwarp/test18.h5 >automated_test/Test2/test18.log

	h5diff -vc automated_test/Test2/cli_results/test19/simulation/db.hdf5 automated_test/testruns/openwarp/test19.h5 >automated_test/Test2/test19.log

	h5diff -vc automated_test/Test2/cli_results/test20/simulation/db.hdf5 automated_test/testruns/openwarp/test20.h5 >automated_test/Test2/test20.log

	h5diff -vc automated_test/Test2/cli_results/test21/simulation/db.hdf5 automated_test/testruns/openwarp/test21.h5 >automated_test/Test2/test21.log

	h5diff -vc automated_test/Test2/cli_results/test22/simulation/db.hdf5 automated_test/testruns/openwarp/test22.h5 >automated_test/Test2/test22.log
	
	h5diff -vc automated_test/Test2/cli_results/test22a/simulation/db.hdf5 automated_test/testruns/openwarp/test22a.h5 >automated_test/Test2/test22a.log

	

	WAMIT_H5 = '/automated_test/testruns/openwarp/'
	OPENWARP_H5 = '/automated_test/testruns/openwarp/openwarpH5'
	
	# Compare hdf5 file using H5DIFF tool : 
	#h5diff -vc ${WAMIT_H5}/test01.h5 ${OPENWARP_H5}/db_test01.h5
	#pwd	
	#cd automated_test/testruns/openwarp
	#pwd

	
	
echo "---------------------------finishing Test2--------------------------------"
   
#---------------------------------------------------------------------------------------------------------------------
# Test3. read test3.json to that enables USE_HIGHER_ORDER ad TRUE with B_SPLINE_ORDER as 1. Later, it converts 
echo "---------------------------Starting  Test3--------------------------------"
	python openwarpgui/openwarp_cli.py automated_test/Test3/test3.json
	
	h5diff -vc automated_test/Test3/cli_results/test11/simulation/db.hdf5 automated_test/testruns/openwarp/test11.h5 >automated_test/Test3/test11.log

	h5diff -vc automated_test/Test3/cli_results/test11a/simulation/db.hdf5 automated_test/testruns/openwarp/test11a.h5 >automated_test/Test3/test11a.log

	h5diff -vc automated_test/Test3/cli_results/test012/simulation/db.hdf5 automated_test/testruns/openwarp/test12.h5 >automated_test/Test3/test12.log

	h5diff -vc automated_test/Test3/cli_results/test11b/simulation/db.hdf5 automated_test/testruns/openwarp/test11b.h5 >automated_test/Test3/test11b.log

	h5diff -vc automated_test/Test3/cli_results/test11c/simulation/db.hdf5 automated_test/testruns/openwarp/test11c.h5 >automated_test/Test3/test11c.log
	
	h5diff -vc automated_test/Test3/cli_results/test12/simulation/db.hdf5 automated_test/testruns/openwarp/test12.h5 >automated_test/Test3/test12.log

	h5diff -vc automated_test/Test3/cli_results/test13/simulation/db.hdf5 automated_test/testruns/openwarp/test13.h5 >automated_test/Test3/test13.log

	h5diff -vc automated_test/Test3/cli_results/test13a/simulation/db.hdf5 automated_test/testruns/openwarp/test13a.h5 >automated_test/Test3/test13a.log

	h5diff -vc automated_test/Test3/cli_results/test14/simulation/db.hdf5 automated_test/testruns/openwarp/test14.h5 >automated_test/Test3/test14.log

	h5diff -vc automated_test/Test3/cli_results/test14a/simulation/db.hdf5 automated_test/testruns/openwarp/test14a.h5 >automated_test/Test3/test14a.log

	h5diff -vc automated_test/Test3/cli_results/test15/simulation/db.hdf5 automated_test/testruns/openwarp/test15.h5 >automated_test/Test3/test15.log

	h5diff -vc automated_test/Test3/cli_results/test16/simulation/db.hdf5 automated_test/testruns/openwarp/test16.h5 >automated_test/Test3/test16.log

	h5diff -vc automated_test/Test3/cli_results/test16a/simulation/db.hdf5 automated_test/testruns/openwarp/test16a.h5 >automated_test/Test3/test16a.log

	h5diff -vc automated_test/Test3/cli_results/test17/simulation/db.hdf5 automated_test/testruns/openwarp/test17.h5 >automated_test/Test3/test17.log

	h5diff -vc automated_test/Test3/cli_results/test17a/simulation/db.hdf5 automated_test/testruns/openwarp/test17a.h5 >automated_test/Test3/test17a.log

	h5diff -vc automated_test/Test3/cli_results/test17b/simulation/db.hdf5 automated_test/testruns/openwarp/test17b.h5 >automated_test/Test3/test17b.log

	h5diff -vc automated_test/Test3/cli_results/test17c/simulation/db.hdf5 automated_test/testruns/openwarp/test17c.h5 >automated_test/Test3/test17c.log

	h5diff -vc automated_test/Test3/cli_results/test18/simulation/db.hdf5 automated_test/testruns/openwarp/test18.h5 >automated_test/Test3/test18.log

	h5diff -vc automated_test/Test3/cli_results/test19/simulation/db.hdf5 automated_test/testruns/openwarp/test19.h5 >automated_test/Test3/test19.log

	h5diff -vc automated_test/Test3/cli_results/test20/simulation/db.hdf5 automated_test/testruns/openwarp/test20.h5 >automated_test/Test3/test20.log

	h5diff -vc automated_test/Test3/cli_results/test21/simulation/db.hdf5 automated_test/testruns/openwarp/test21.h5 >automated_test/Test3/test21.log

	h5diff -vc automated_test/Test3/cli_results/test22/simulation/db.hdf5 automated_test/testruns/openwarp/test22.h5 >automated_test/Test3/test22.log
	
	h5diff -vc automated_test/Test3/cli_results/test22a/simulation/db.hdf5 automated_test/testruns/openwarp/test22a.h5 >automated_test/Test3/test22a.log

	h5diff -vc automated_test/Test3/cli_results/test23/simulation/db.hdf5 automated_test/testruns/openwarp/test23.h5 >automated_test/Test3/test23.log

	h5diff -vc automated_test/Test3/cli_results/test24/simulation/db.hdf5 automated_test/testruns/openwarp/test24.h5 >automated_test/Test3/test24.log

	h5diff -vc automated_test/Test3/cli_results/test25/simulation/db.hdf5 automated_test/testruns/openwarp/test25.h5 >automated_test/Test3/test25.log

#	python compare.py -args	 testruns/out/test11.out
#	python compare.py -args	 testruns/out/test13.out	
echo "---------------------------finishing Test3--------------------------------"

#------------------------------------------------------------------------------------------------------------------
# Test4. Testing the Dipoles/Thin Panels Implementation in OpenWarp
#	  Added Dipoles values as [673, 960] in json
echo "---------------------------Starting  Test4--------------------------------"
 	 python openwarpgui/openwarp_cli.py automated_test/Test4/test4.json
	 
	h5diff -vc automated_test/Test4/cli_results/test21/simulation/db.hdf5 automated_test/testruns/openwarp/test21.h5 >automated_test/Test4/test21.log

	h5diff -vc automated_test/Test4/cli_results/test09/simulation/db.hdf5 automated_test/testruns/openwarp/test09.h5 >automated_test/Test4/test09.log
#	 python compare.py -args testruns/out/test9.out		 
echo "---------------------------finishing Test4--------------------------------"
#--------------------------------------------------------------------------------------------------
# Test5. Testing the Irregular Frequency Removal in OpenWarp
echo "---------------------------Starting  Test5--------------------------------"
	  python openwarpgui/openwarp_cli.py automated_test/Test5/test5.json
	 
	h5diff -vc automated_test/Test5/cli_results/test02/simulation/db.hdf5 automated_test/testruns/openwarp/test02.h5 >automated_test/Test6/test02.log

	
echo "---------------------------finishing Test5--------------------------------"
# -------------------------------------------------------------------------------------------------------------
echo "---------------------------Starting  Test6--------------------------------"
# Test6. Testing the Mean Drift forces and Yaw moment implementation in OpenWarp
# 
	  python openwarpgui/openwarp_cli.py automated_test/Test6/test6.json
	 
	h5diff -vc automated_test/Test6/cli_results/test01/simulation/db.hdf5 automated_test/testruns/openwarp/test01.h5 >automated_test/Test6/test01.log

echo "---------------------------finishing Test6--------------------------------"

