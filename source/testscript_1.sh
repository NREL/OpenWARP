#! /bin/bash
# [Topcoder]
# This script consists of Test-case from 7 to 11 
exec >> automated_test.log

# to stop script when shell script returns zero
# set -e

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

echo "Starting Test-Script 1 "

  echo "-------------------------starting Test1--------------------------------"	 

# These will create a new file named filename_bemio_output in the test run folder   
  echo "Using Bemio"
python automated_test/Test1/GDFtoDAT.py 'automated_test/testruns/test01.gdf'
python automated_test/Test1/GDFtoDAT.py 'automated_test/testruns/test02.gdf'
python automated_test/Test1/GDFtoDAT.py 'automated_test/testruns/test09.gdf'
python automated_test/Test1/GDFtoDAT.py 'automated_test/testruns/test11.gdf'

# It has errors,
#python automated_test/Test1/GDFtoDAT.py 'automated_test/testruns/test13ac.gdf'


echo "--------------------------finishing Test1 ---------------------------------"

#--------------------------------------------------------------------------------------------------------------------
echo "--------------------------starting Test2 -------------------------------"

	chmod u+x automated_test/testruns/openwarp/
	#Provide filename, not extension '.out'
	# Create .H5 files from all Wamit's .out file 
	python automated_test/Test2/OUTtoHDF5.py 'automated_test/testruns/out/test01'
	python automated_test/Test2/OUTtoHDF5.py 'automated_test/testruns/out/test02'
	python automated_test/Test2/OUTtoHDF5.py 'automated_test/testruns/out/test09'
	python automated_test/Test2/OUTtoHDF5.py 'automated_test/testruns/out/test11'
	python automated_test/Test2/OUTtoHDF5.py 'automated_test/testruns/out/test13'
	
	
	# Test01
	echo "testing test01"
	python openwarpgui/openwarp_cli.py automated_test/Test2/test01.json
	h5diff -vc automated_test/Test2/test01/test01/db.hdf5 automated_test/testruns/openwarp/test01.h5 >automated_test/Test2/test01.log
	echo "test01 worked fine"
	
				
	
echo "---------------------------finishing Test2--------------------------------"
   
#---------------------------------------------------------------------------------------------------------------------
# Test3. read test3.json to that enables USE_HIGHER_ORDER ad TRUE with B_SPLINE_ORDER as 1. Later, it converts 
echo "---------------------------Starting  Test3--------------------------------"
# Running Test 11
	python openwarpgui/openwarp_cli.py automated_test/Test3/test11.json
	
	h5diff -vc automated_test/Test3/test11/test11/db.hdf5 automated_test/testruns/openwarp/test11.h5 >automated_test/Test3/test11.log

# Running Test 13
	python openwarpgui/openwarp_cli.py automated_test/Test3/test13.json

	h5diff -vc automated_test/Test3/test13/test13/db.hdf5 automated_test/testruns/openwarp/test13.h5 >automated_test/Test3/test13.log

	
#	python compare.py -args	 testruns/out/test11.out
#	python compare.py -args	 testruns/out/test13.out	
echo "---------------------------finishing Test3--------------------------------"

#------------------------------------------------------------------------------------------------------------------
# Test4. Testing the Dipoles/Thin Panels Implementation in OpenWarp
#	  Added Dipoles values as [673, 960] in json
echo "---------------------------Starting  Test4--------------------------------"
 	 python openwarpgui/openwarp_cli.py automated_test/Test4/test4.json
	 
	h5diff -vc automated_test/Test4/cli_results/test09/db.hdf5 automated_test/testruns/openwarp/test09.h5 >automated_test/Test4/test09.log
#	 python compare.py -args testruns/out/test9.out		 
echo "---------------------------finishing Test4--------------------------------"
#--------------------------------------------------------------------------------------------------
# Test5. Testing the Irregular Frequency Removal in OpenWarp
echo "---------------------------Starting  Test5--------------------------------"
	  python openwarpgui/openwarp_cli.py automated_test/Test5/test5.json
	 
	h5diff -vc automated_test/Test5/cli_results/test02/db.hdf5 automated_test/testruns/openwarp/test02.h5 >automated_test/Test5/test02.log

	
echo "---------------------------finishing Test5--------------------------------"
# -------------------------------------------------------------------------------------------------------------
echo "---------------------------Starting  Test6--------------------------------"
# Test6. Testing the Mean Drift forces and Yaw moment implementation in OpenWarp
# 
	  python openwarpgui/openwarp_cli.py automated_test/Test6/test6.json
	 
	h5diff -vc automated_test/Test6/cli_results/test01/db.hdf5 automated_test/testruns/openwarp/test01.h5 >automated_test/Test6/test01.log

echo "---------------------------finishing Test6--------------------------------"

