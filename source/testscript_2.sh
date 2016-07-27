#! /bin/bash
# [Topcoder]
# This script consists of Test-case from 7 to 11 
exec >> automated_test2.log

# Enable set -e to stop script in case of exception 
#set -e
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

echo "###  [Test Script 2] ###"


echo "------------------------"
echo "-->Starting Test 7-"
echo "------------------------"

# Workspace is automated_test/Test7/exp/2Bodies_true/
python openwarpgui/openwarp_cli.py automated_test/Test7/exp/2bodies_true.json
echo " --> completed test7 with USE_ODE_INFLUENCE_COEFFICIENTS as true <--"

# Workspace is automated_test/Test7/exp/2Bodies_false/
python openwarpgui/openwarp_cli.py automated_test/Test7/exp/2bodies_false.json
echo "--> completed  test7 with USE_ODE_INFLUENCE_COEFFICIENTS as false  "

# Diff the two generated files 
h5diff -vc automated_test/Test7/exp/2Bodies_true/db.hdf5 automated_test/Test7/exp/2Bodies_false/db.hdf5 >automated_test/Test7/exp/test07_2bodies.log &
echo "--> Comparing both the outputs .."
echo "--> H5DIFF Log at: /automated_test/Test7/test07.log"
echo "--> Test 7 finished successfully "

echo "------------------------"
echo "-->Starting Test8 "
echo "------------------------"

# Workspace is automated_test/Test8/exp/2Bodies_true/
python openwarpgui/openwarp_cli.py automated_test/Test8/exp/2bodies_true.json
echo " --> completed test7 with USE_ODE_INFLUENCE_COEFFICIENTS as true <--"

# Workspace is automated_test/Test8/2Bodies_false/
python openwarpgui/openwarp_cli.py automated_test/Test8/exp/2bodies_false.json
echo "--> completed  test7 with USE_ODE_INFLUENCE_COEFFICIENTS as false  "

# Diff the two generated files 
h5diff -vc automated_test/Test8/exp/2Bodies_true/db.hdf5 automated_test/Test8/exp/2Bodies_false/db.hdf5 >automated_test/Test8/exp/test07_2bodies.log &
echo "--> Comparing both the outputs .."
echo "--> H5DIFF Log at: /automated_test/Test7/test07.log"
echo "--> Test 7/exp finished successfully "


# Test 09
# ---------------------------------------------------------------
echo "----------------------"
echo "--> Starting Test 9 -"
echo "----------------------"
	python openwarpgui/openwarp_cli.py automated_test/Test9/2bodies.json
echo "--> simulation workspace at : automated_test/Test9/2bodies"
echo "--> Test 9 finished successfully "


# Test 10 
# ---------------------------------------------------------------
echo "----------------------"
echo "--> Starting Test 10 -"
echo "----------------------"
	python openwarpgui/openwarp_cli.py automated_test/Test10/cylinder.json
echo "--> simulation workspace at : automated_test/Test10/Cylinder"
echo "--> Test 10 finished successfully"

# Test 11 
# ---------------------------------------------------------------
echo "----------------------"
echo "--> Starting Test 11 -"
echo "----------------------"
	python openwarpgui/openwarp_cli.py automated_test/Test11/nonsymmetrical.json
echo "--> simulation workspace at : automated_test/Test11/NonSymmetrical"
echo "--> Test 11 finished successfully"
