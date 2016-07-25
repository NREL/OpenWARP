#! /bin/bash
# [Topcoder]
# This script consists of Test-case from 7 to 11 
exec >> automated_test2.log

exec >> automated_test2.log

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

echo "Starting Test Script 2 "



# 
# Test 07 
# ---------------------------------------------------------------
echo "------------------------"
echo "-->Starting Test 7 -"
echo "------------------------"
# Workspace is automated_test/Test7/flap-meshed-quad_true/
#python openwarpgui/openwarp_cli.py automated_test/Test7/test7_true.json
echo " --> completed test7 with USE_ODE_INFLUENCE_COEFFICIENTS as true <--"
# Workspace is automated_test/Test7/flap-meshed-quad_false/
#python openwarpgui/openwarp_cli.py automated_test/Test7/test7_false.json
echo "--> completed  test7 with USE_ODE_INFLUENCE_COEFFICIENTS as false  "
# Diff the two generated files 
#h5diff -vc automated_test/Test7/flap-meshed-quad_true/db.hdf5 automated_test/Test7/flap-meshed-quad_false/db.hdf5 >automated_test/Test7/test07.log &
echo "--> Comparing both the outputs .."
echo "--> H5DIFF Log at: /automated_test/Test7/test07.log"
echo "--> Test 7 finished successfully "



# Test 08
# ---------------------------------------------------------------
echo "----------------------"
echo "--> Starting Test 8 -"
echo "----------------------"
# Workspace is automated_test/Test7/flap-meshed-quad_negative/
#python openwarpgui/openwarp_cli.py automated_test/Test8/test8_negative.json
echo "--> Testing Green function, tabulation with null values  "
# Workspace is automated_test/Test7/flap-meshed-quad_default/
#python openwarpgui/openwarp_cli.py automated_test/Test8/test8_default.json
echo "--> Testing Green function, tabulation with default values "
# Diff the two generated files 
#h5diff -vc automated_test/Test8/flap-meshed-quad_negative/db.hdf5 automated_test/Test8/flap-meshed-quad_default/db.hdf5 >automated_test/Test8/test08.log &
echo "--> Comparing both the outputs .."
echo "--> H5DIFF Log at: /automated_test/Test8/test08.log"
echo "--> Test 8 finished successfully "


# Test 09
# ---------------------------------------------------------------
echo "----------------------"
echo "--> Starting Test 9 -"
echo "----------------------"
	#python openwarpgui/openwarp_cli.py automated_test/Test9/2bodies.json
echo "--> simulation workspace at : automated_test/Test9/2bodies"
echo "--> Test 9 finished successfully "


# Test 10 
# ---------------------------------------------------------------
echo "----------------------"
echo "--> Starting Test 10 -"
echo "----------------------"
	#python openwarpgui/openwarp_cli.py automated_test/Test10/cylinder.json
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
