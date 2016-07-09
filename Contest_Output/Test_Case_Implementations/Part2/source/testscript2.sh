#! /bin/bash
# [Topcoder]
# This script consists of Test-case from 7 to 11 
exec >> automated_test.log
# 
# Test 07 
# ---------------------------------------------------------------
echo "------------------------"
echo "-->Starting Test 7 -"
echo "------------------------"
# Workspace is automated_test/Test7/flap-meshed-quad_true/
python openwarpgui/openwarp_cli.py automated_test/Test7/test7_true.json
echo " --> completed test7 with USE_ODE_INFLUENCE_COEFFICIENTS as true <--"
# Workspace is automated_test/Test7/flap-meshed-quad_false/
python openwarpgui/openwarp_cli.py automated_test/Test7/test7_false.json
echo "--> completed  test7 with USE_ODE_INFLUENCE_COEFFICIENTS as false  "
# Diff the two generated files 
h5diff -vc automated_test/Test7/flap-meshed-quad_true/simulation/db.hdf5 automated_test/Test7/flap-meshed-quad_false/simulation/db.hdf5 >automated_test/Test7/test07.log
echo "--> Comparing both the outputs .."
echo "--> H5DIFF Log at: /automated_test/Test7/test07.log"
echo "--> Test 7 finished successfully "



# Test 08
# ---------------------------------------------------------------
echo "----------------------"
echo "--> Starting Test 8 -"
echo "----------------------"
# Workspace is automated_test/Test7/flap-meshed-quad_negative/
python openwarpgui/openwarp_cli.py automated_test/Test8/test8_negative.json
echo "--> Testing Green function, tabulation with null values  "
# Workspace is automated_test/Test7/flap-meshed-quad_default/
python openwarpgui/openwarp_cli.py automated_test/Test8/test8_default.json
echo "--> Testing Green function, tabulation with default values "
# Diff the two generated files 
h5diff -vc automated_test/Test8/flap-meshed-quad_negative/simulation/db.hdf5 automated_test/Test8/flap-meshed-quad_default/simulation/db.hdf5 >automated_test/Test8/test08.log
echo "--> Comparing both the outputs .."
echo "--> H5DIFF Log at: /automated_test/Test8test08.log"
echo "--> Test 8 finished successfully "



# Test 09
# ---------------------------------------------------------------
echo "----------------------"
echo "--> Starting Test 9 -"
echo "----------------------"
	python openwarpgui/openwarp_cli.py automated_test/Test9/2bodies.json
echo "--> simulation workspace at : automated_test/Test9/2Bodies"
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
