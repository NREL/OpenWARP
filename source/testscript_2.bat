@echo off

:: --------------------------------
:: Find Current Working Directory 
:: --------------------------------

SET DIR=%~dp0
ECHO Current Working Directory %DIR%

:: Strip trailing backslash 
set DIR1=%DIR:~0,-1%

::  ~dp does not work for regular environment variables:
::  set ParentDirectory=%Directory:~dp%  set ParentDirectory=%Directory:~dp%
::  ~dp only works for batch file parameters and loop indexes
  for %%d in (%DIR1%) do set PARENTDIR=%%~dpd


ECHO Parent Directory is %PARENTDIR%

SET ROOT=%PARENTDIR%source\openwarpgui\
ECHO Root Directory is %ROOT%

IF EXIST C:\Anaconda (
	SET "ANACONDA=C:\Anaconda"
	SET "ANACONDASCRIPTS=C:\Anaconda\Scripts"
	SET "PATH=%PATH%;%ANACONDA%;%ANACONDASCRIPTS%"
	)
		
	IF EXIST %UserProfile%\Anaconda (
	SET "ANACONDA=%UserProfile%\Anaconda"
	SET "ANACONDASCRIPTS=%UserProfile%\Anaconda\Scripts"
	SET "PATH=%PATH%;%ANACONDA%;%ANACONDASCRIPTS%"
	)

SET "MINGW_ROOT=C:\mingw64\"
SET "PATH=%PATH%;%MINGW_ROOT%bin;%MINGW_ROOT%lib"


echo "###  [Test Script 2] ###"


echo "------------------------"
echo "-->Starting Test 7-"
echo "------------------------"

:: Workspace is automated_test\Test7\exp\2Bodies_true\
python openwarpgui\openwarp_cli.py automated_test\Test7\exp\2bodies_true.json
echo " --> completed test7 with USE_ODE_INFLUENCE_COEFFICIENTS as true <--"

:: Workspace is automated_test\Test7\exp\2Bodies_false\
python openwarpgui\openwarp_cli.py automated_test\Test7\exp\2bodies_false.json
echo "--> completed  test7 with USE_ODE_INFLUENCE_COEFFICIENTS as false  "

:: Diff the two generated files 
h5diff -vc automated_test\Test7\exp\2Bodies_true\db.hdf5 automated_test\Test7\exp\2Bodies_false\db.hdf5 >automated_test\Test7\exp\test07_2bodies.log &
echo "--> Comparing both the outputs .."
echo "--> H5DIFF Log at: \automated_test\Test7\test07.log"
echo "--> Test 7 finished successfully "

echo "------------------------"
echo "-->Starting Test8 "
echo "------------------------"

:: Workspace is automated_test\Test8\exp\2Bodies_true\
python openwarpgui\openwarp_cli.py automated_test\Test8\exp\2bodies_true.json
echo " --> completed test7 with USE_ODE_INFLUENCE_COEFFICIENTS as true <--"

:: Workspace is automated_test\Test8\2Bodies_false\
python openwarpgui\openwarp_cli.py automated_test\Test8\exp\2bodies_false.json
echo "--> completed  test7 with USE_ODE_INFLUENCE_COEFFICIENTS as false  "

:: Diff the two generated files 
h5diff -vc automated_test\Test8\exp\2Bodies_true\db.hdf5 automated_test\Test8\exp\2Bodies_false\db.hdf5 >automated_test\Test8\exp\test08_2bodies.log &
echo "--> Comparing both the outputs .."
echo "--> H5DIFF Log at: \automated_test\Test7\test07.log"
echo "--> Test 7\exp finished successfully "


:: Test 09
:: ---------------------------------------------------------------
echo "----------------------"
echo "--> Starting Test 9 -"
echo "----------------------"
	python openwarpgui\openwarp_cli.py automated_test\Test9\2bodies.json
echo "--> simulation workspace at : automated_test\Test9\2bodies"
echo "--> Test 9 finished successfully "


:: Test 10 
:: ---------------------------------------------------------------
echo "----------------------"
echo "--> Starting Test 10 -"
echo "----------------------"
	python openwarpgui\openwarp_cli.py automated_test\Test10\cylinder.json
echo "--> simulation workspace at : automated_test\Test10\Cylinder"
echo "--> Test 10 finished successfully"
:: Test 11 
:: ---------------------------------------------------------------
echo "----------------------"
echo "--> Starting Test 11 -"
echo "----------------------"
	python openwarpgui\openwarp_cli.py automated_test\Test11\nonsymmetrical.json
echo "--> simulation workspace at : automated_test\Test11\NonSymmetrical"
echo "--> Test 11 finished successfully"
