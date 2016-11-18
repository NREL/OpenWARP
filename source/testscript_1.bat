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



echo "Starting Test-Script 1 "

  echo "-------------------------starting Test1--------------------------------"	 

:: These will create a new file named filename_bemio_output in the test run folder   
  echo "Using Bemio"
	python automated_test\Test1\GDFtoDAT.py automated_test\testruns\test01.gdf
	python automated_test\Test1\GDFtoDAT.py automated_test\testruns\test02.gdf
	python automated_test\Test1\GDFtoDAT.py automated_test\testruns\test09.gdf
	python automated_test\Test1\GDFtoDAT.py automated_test\testruns\test11.gdf


:: It has errors,
::python automated_test\Test1\GDFtoDAT.py 'automat\d_test\testruns\test13ac.gdf'


echo "--------------------------finishing Test1 ---------------------------------"

::--------------------------------------------------------------------------------------------------------------------
echo "--------------------------starting Test2 -------------------------------"

:: -- updating numpy+mtk to resolve numpy errors
	::%ANACONDA%\Scripts\pip uninstall numpy
	::%ANACONDA%\Scripts\pip uninstall scipy
	::%ANACONDA%\Scripts\pip install --upgrade pip 
	::%ANACONDA%\Scripts\pip install %DIR%numpy-1.9.2+mkl-cp27-none-win_amd64.whl


	:: Provide filename, not extension '.out'
	:: Create .H5 files from all Wamit's .out file

	python automated_test\Test2\OUTtoHDF5.py automated_test\testruns\out\test01
	python automated_test\Test2\OUTtoHDF5.py automated_test\testruns\out\test02
	python automated_test\Test2\OUTtoHDF5.py automated_test\testruns\out\test09
	python automated_test\Test2\OUTtoHDF5.py automated_test\testruns\out\test11
	python automated_test\Test2\OUTtoHDF5.py automated_test\testruns\out\test13
	
	
	:: Test01
	echo "testing test01"
	python openwarpgui\openwarp_cli.py automated_test\Test2\test01.json
	h5diff -vc automated_test\Test2\test01\test01\db.hdf5 automated_test\testruns\openwarp\test01.h5 >automated_test\Test2\test01.log
	echo "test01 worked fine"
	
				
	
echo "---------------------------finishing Test2--------------------------------"
   
::---------------------------------------------------------------------------------------------------------------------
:: Test3. read test3.json to that enables USE_HIGHER_ORDER ad TRUE with B_SPLINE_ORDER as 1. Later, it converts 
echo "---------------------------Starting  Test3--------------------------------"
:: Running Test 11
	python openwarpgui\openwarp_cli.py automated_test\Test3\test11.json
	
	h5diff -vc automated_test\Test3\test11\test11\db.hdf5 automated_test\testruns\openwarp\test11.h5 >automated_test\Test3\test11.log

:: Running Test 13
	python openwarpgui\openwarp_cli.py automated_test\Test3\test13.json

	h5diff -vc automated_test\Test3\test13\test13\db.hdf5 automated_test\testruns\openwarp\test13.h5 >automated_test\Test3\test13.log

	
::	python compare.py -args	 testruns\out\test11.out
::	python compare.py -args	 testruns\out\test13.out	
echo "---------------------------finishing Test3--------------------------------"

::------------------------------------------------------------------------------------------------------------------
:: Test4. Testing the Dipoles\Thin Panels Implementation in OpenWarp
::	  Added Dipoles values as [673, 960] in json
echo "---------------------------Starting  Test4--------------------------------"
 	 python openwarpgui\openwarp_cli.py automated_test\Test4\test4.json
	 
	h5diff -vc automated_test\Test4\cli_results\test09\db.hdf5 automated_test\testruns\openwarp\test09.h5 >automated_test\Test4\test09.log
::	 python compare.py -args testruns\out\test9.out		 
echo "---------------------------finishing Test4--------------------------------"
::--------------------------------------------------------------------------------------------------
:: Test5. Testing the Irregular Frequency Removal in OpenWarp
echo "---------------------------Starting  Test5--------------------------------"
	  python openwarpgui\openwarp_cli.py automated_test\Test5\test5.json
	 
	h5diff -vc automated_test\Test5\cli_results\test02\db.hdf5 automated_test\testruns\openwarp\test02.h5 >automated_test\Test5\test02.log

	
echo "---------------------------finishing Test5--------------------------------"
:: -------------------------------------------------------------------------------------------------------------
echo "---------------------------Starting  Test6--------------------------------"
:: Test6. Testing the Mean Drift forces and Yaw moment implementation in OpenWarp
:: 
	  python openwarpgui\openwarp_cli.py automated_test\Test6\test6.json
	 
	h5diff -vc automated_test\Test6\cli_results\test01\db.hdf5 automated_test\testruns\openwarp\test01.h5 >automated_test\Test6\test01.log

echo "---------------------------finishing Test6--------------------------------"

