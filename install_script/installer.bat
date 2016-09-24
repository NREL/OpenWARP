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

::---------------------------------------------------------------------------------------------------------------
:: Find Platform Architecture (32bit or 64bit)
::--------------------------------------------------------------------------------------------------------------

ECHO Checking Platform Architecture 
SET "flag="
if defined ProgramFiles(x86) (
    
    @echo Some 64-bit work
    SET "flag=64"
    SET CURL=%DIR%curl\x64
      
) else (
    
    @echo Some 32-bit work
    SET "flag=32"
    SET	CURL=%DIR%curl\x86
	ECHO "32-bit is not certified yet!"
	goto :endProcess
)
ECHO OS Architecture %FLAG%
ECHO CURL path is: %CURL%
ECHO Setting curl to PATH
SET "PATH=%PATH%;%CURL%"

:: -------------------------------------------------------------------------
:: Use 7Zip, (To extract other software during installation)
:: -------------------------------------------------------------------------
	:7zip	
	ECHO "Setting Path for 7za"
	SET "PATH=%PATH%;%DIR%7za"
	goto ANACONDA

	:7ZInstalled
	ECHO 7Z is installed already 
	
:: -------------------------------------
:: Download and install Anaconda 
:: -------------------------------------	
	:ANACONDA
	
	IF EXIST C:\Anaconda (
	SET "ANACONDA=C:\Anaconda"
	SET "ANACONDASCRIPTS=C:\Anaconda\Scripts"
	SET "PATH=%PATH%;%ANACONDA%;%ANACONDASCRIPTS%"
	goto anacondaInstalled
	)
		
	IF EXIST %UserProfile%\Anaconda (
	SET "ANACONDA=%UserProfile%\Anaconda"
	SET "ANACONDASCRIPTS=%UserProfile%\Anaconda\Scripts"
	SET "PATH=%PATH%;%ANACONDA%;%ANACONDASCRIPTS%"
	goto anacondaInstalled
	)
	
	ECHO Downloading Anaconda-2.1.0-Windows-x86_64
	ECHO %DIR%Anaconda-2.1.0-Windows-x86_64.exe
	IF NOT EXIST %DIR%Anaconda-2.1.0-Windows-x86_64.exe (
	curl -O https://repo.continuum.io/archive/Anaconda-2.1.0-Windows-x86_64.exe
	
	ECHO Installing Anaconda and registering as System's Python 
	Anaconda-2.1.0-Windows-x86_64.exe /RegisterPython=1
	
	goto :ANACONDA
	)
		
	:anacondaInstalled
	ECHO "Anaconda is installed already"
	SET "PATH=%PATH%;%ANACONDA%;%ANACONDASCRIPTS%"
	:: Delete the installer 
	DEL Anaconda-2.1.0-Windows-x86_64.exe /f
	
::-----------------------------------
:: Download and Extract MINGW 4.8.1
::-----------------------------------
	:MinGw
	SET "MINGW_ROOT=C:\mingw64\"
	:: IF C:\mingw64\ exists, Assume mingw is already installed 
	IF EXIST C:\mingw64\ (goto MinGwInstalled)
	
	ECHO Downloading MINGW 4.8.1
	IF NOT EXIST %DIR%x64-4.8.1-release-posix-sjlj-rev5.7z (
	curl -L -O --retry 20 --retry-max-time 6000 https://sourceforge.net/projects/mingwbuilds/files/host-windows/releases/4.8.1/64-bit/threads-posix/sjlj/x64-4.8.1-release-posix-sjlj-rev5.7z
	
	ECHO Extracting MINGW USING 7Z
	7za x x64-4.8.1-release-posix-sjlj-rev5.7z -oC:\
	)
			
	ECHO Setting Environment Variables
	SET "PATH=%PATH%;%MINGW_ROOT%bin;%MINGW_ROOT%lib"
	:: Deleting the downloaded file to avoid corruption case in next run
	DEL x64-4.8.1-release-posix-sjlj-rev5.7z /f
	goto cmake

	:MinGwInstalled
	ECHO MinGw is installed already 
	SET "PATH=%PATH%;%MINGW_ROOT%bin;%MINGW_ROOT%lib"
	
:: -----------------------------------------------
:: Download and Install CMAKE
:: -----------------------------------------------
	:cmake		
	SET "CMAKE=C:\cmake-2.8.12.2-win32-x86\bin"
	:: IF CMake\bin exists, Assume Cmake is installed ! 
	IF EXIST C:\cmake-2.8.12.2-win32-x86\bin (goto CMAKEInstalled)

	ECHO Downloading CMAKE 2.8
	IF NOT EXIST %DIR%cmake-2.8.12.2-win32-x86.zip (
	curl -O https://cmake.org/files/v2.8/cmake-2.8.12.2-win32-x86.zip

	ECHO Extracting CMAKE
	7za x cmake-2.8.12.2-win32-x86.zip -oC:\
	)
	
	SET CMAKE=C:\cmake-2.8.12.2-win32-x86\bin
	ECHO Setting Environment Variables for CMAKE
	SET "PATH=%PATH%;%CMAKE%"
	goto gfortranBuild

	:CMAKEInstalled
	ECHO CMake is installed already
	SET "PATH=%PATH%;%CMAKE%"
	DEL cmake-2.8.12.2-win32-x86.zip /f
	ECHO PATH IS: %PATH%	

:: ----------------------------------
:: Making GFortran Build Folder (create an empty folder to create makefiles ! )
:: ----------------------------------
	:gfortranBuild
	CD %PARENTDIR%
	ECHO Creating new GFortran build folder 
	RMDIR source\NemohImproved\Nemoh\gFortranBuild /S /Q
	MKDIR source\NemohImproved\Nemoh\gFortranBuild
	CD %PARENTDIR%source\NemohImproved\Nemoh\
	goto build
		

:: -----------------------------------------------
:: Build libnemoh.dll and libnemoh.dll.a 
:: -----------------------------------------------
	:build
	cmake -DCMAKE_Fortran_COMPILER="gfortran" "%PARENTDIR%source\NemohImproved\Nemoh" -G "MinGW Makefiles"
	mingw32-make
	CD %PARENTDIR%

:: --------------------------------------------------------------------
:: Copy libnemoh.dll to Anaconda/Dlls, Anaconda/Libs and MinGW/libs	
:: --------------------------------------------------------------------
	ECHO Copying libnemoh.dll to Anaconda\Dlls
	echo %ANACONDA%\DLLs
	copy %PARENTDIR%source\NemohImproved\Nemoh\libnemoh.dll %ANACONDA%\DLLs
	copy %PARENTDIR%source\NemohImproved\Nemoh\libnemoh.dll.a %ANACONDA%\DLLs
	
	ECHO Copying libnemoh.dll to Anaconda\libs
	copy %PARENTDIR%source\NemohImproved\Nemoh\libnemoh.dll %ANACONDA%\libs
	copy %PARENTDIR%source\NemohImproved\Nemoh\libnemoh.dll.a %ANACONDA%\libs
	
	ECHO Copying libnemoh.dll to Anaconda
	copy %PARENTDIR%source\NemohImproved\Nemoh\libnemoh.dll %MINGW_ROOT%\lib
	copy %PARENTDIR%source\NemohImproved\Nemoh\libnemoh.dll.a %MINGW_ROOT%\lib
	
	
:: ---------------------------------------------------------------------------
:: Copy Blas and Lapack from install_script folder to bundled/simulations/libs
:: ---------------------------------------------------------------------------
	
	ECHO Copying BlAS and LAPACK to MinGW/libs
	copy %PARENTDIR%install_script\dlls\libblas.dll %MINGW_ROOT%lib% 
	copy %PARENTDIR%install_script\dlls\liblapack.dll %MINGW_ROOT%lib% 
	
	ECHO Copying 19 necessary Dlls inside Anaconda\Dlls
	xcopy /s %PARENTDIR%install_script\dlls %ANACONDA%\DLLs
	
	IF NOT EXIST %MINGW_ROOT%\lib\libnemoh.dll (
	ECHO libnemoh.dll was not copied in the last step so copying it from repo
	copy %PARENTDIR%install_script\dlls\libnemoh.dll %ANACONDA%\DLLs
	copy %PARENTDIR%install_script\dlls\libnemoh.dll.a %ANACONDA%\DLLs
	copy %PARENTDIR%install_script\dlls\libnemoh.dll %ANACONDA%\libs
	copy %PARENTDIR%install_script\dlls\libnemoh.dll.a %ANACONDA%\libs
	copy %PARENTDIR%install_script\dlls\libnemoh.dll %MINGW_ROOT%\lib
	copy %PARENTDIR%install_script\dlls\libnemoh.dll.a %MINGW_ROOT%\lib
	)
	
	CD %PARENTDIR%install_script
	
:: ------------------------------------------------
:: Installing VCREDIST x64 , it installs MSVSCRT10.DLL
:: ---------------------------------------------
	%PARENTDIR%install_script\vsredistributable12\vcredist_x64.exe	
	
	
:: ----------------------------------------------------------
:: Installing HDF5-tools (H5DIFF are used in testscripts)
:: ---------------------------------------------------------
	%PARENTDIR%install_script\hdf5\HDF5-1.8.17-win64.msi

:: ----------------------------------------
:: Installing Python libraries using pip
:: -----------------------------------------

	ECHO Installing required python libraries using pip
	ECHO %ROOT%
	ECHO %ANACONDA%\Scripts\pip
		
		
	%ANACONDA%\Scripts\pip install -r %ROOT%requirements.txt
	%ANACONDA%\Scripts\pip install --upgrade numpy

	::---- for Automated Test 2 (Bemio-Wamit)
	%ANACONDA%\Scripts\pip install %DIR%numpy-1.9.2+mkl-cp27-none-win_amd64.whl
	%ANACONDA%\Scripts\pip install progressbar
	
		
:: --------------------------------
:: Installing ParaView
:: ---------------------------------
	ECHO Downloading Paraview 
	ECHO Download Paraview from dropbox !
	:: Downloadig from dropbox as curl is throwing some issues with Parview sever
	IF NOT EXIST %DIR%ParaView-4.1.0-Windows-64bit.zip (
	
	curl -L -o %DIR%ParaView-4.1.0-Windows-64bit.zip https://www.dropbox.com/s/7h6zfjhm23bizjb/ParaView-4.1.0-Windows-64bit.zip?dl=1
	ECHO Installing Paraview
	7za x ParaView-4.1.0-Windows-64bit.zip -o%ROOT%\bundled\
	
	)
	xcopy /s %ROOT%\bundled\ParaView-4.1.0-Windows-64bit %ROOT%\bundled\paraview
	rmdir %ROOT%\bundled\ParaView-4.1.0-Windows-64bit /s /q

	::ParaView-4.1.0-Windows-64bit
	::SET "PATH=%PATH%;C:\Program Files (x86)\ParaView 4.1.0\bin"
	
	
:: --------------------------------
:: Building Solver Fortran 
:: -------------------------------

	ECHO Building setup.py 
	ECHO %ROOT%
	python %ROOT%\nemoh\setup.py cleanall
	python %ROOT%\nemoh\setup.py build_ext --inplace

	ECHO "OpenWarp Installation Steps are complete !"
	ECHO "Execute run.bat for GUI and testscript.bat for CLI tests !"

:: ------------------------------
:: Execute UI
:: ------------------------------
	:: ECHO "RUNNING OPENWARP"
	::CD %ROOT%
	::python main.py
	
	
:endProcess
ECHO "Ending Process !"	