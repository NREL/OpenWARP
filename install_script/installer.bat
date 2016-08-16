
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
)
ECHO OS Architecture %FLAG%
ECHO CURL path is: %CURL%
ECHO Setting curl to PATH
SET PATH=%PATH%;%CURL%


:: -------------------------------------------------------------------------
:: Download and Install 7Zip (To extract other software during installation) only installing x86 version
:: -------------------------------------------------------------------------
	:7zip
	IF EXIST C:\Program Files\7-zip\ (goto :7ZInstalled)	
	:: if 7zip is installed already, don't download it.
	
	ECHO If not installed, use 7za
	SET 7ZA=%DIR%7za
	SET PATH=%PATH%;%7ZA%
	goto MinGw

	:7ZInstalled
	ECHO 7Z is installed already 
::-----------------------------------
:: Download and Extract MINGW 4.8.1
::-----------------------------------
	:MinGw
	SET MINGW_ROOT=C:\mingw64\
	IF EXIST C:\mingw64\ (goto MinGwInstalled)
	
	ECHO Downloading MINGW 4.8.1
	curl -O https://sourceforge.net/projects/mingwbuilds/files/host-windows/releases/4.8.1/64-bit/threads-posix/sjlj/x64-4.8.1-release-posix-sjlj-rev5.7z/download
	
	ECHO Extracting MINGW 
	7z x x64-4.8.1-release-posix-sjlj-rev5.7z -oC:\
	
	ECHO Setting Environment Variables
	SET PATH=%PATH%;%MINGW_ROOT%bin;%MINGW_ROOT%lib
	goto cmake

	:MinGwInstalled
	ECHO MinGw is installed already 
	ECHO Setting Environment Variables
	SET PATH=%PATH%;%MINGW_ROOT%bin;%MINGW_ROOT%lib

:: ---------------------------------------------------------------------------
:: Copy Blas and Lapack from install_script folder to bundled/simulations/libs
:: ---------------------------------------------------------------------------
	
	ECHO Copying BlAS and LAPACK
	
	:: Is this step required ? 
	IF %flag%==32 ( copy %DIR%\blas\win32\libblas.dll %ROOT%bundled\simulation\libs )
	IF %flag%==32 (copy %DIR%\lapack\win32\liblapack.dll %ROOT%bundled\simulation\libs )
	
	:: copy to MinGw Folder
	IF %flag%==32 copy %DIR%\blas\win32\libblas.dll %MINGW_ROOT%lib
	IF %flag%==32 copy %DIR%\lapack\win32\liblapack.dll %MINGW_ROOT%lib

	:: Is this step required ? 
	IF %flag%==64 copy %DIR%\blas\x64\libblas.dll %ROOT%bundled\simulation\libs
	IF %flag%==64 copy %DIR%\lapack\x64\liblapack.dll %ROOT%bundled\simulation\libs
	
	:: copy to MinGw Folder 
	IF %flag%==64 copy %DIR%\blas\x64\libblas.dll %MINGW_ROOT%lib
	IF %flag%==64 copy %DIR%\lapack\x64\liblapack.dll %MINGW_ROOT%lib
:: -----------------------------------------------
:: Download and Install CMAKE
:: -----------------------------------------------
	:cmake
	::cmake --version 2>NUL
	::if errorlevel 0 goto CMAKEInstalled
	
	SET CMAKE=C:\cmake-2.8.12.2-win32-x86\bin
	IF EXIST C:\cmake-2.8.12.2-win32-x86\ (goto CMAKEInstalled)

	ECHO Downloading CMAKE 2.8
	curl -O https://cmake.org/files/v2.8/cmake-2.8.12.2-win32-x86.zip

	ECHO Extracting CMAKE
	7z x cmake-2.8.12.2-win32-x86.zip -oC:\
	
	SET CMAKE=C:\cmake-2.8.12.2-win32-x86\bin
	ECHO Setting Environment Variables for CMAKE
	SET PATH=%PATH%;%CMAKE%
	goto gfortranBuild

	:CMAKEInstalled
	ECHO CMake is installed already
	SET PATH=%PATH%;%CMAKE%
	ECHO PATH IS %PATH%

:: ----------------------------------
:: Making GFortran Build Folder
:: ----------------------------------
	:gfortranBuild
	CD %PARENTDIR%
	ECHO Creating new GFortran build folder 
	RMDIR source\NemohImproved\Nemoh\gFortranBuildWindows /S /Q
	MKDIR source\NemohImproved\Nemoh\gFortranBuildWindows
	CD %PARENTDIR%source\NemohImproved\Nemoh\gFortranBuildWindows
	::CD source
	goto build
	:: delete this folder on new installation
:: -----------------------------------------------
:: Build libnemoh.dll and libnemoh.dll.a 
:: -----------------------------------------------
	:build
	cmake -DCMAKE_Fortran_COMPILER="gfortran" "%PARENTDIR%source\NemohImproved\Nemoh" -G "MinGW Makefiles"
	::cmake -DCMAKE_Fortran_COMPILER="gfortran" "%PARENTDIR%source/NemohImproved/Nemoh/" -G "MinGW MakeFiles"
	::cmake -DCMAKE_Fortran_COMPILER="gfortran" "%PARENTDIR%source/NemohImproved/Nemoh/" 
	mingw32-make
	
	CD %DIR%
:: ------------------------------------------------------------------------
:: Copy libnemoh.dll and libnemoh.dll.a from simulation\libs to MinGw\lib
:: -----------------------------------------------------------------------
	ECHO "Copying libnemoh.dll"
	copy %ROOT%bundled\simulation\libs\libnemoh.dll %MINGW_ROOT%lib
	copy %ROOT%bundled\simulation\libs\libnemoh.dll %MINGW_ROOT%lib
	copy %ROOT%bundled\simulation\libs\libnemoh.dll %MINGW_ROOT%lib
	copy %ROOT%bundled\simulation\libs\libnemoh.dll %MINGW_ROOT%lib

:: Download and Install Anaconda 2
:: ----------------------------------------------------
	
	ECHO Downloading Anaconda
	SET ANACONDA=C:\Anaconda2
	IF EXIST C:\Anaconda2 goto :AnacondaInstalled
	
	::curl -O https://repo.continuum.io/archive/.winzip/Anaconda2-4.1.1-Windows-x86_64.zip
	ECHO Extracting Anaconda
	7z x %DIR%Anaconda2-4.1.1-Windows-x86_64.zip -oC:\
	
	ECHO Installing Anaconda
	Anaconda2-4.1.1-Windows-x86_64.exe
	
	ECHO Setting Env for Anaconda
	SET ANACONDA=C:\Anaconda2

	:AnacondaInstalled
	ECHO "Anaconda Installed already "
	SET PATH=%PATH%;%ANACONDA%;%ANACONDA%\Scripts
:: ----------------------------------------
:: Install CherryPy and other requirements
:: -----------------------------------------

	ECHO Instally CherryPy and other requirements
	
	%ANACONDA%\Scripts\pip install -r %ROOT%requirements.txt --upgrade

	:: Remove this later on, I was getting vrvasall.bat missing 
	:: pip install –upgrade setuptools "<-- this command is not working "
	:: pip install setuptools --no-use-wheel --upgrade (<-- not helping )

	:: psycopyg windows (need to have Git installed for it ! )
	::pip install git+https://github.com/nwcell/psycopg2-windows.git@win64-py27#egg=psycopg2
:: --------------------------------
:: Install ParaView
:: ---------------------------------



:: --------------------------------
:: Build Setup.py
:: -------------------------------

	ECHO Building setup.py 
	python %ROOT%\nemoh\setup.py cleanall
	python %ROOT%\nemoh\setup.py build_ext --inplace

:: -------------------------------
:: Set Path variable
:: -----------------------------

:: Finally SET the path so that it can be used later !
::SETX PATH "%PATH%;%ANACONDA%;%ANACONDA%\Scripts /m

ECHO "OpenWarp Installation Steps are complete, Check the logs for possible errors ! "


:: ------------------------------
:: Execute UI
:: ------------------------------
	ECHO "RUNNING OPENWARP"
	CD %ROOT%
	python main.py