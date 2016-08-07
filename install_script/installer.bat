
@echo off

:: --------------------------------
:: Find Current Working Directory 
:: --------------------------------

SET DIR=%~dp0
ECHO Current Working Directory %DIR%

::---------------------------------------------------------------------------------------------------------------
:: Find Platform Architecture (32bit or 64bit)
::--------------------------------------------------------------------------------------------------------------

ECHO Checking Platform Architecture 
SET "flag="
if defined ProgramFiles(x86) (
    
    @echo Some 64-bit work
    SET "flag=64"
    SET CURL=%DIR%\curl\x64\curl.exe
      
) else (
    
    @echo Some 32-bit work
    SET "flag=32"
    SET	CURL=%DIR%\curl\win32\curl.exe
)

ECHO CURL path is: %CURL%

:: -----------------------------------
:: Download and installing  Python 2.7 
:: -----------------------------------
	ECHO Downloading Python 2.7 
	curl -O https://www.python.org/ftp/python/2.7.12/python-2.7.12.msi
	ECHO Installing Python quietly
	python-2.7.12.msi /quiet

:: -------------------------------------------------------------------------
:: Download and Install 7Zip (To extract other software during installation) only installing x86 version
:: -------------------------------------------------------------------------
	ECHO Downloading 7zip
	curl -O http://www.7-zip.org/a/7z1602.exe
	ECHO Installing 7zip quietly
	start /wait 7z1602.exe
::----------------------
:: Download and Extract MINGW 4.8.1
::----------------------
	ECHO Downloading MINGW 4.8.1
	curl -O https://sourceforge.net/projects/mingwbuilds/files/host-windows/releases/4.8.1/64-bit/threads-posix/sjlj/x64-4.8.1-release-posix-sjlj-rev5.7z/download
	ECHO Extracting MINGW 
	


:: -------------------------------
:: Set Path variable
:: -----------------------------
setx PATH "%CURL%" /m


rem Download the dependent software inside folder win_depedency

rem Extracting the downloaded software using 7zip 

rem How do I check if whether it is installed already ! (registry settings,)
rem How to verify the download checksum etc. 

rem Installing the dependent software inside folder win_dependency

rem Add to path 

rem use pip to install Anaconda and other libraries 