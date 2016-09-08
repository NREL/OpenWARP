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

ECHO "RUNNING OPENWARP"
	CD %ROOT%
	python main.py