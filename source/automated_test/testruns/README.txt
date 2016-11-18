This directory contains the input files described in Appendix A of the user manual.
It also contains the following files:

config.wam	sample config.wam file

runtests.bat	sample batch file for running all tests within a directory

runwamit.bat	sample batch file for running a single test within a directory


##############################

To run either set of tests, perform the following from a DOS prompt

1) edit config.wam as appropriate for your system
2) To execute all of the tests simply execute runtests.bat
	Note that the batch file must be modified to the appropriate path for wamit
	if \wamitv7 was not used at installation
3) To execute a specific test execute runwamit.bat (testname)
	Note that the batch file must be modified to the appropriate path for wamit
	if \wamitv7 was not used at installation
	For example:  runwamit.bat test12
	will execute test12
 