- The merged Nemoh code is available at NemohMerged/

- The deployment guide doc of the merged Nemoh code is available at NemohMerged/docs/Dipoles Implementation in NEMOH.pdf
Note that there is minimal changes to this documentation as the merged nemoh code can still be
run as before independant of the openwarp gui


- The updated openwarp gui application is located at the root of this submission (The directory containing the folder src/)


- The updated deployment guide doc of the OpenWARP GUI application is available at docs/Merge Code and Update GUI Deployment Guide.pdf

- There is a complete quick installation guide for Windows available at QuickWindowsInstallStep.txt and in section 3.4 of docs/Merge Code and Update GUI Deployment Guide.pdf


- The new Nemoh code do not lib to .exe anymore. Instead it leads to .dll which are available in

src\bundled\simulation\libs

- You might need to edit src/openwarp/settings.py to change the default port from 80 to another if you run into permission issue. The settings name is WEB_SERVER_PORT


- For testing the mesh, simulation and post processing we have make sure that the default values are correct.
All you need to do is to provide the path to the mesh file when required



- Don't forget to download paraview and copy it to the relevant path if you want to test the visualize



- If you want to recompile mesh-nglib and do not have access to http://apps.topcoder.com/forums/?module=Thread&threadID=829570&start=0  you can request it.

Note that in the guide provided in above link, you should skip any instructions concerning Nemoh including any code changes.
Just execute steps for compiling mesh-nglib
