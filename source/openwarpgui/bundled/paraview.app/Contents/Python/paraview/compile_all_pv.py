import compileall
compileall.compile_dir('/Users/kitware/Dashboards/MyTests/NightlyMaster/ParaViewSuperbuild-Release-Python27/paraview/src/paraview-build/lib/site-packages/paraview')
file = open('/Users/kitware/Dashboards/MyTests/NightlyMaster/ParaViewSuperbuild-Release-Python27/paraview/src/paraview-build/lib/site-packages/paraview/pv_compile_complete', 'w')
file.write('Done')
