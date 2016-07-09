# Import the bemio.mesh_utilities module
import mesh
import numpy as np
import sys 

gdfFile = sys.argv[1]
# Read WAMIT mesh
buoy = mesh.read(file_name=gdfFile)

# Save to a NEMOH mesh
buoy.write(mesh_format='NEMOH')

