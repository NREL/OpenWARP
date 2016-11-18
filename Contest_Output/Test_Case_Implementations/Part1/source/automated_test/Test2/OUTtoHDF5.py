# Load the needed modules from bemio
from bemio.io.wamit import read
from bemio.io.output import write_hdf5
from bemio.utilities.hdf_utilities import combine_h5
from bemio.utilities.hdf_utilities import create_hydro_data

import sys 
import os

# .out files
inFile = sys.argv[1]
inDir = os.path.dirname(inFile)
head, tail = os.path.split(inFile)
head = os.path.dirname(head)
# H5 files 
#outFile = os.path.splitext(inFile)[0]

# Load the data using the wamit module.

wamit_1 = read(out_file=inFile+'.out')

write_hdf5(wamit_1,out_file = head + '/openwarp/' + tail + '.h5' )


