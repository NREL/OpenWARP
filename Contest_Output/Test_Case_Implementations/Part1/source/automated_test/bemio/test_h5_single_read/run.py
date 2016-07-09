# Load the needed modules from bemio
from bemio.io.wamit import read
from bemio.io.output import write_hdf5
from bemio.utilities.hdf_utilities import combine_h5
from bemio.utilities.hdf_utilities import create_hydro_data

# Load the data using the wamit module.
wamit_1 = read(out_file='wamit_data/coer_comp_f.out')
write_hdf5(wamit_1,out_file='wamit_1.h5')

# Load the data using the wamit module.
wamit_2 = read(out_file='wamit_data/coer_comp_f.out')
write_hdf5(wamit_2,out_file='wamit_2.h5')

hdf5_data = combine_h5(['wamit_1.h5','wamit_2.h5'])
hydro_data = create_hydro_data(hdf5_data)
