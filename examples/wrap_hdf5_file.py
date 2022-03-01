#!/bin/env python

from metadata_wrapper.hdf5 import HDF5MetadataWrapper
from metadata_wrapper.read_binary import BinaryFile

import numpy as np
import unyt
import h5py

def wrapper_test():
    
    # Create a HDF5 file with dummy data
    f = h5py.File("test.hdf5", "w")
    x = np.arange(100, dtype=np.int32)
    f["test_data"] = x
    f.close()

    # Open the input HDF5 file
    input_hdf5_file = h5py.File("test.hdf5", "r")
    
    # Create the output HDF5 wrapper file
    wrapper = HDF5MetadataWrapper("wrapped_test_hdf5.hdf5")

    # Open the dataset in the target file to be wrapped
    dataset = input_hdf5_file["test_data"]

    # Add the new dataset to the wrapper file:
    # This creates a HDF5 virtual dataset that refers to the raw data in the input HDF5 file
    # and adds attributes with the description and units following the convention used by SWIFT.
    wrapper.add_dataset("test_data", dataset, unit=unyt.km/unyt.s, description="Speed in km/s")


if __name__ == "__main__":
    wrapper_test()
