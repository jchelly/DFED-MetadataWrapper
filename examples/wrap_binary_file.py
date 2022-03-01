#!/bin/env python

from metadata_wrapper.hdf5 import HDF5MetadataWrapper
from metadata_wrapper.read_binary import BinaryFile

import numpy as np
import unyt

class TestFile(BinaryFile):
    """
    Class to read a simple binary file containing a single numpy array.
    The array is assumed to be 1D with 100 int32 elements.
    """
    def __init__(self, fname, *args):
        BinaryFile.__init__(self, fname, *args)
        self.add_dataset(name="test_data", dtype=np.int32, shape=(100,))

def wrapper_test():
    
    # Create a binary file with dummy data
    f = open("test.dat", "wb")
    x = np.arange(100, dtype=np.int32)
    x.tofile(f)
    f.close()

    # Open the input binary file
    input_binary_file = TestFile("test.dat")
    
    # Create the HDF5 wrapper file
    wrapper = HDF5MetadataWrapper("wrapped_test_dat.hdf5")

    # Open the dataset in the target file to be wrapped
    dataset = input_binary_file["test_data"]

    # Add the new dataset to the wrapper file:
    # This creates a HDF5 external dataset that refers to the raw data in the binary file
    # and adds attributes with the description and units following the convention used by SWIFT.
    wrapper.add_dataset("test_data", dataset, unit=unyt.km/unyt.s, description="Speed in km/s")


if __name__ == "__main__":
    wrapper_test()
