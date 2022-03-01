#!/bin/env python

import unyt
import swift_units

import h5py
import read_binary

class HDF5MetadataWrapper:
    """
    Class to create a HDF5 'wrapper' file which describes external
    datasets with additional metadata attached.
    """
    def __init__(self, filename, truncate=True):

        # Create or open the output HDF5 file
        mode = "w" if truncate else "r+"
        self.file = h5py.File(filename, mode)

    def add_binary_dataset(self, name, dataset, description=None, unit=None,
                           extra_attributes=None):
        """
        Add a dataset which describes binary data in another file.
        The dataset parameter should be an instance of the
        read_binary.BinaryDataset class.
        """
        
        # Create a HDF5 external dataset describing the input binary dataset
        dataset = self.file.create_dataset(name, shape=dataset.shape, dtype=dataset.dtype,
                                           external=((dataset.fname, dataset.offset, dataset.nbytes),))
        
        # Add description, if we have one
        if description is not None:
            dataset.attrs["Description"] = description

        # Add units, if specified
        swift_units.write_unit_attributes(dataset, unit)
        
        # Add any additional attributes
        if extra_attributes is not None:
            for name, value in extra_attributes.iter():
                dataset.attrs[name] = value

    def close(self):
        if self.file is not None:
            self.file.close()
            self.file = None

    def __del__(self):
        self.close()
