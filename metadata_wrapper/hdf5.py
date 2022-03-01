#!/bin/env python

import unyt
import h5py

import metadata_wrapper.swift_units
from metadata_wrapper.read_binary import BinaryDataset

class HDF5MetadataWrapper:
    """
    Class to create a HDF5 'wrapper' file which describes external
    datasets with additional metadata attached.
    """
    def __init__(self, filename, truncate=True):

        # Create or open the output HDF5 file
        mode = "w" if truncate else "r+"
        self.file = h5py.File(filename, mode)

    def add_dataset(self, name, dataset, description=None, unit=None,
                    extra_attributes=None):
        """
        Add a dataset which describes HDF5 or binary data in another file.
        The dataset parameter should be an instance of read_binary.BinaryDataset
        or h5py.Dataset.
        """
        
        if isinstance(dataset, BinaryDataset):
            # Create a HDF5 external dataset describing the input binary dataset
            dataset = self.file.create_dataset(name, shape=dataset.shape, dtype=dataset.dtype,
                                               external=((dataset.fname, dataset.offset, dataset.nbytes),))
        elif isinstance(dataset, h5py.Dataset):
            # Create a HDF5 virtual dataset which gets its data from the input HDF5 dataset
            layout = h5py.VirtualLayout(shape=dataset.shape, dtype=dataset.dtype)
            layout[...] = h5py.VirtualSource(dataset)
            dataset = self.file.create_virtual_dataset(name, layout)
        else:
            raise ValueError("Input dataset must be a metadata_wrapper.read_binary.BinaryDataset or a h5py.Dataset")

        # Add description, if we have one
        if description is not None:
            dataset.attrs["Description"] = description

        # Add units, if specified
        if unit is not None:
            metadata_wrapper.swift_units.write_unit_attributes(dataset, unit)
        
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
