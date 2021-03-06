#!/bin/env python

import unyt
import h5py

import metadata_wrapper.swift_units
from metadata_wrapper.read_binary import BinaryDataset, big_or_little

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
                    extra_attributes=None, cosmological_factors=None):
        """
        Add a dataset which describes HDF5 or binary data in another file.
        The dataset parameter should be an instance of read_binary.BinaryDataset
        or h5py.Dataset.

        For cosmology simulations, cosmological_factors should be a tuple
        (a, a_exponent, h, h_exponent) where a is the expansion factor, h is
        the Hubble parameter, and the exponents indicate the factors of
        a and h which are included in the units.
        """
        
        if isinstance(dataset, BinaryDataset):
            # Need to be careful in case of non-native endian input data
            dtype = dataset.dtype
            if dtype.byteorder != "|":
                data_endian = big_or_little(dataset.endian)
                dtype = dataset.dtype.newbyteorder(data_endian)
            # Create a HDF5 external dataset describing the input binary dataset
            new_dataset = self.file.create_dataset(name, shape=dataset.shape, dtype=dtype,
                                                   external=((dataset.fname, dataset.offset, dataset.nbytes),))
        elif isinstance(dataset, h5py.Dataset):
            # Create a HDF5 virtual dataset which gets its data from the input HDF5 dataset
            layout = h5py.VirtualLayout(shape=dataset.shape, dtype=dataset.dtype)
            layout[...] = h5py.VirtualSource(dataset)
            new_dataset = self.file.create_virtual_dataset(name, layout)
        else:
            raise ValueError("Input dataset must be a metadata_wrapper.read_binary.BinaryDataset or a h5py.Dataset")

        # Add description, if we have one
        if description is not None:
            new_dataset.attrs["Description"] = description

        # Add units, if specified
        if unit is not None:
            metadata_wrapper.swift_units.write_unit_attributes(new_dataset, unit, cosmological_factors)
        
        # Add any additional attributes
        if extra_attributes is not None:
            for name, value in extra_attributes.iter():
                new_dataset.attrs[name] = value

    def close(self):
        if self.file is not None:
            self.file.close()
            self.file = None

    def __del__(self):
        self.close()
