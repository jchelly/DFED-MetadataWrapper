#!/bin/env python

import sys
import re
import numpy as np

from metadata_wrapper.read_binary import BinaryFile
from metadata_wrapper.hdf5 import HDF5MetadataWrapper

class GadgetBinarySnapshotFile(BinaryFile):
    """
    Class which provides a h5py-like interface to a binary snapshot file.
    Gadget binary snapshots are Fortran unformatted binary files.

    By default this can read the following quantities:

      Coordinates
      Velocities
      Masses
      ParticleIDs
      InternalEnergy
      Density
      SmoothingLength

    Extra datasets can be specified via the "extra" parameter,
    which should be a sequence of tuples. Each tuple consists of

    (name, typestr, shape, ptypes)
    
    where the components are

    name    : name of the dataset
    typestr : string, either "float" or "int" depending on type of data.
              Number of bytes per quantity is determined from record markers.
    shape   : shape of the data for ONE particle:
              should be () for scalar quanities, (3,) for vector quantities
    ptypes  : sequence of six booleans, true for particle types which have this quantity

    These datasets should be specified in the order in which they
    appear in the file.

    """
    def __init__(self, fname, extra=None):
        BinaryFile.__init__(self, fname)

        # Read the header record marker and establish endian-ness
        irec = self.read_and_skip(np.uint32)
        if irec == 256:
            self.enable_byteswap(False)
        elif irec == 65536:
            self.enable_byteswap(True)
        else:
            raise IOError("Header record length is incorrect!")
            
        # Define header blocks
        self.add_attribute("Header/NumPart_ThisFile", np.int32,   (6,))
        self.add_attribute("Header/MassTable",        np.float64, (6,))
        self.add_attribute("Header/Time",             np.float64)
        self.add_attribute("Header/Redshift",         np.float64)
        self.add_attribute("Header/Flag_Sfr",         np.int32)
        self.add_attribute("Header/Flag_Feedback",    np.int32)
        self.add_attribute("Header/NumPart_Total",    np.uint32,  (6,))
        self.add_attribute("Header/Flag_Cooling",     np.int32)
        self.add_attribute("Header/NumFilesPerSnapshot", np.int32)
        self.add_attribute("Header/BoxSize",          np.float64)
        self.add_attribute("Header/Omega0",           np.float64)
        self.add_attribute("Header/OmegaLambda",      np.float64)
        self.add_attribute("Header/HubbleParam",      np.float64)
        self.add_attribute("Header/Flag_StellarAge",  np.int32)
        self.add_attribute("Header/Flag_Metals",      np.int32)
        self.add_attribute("Header/NumPart_Total_HighWord", np.uint32,  (6,))
        self.skip_bytes(256+4-self.offset)

        # Get total number of particles in this file
        npart_type = self["Header"].attrs["NumPart_ThisFile"][...]
        masstable  = self["Header"].attrs["MassTable"][...]
        npart      = sum(npart_type)

        # Check end of header marker
        irec = self.read_and_skip(np.uint32)
        if irec != 256:
            raise IOError("Header end of record marker is incorrect!")

        # Make full list of datasets to read
        all_datasets = (
            ("Coordinates",     "float", (3,), (True,)*6),
            ("Velocities",      "float", (3,), (True,)*6),
            ("ParticleIDs",     "int",   (),   (True,)*6),
            ("Masses",          "float", (),   masstable==0),
            ("InternalEnergy",  "float", (),   (True, False, False, False, False, False)),
            ("Density",         "float", (),   (True, False, False, False, False, False)),
            ("SmoothingLength", "float", (),   (True, False, False, False, False, False)),
            )
        # Add any user specified fields
        if extra is not None:
            all_datasets += extra
        
        # Determine what datasets are present in this file
        count_records = 0
        for(name, typestr, shape, ptypes) in all_datasets:

            # Calculate number of particles we expect in this dataset
            nextra = sum(npart_type[np.asarray(ptypes, dtype=bool)])
            if nextra > 0:

                # Calculate number of numbers per particle
                n_per_part = 1
                for s in shape:
                    n_per_part *= s

                # Check if there's another record in the file
                try:
                    irec = self.read_and_skip(np.uint32)
                except IOError:
                    if count_records < 3:
                        # pos, vel, ids should always be present
                        raise
                    else:
                        break
                else:
                    count_records += 1

                # Determine bytes per quantitiy
                nbytes = np.int64(irec) // (n_per_part*nextra)
                if (nbytes != 4 and nbytes != 8) or nbytes*n_per_part*nextra != irec:
                    raise IOError("%s record has unexpected length!" % name)

                # Determine data type for this record
                if typestr == "int":
                    if nbytes==4:
                        dtype = np.int32
                    else:
                        dtype = np.int64
                elif typestr == "float":
                    if nbytes==4:
                        dtype = np.float32
                    else:
                        dtype = np.float64
                else:
                    raise ValueError("typestr parameter should be 'int' or 'float'")

                # Loop over particle types and add datasets
                for i in range(6):
                    if ptypes[i] and npart_type[i] > 0:
                        full_shape = (npart_type[i],)+tuple(shape)
                        self.add_dataset("PartType%i/%s" % (i, name), dtype, full_shape)

                # Read end of record marker
                irec = self.read_and_skip(np.uint32)
                if irec != n_per_part * np.dtype(dtype).itemsize * nextra:
                    raise IOError("%s end of record marker is incorrect!" % name)

                        
def wrap_gadget_snapshot(input_file_name, output_file_name):
    """
    Make a HDF5 wrapper file with metadata for a Gadget binary snapshot
    """
    
    # Open the input file
    infile = GadgetBinarySnapshotFile(input_file_name)

    # Create the output file
    outfile = HDF5MetadataWrapper(output_file_name)

    # Copy the Gadget snapshot header
    group = outfile.file.create_group("Header")
    for name, value in infile["Header"].attrs.items():
        group.attrs[name] = value

    # Find the number of particles of each type
    npart = infile["Header"].attrs["NumPart_ThisFile"]

    # Loop over particle types
    for ptype in range(len(npart)):
        if npart[ptype] > 0:
            
            # Ensure the group exists in the output file
            group_name = "PartType%d" % ptype
            outfile.file.create_group(group_name)

            # Loop over datasets we have for this particle type
            for dataset_name in infile[group_name]:

                # Wrap this dataset
                dataset = infile[group_name][dataset_name]
                outfile.add_binary_dataset(group_name+"/"+dataset_name, dataset=dataset)


if __name__ == "__main__":

    if len(sys.argv) == 3:
        input_file_name = sys.argv[1]
        output_file_name = sys.argv[2]
    else:
        print("Usage: python3 ./wrap_gadget_snapshot.py input_gadget_file output_wrapper_file")
        sys.exit()

    wrap_gadget_snapshot(input_file_name, output_file_name)
