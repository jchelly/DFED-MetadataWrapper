#!/bin/env python

import numpy as np
import unyt
import peano
import h5py

from metadata_wrapper.read_binary import BinaryFile
from metadata_wrapper.hdf5 import HDF5MetadataWrapper

class SnapshotFile(BinaryFile):
    """
    Class for reading Millennium-1 snapshot files
    """
    def __init__(self, fname, *args):
        BinaryFile.__init__(self, fname, *args)

        # Header
        self.start_fortran_record(auto_byteswap=True)
        self.add_attribute("Header/NumPart_ThisFile",  np.uint32,  (6,))
        self.add_attribute("Header/MassTable",         np.float64, (6,))
        self.add_attribute("Header/Time",              np.float64)
        self.add_attribute("Header/Redshift",          np.float64)
        self.add_attribute("Header/Flag_Sfr",          np.int32)
        self.add_attribute("Header/Flag_Feedback",     np.int32)
        self.add_attribute("Header/NumPart_Total",     np.uint32,  (6,))
        self.add_attribute("Header/Flag_Cooling",      np.int32)
        self.add_attribute("Header/NumFilesPerSnapshot", np.int32)
        self.add_attribute("Header/BoxSize",         np.float64)
        self.add_attribute("Header/Omega0",          np.float64)
        self.add_attribute("Header/OmegaLambda",     np.float64)
        self.add_attribute("Header/HubbleParam",     np.float64)
        self.add_attribute("Header/Flag_StellarAge", np.int32)
        self.add_attribute("Header/Flag_Metals",     np.int32)
        self.add_attribute("Header/HashTabSize",     np.int32)
        self.add_attribute("Header/fill",            np.uint8, (84,))
        self.end_fortran_record()

        # Get number of particles in this file
        n = self["Header"].attrs["NumPart_ThisFile"][1]

        # Coordinates
        self.start_fortran_record()
        self.add_dataset("PartType1/Coordinates", np.float32, (n,3))
        self.end_fortran_record()

        # Velocities
        self.start_fortran_record()
        self.add_dataset("PartType1/Velocities", np.float32, (n,3))
        self.end_fortran_record()

         # IDs
        self.start_fortran_record()
        self.add_dataset("PartType1/ParticleIDs", np.uint64, (n,))
        self.end_fortran_record()

        # Range of hash cells in this file
        self.start_fortran_record()
        self.add_dataset("first_hash_cell", np.int32)
        self.add_dataset("last_hash_cell",  np.int32)
        self.end_fortran_record()
        
        # Calculate how many hash cells we have in this file
        nhash = self["last_hash_cell"][...] - self["first_hash_cell"][...] + 1

        # Location of first particle in each cell relative to start of file
        self.start_fortran_record()
        self.add_dataset("blockid", np.int32, (nhash,))
        self.end_fortran_record()


def wrap_millennium_snapshot(input_basename, output_basename):
    """
    Make a set of HDF5 files which wrap the specified Millennium
    simulation snapshot.

    This takes fortran unformatted binary outputs from L-Gadget2 and makes
    them compatible with the swiftsimio package or other tools that can read
    SWIFT snapshots.

    input_basename: name of the input snapshot files, minus the trailing .N
    output_basename: name of the output snapshot files, minus the trailing .N
    """

    hashtable = []

    # Create a wrapper file with references to the particle data for each snapshot
    nr_files = 1
    file_nr = 0
    while file_nr < nr_files:
        
        # Open the snapshot file
        input_filename = f"{input_basename}.{file_nr}"
        infile = SnapshotFile(input_filename)

        # Extract information we need from the header
        header_data = {}
        for name in infile["Header"].attrs:
            header_data[name] = infile["Header"].attrs[name]
        a = header_data["Time"]        # Cosmological expansion factor
        h = header_data["HubbleParam"] # Hubble parameter
        z = header_data["Redshift"]
        boxsize = header_data["BoxSize"]
        nr_files = header_data["NumFilesPerSnapshot"]
        print(f"Wrapping file {file_nr} of {nr_files}")

        # Determine hash table size
        hashtabsize = infile["Header"].attrs["HashTabSize"]
        hashbits = 1
        while hashtabsize > 2:
            hashtabsize //= 2
            hashbits += 1
        hashbits //= 3 # bits per dimension

        # Read the hash table information from this file
        first_hash_cell = infile["first_hash_cell"][()]
        last_hash_cell  = infile["last_hash_cell"][()]
        blockid         = infile["blockid"][...]
        nr_parts        = infile["Header"].attrs["NumPart_ThisFile"][1]
        hashtable.append((first_hash_cell, last_hash_cell, blockid, nr_parts))

        # Create the wrapper file
        output_filename = f"{output_basename}.{file_nr}.hdf5"
        wrapper = HDF5MetadataWrapper(output_filename)
        wrapper.file.create_group("PartType1")

        # Wrap the particle coordinates
        a_exponent =  1.0 # Positions are in comoving coordinates
        h_exponent = -1.0 # Positions are in 1/h units
        wrapper.add_dataset("PartType1/Coordinates", infile["PartType1/Coordinates"], unit=unyt.Mpc,
                            cosmological_factors=(a, a_exponent, h, h_exponent))

        # Wrap the particle velocities
        a_exponent = 0.5 # Velocities are (peculiar velocity / sqrt(a))
        h_exponent = 0.0
        wrapper.add_dataset("PartType1/Velocities", infile["PartType1/Velocities"], unit=unyt.km/unyt.s,
                            cosmological_factors=(a, a_exponent, h, h_exponent))

        # Wrap the particle IDs
        a_exponent = 0.0
        h_exponent = 0.0
        wrapper.add_dataset("PartType1/ParticleIDs", infile["PartType1/ParticleIDs"], unit=unyt.dimensionless,
                            cosmological_factors=(a, a_exponent, h, h_exponent))

        # Store data types
        pos_dtype = infile["PartType1/Coordinates"].dtype
        vel_dtype = infile["PartType1/Velocities"].dtype
        ids_dtype = infile["PartType1/ParticleIDs"].dtype

        # Advance to the next file
        wrapper.close()
        infile.close()
        file_nr += 1

    # In L-Gadget2 particles are stored in order of which 'hash cell' they belong to.
    # Will try to translate this into the top level cell grid which SWIFT uses.

    # Type to store information about a SWIFT cell
    swift_cell_t = np.dtype([
        ("centre", np.float64, 3), # coordinates of cell centre
        ("count",  np.int64),      # number of particles in the cell
        ("offset", np.int64),      # offset to first particle
        ("file",   np.int32),      # file containing this cell
        ("order",  np.int32),      # ordering of the cells in the snapshot file(s)
    ])

    # Create the array of cells
    total_nr_cells = np.sum([len(h[2]) for h in hashtable])
    assert total_nr_cells == (2**hashbits)**3
    cell = np.ndarray(total_nr_cells, dtype=swift_cell_t)

    # Determine size of a cell
    cells_per_dimension = 2**hashbits
    cell_width = boxsize / cells_per_dimension

    # Loop over files and populate array of cells
    for file_nr in range(nr_files):

        # Find the Gadget hashtable info for this file
        first_hash_cell, last_hash_cell, blockid, nr_parts = hashtable[file_nr]
        nr_cells = last_hash_cell - first_hash_cell + 1

        # Find the corresponding SWIFT cells
        cells_in_file = cell[first_hash_cell:last_hash_cell+1]
        cells_in_file["file"] = file_nr
        cells_in_file["count"][:-1] = blockid[1:] - blockid[:-1]
        cells_in_file["count"][-1] = nr_parts - blockid[-1]
        assert np.sum(cells_in_file["count"]) == nr_parts
        cells_in_file["offset"] = blockid
        
        # Compute cell centres
        keys = np.arange(first_hash_cell, last_hash_cell+1, dtype=int)
        ix, iy, iz = peano.peano_hilbert_key_inverses(keys, hashbits)
        cells_in_file["centre"][:,0] = (ix+0.5)*cell_width
        cells_in_file["centre"][:,1] = (iy+0.5)*cell_width
        cells_in_file["centre"][:,2] = (iz+0.5)*cell_width

    total_nr_parts = sum([int(h[3]) for h in hashtable])
 
    # Now go back and write metadata to the output files.
    # If we're only going to use the single 'virtual' file then
    # we only need to add the metadata to the first file.
    for file_nr in (0,): #range(nr_files):

        print(f"Writing metadata to wrapper file {file_nr}")
        output_filename = f"{output_basename}.{file_nr}.hdf5"
        outfile = h5py.File(output_filename, "r+")

        # Write the cells
        outfile.create_group("Cells")
        outfile["Cells"]["Centres"] = cell["centre"]
        outfile.create_group("Cells/Counts")
        outfile["Cells"]["Counts"]["PartType1"] = cell["count"]
        outfile.create_group("Cells/Files")
        outfile["Cells"]["Files"]["PartType1"] = cell["file"]
        outfile.create_group("Cells/OffsetsInFile")
        outfile["Cells"]["OffsetsInFile"]["PartType1"] = cell["offset"]
        outfile.create_group("Cells/Meta-data")
        outfile["Cells"]["Meta-data"].attrs["nr_cells"] = (total_nr_cells,)
        outfile["Cells"]["Meta-data"].attrs["dimension"] = (cells_per_dimension,)*3
        outfile["Cells"]["Meta-data"].attrs["size"] = (cell_width,)*3

        # Add particle number etc
        nr_parts = hashtable[file_nr][3]
        outfile["PartType1"].attrs["NumberOfParticles"] = (nr_parts,)
        outfile["PartType1"].attrs["NumberOfFields"] = (3,)
        
        # Add system of units
        units = outfile.create_group("Units")
        units.attrs["Unit length in cgs (U_L)"] = (unyt.Mpc.get_conversion_factor(unyt.cm)[0],)
        units.attrs["Unit mass in cgs (U_M)"]   = (float(unyt.solar_mass.to(unyt.g).value),)
        units.attrs["Unit time in cgs (U_t)"]   = ((unyt.Mpc/unyt.km*unyt.s).get_conversion_factor(unyt.s)[0],)
        units.attrs["Unit temperature in cgs (U_T)"] = (1.0,)
        units.attrs["Unit current in cgs (U_I)"] = (1.0,)

        # Add the header
        header = outfile.create_group("Header")
        header.attrs["BoxSize"] = (boxsize,)*3
        header.attrs["Code"] = "L-Gadget2"
        header.attrs["Dimension"] = (3,)
        header.attrs["Flag_Entropy_ICs"] = np.zeros(6, dtype=int)
        header.attrs["MassTable"] = np.zeros(6, dtype=float)
        header.attrs["NumFilesPerSnapshot"] = nr_files
        header.attrs["NumPartTypes"] = (6,)
        numpart_thisfile = np.zeros(6, dtype=np.int32)
        numpart_thisfile[1] = nr_parts
        header.attrs["NumPart_ThisFile"] = numpart_thisfile
        numpart_total_hw = np.zeros(6, dtype=np.uint32)
        numpart_total_hw[1] = (total_nr_parts >> 32)
        header.attrs["NumPart_Total_HighWord"] = numpart_total_hw
        numpart_total = np.zeros(6, dtype=np.uint32)
        numpart_total[1] = total_nr_parts - ((total_nr_parts >> 32) << 32)
        header.attrs["NumPart_Total"] = numpart_total
        header.attrs["Redshift"] = (z,)
        header.attrs["Scale-factor"] = (a,)
        header.attrs["ThisFile"] = (file_nr,)

        # Add cosmology
        cosmo = outfile.create_group("Cosmology")
        cosmo.attrs["Cosmological run"] = (1,)
        cosmo.attrs["H0"] = (100.0*header_data["HubbleParam"],)
        cosmo.attrs["H0 [internal units]"] = (100.0*header_data["HubbleParam"],)
        cosmo.attrs["h"]  = (header_data["HubbleParam"],)
        cosmo.attrs["Omega_b"] = (0.0,)
        cosmo.attrs["Omega_cdm"] = (header_data["Omega0"],)
        cosmo.attrs["Omega_m"] = (header_data["Omega0"],)
        cosmo.attrs["Omega_g"] = (0.0,)
        cosmo.attrs["Omega_k"] = (0.0,)
        cosmo.attrs["Omega_r"] = (0.0,)
        cosmo.attrs["Omega_lambda"] = (header_data["OmegaLambda"],)
        cosmo.attrs["Omega_nu"] = (0.0,)
        cosmo.attrs["Omega_nu_0"] = (0.0,)
        cosmo.attrs["Redshift"] = (z,)
        cosmo.attrs["Scale-factor"] = (a,)
        cosmo.attrs["w"]   = ( 1.0,)
        cosmo.attrs["w_0"] = (-1.0,)
        cosmo.attrs["w_a"] = ( 0.0,)
        outfile.close()

    # Reopen the first file
    source_filename = f"{output_basename}.0.hdf5"
    sourcefile = h5py.File(source_filename,"r")

    #
    # Since swiftsimio can't handle multi-file snapshots, here we create
    # a single file containing virtual datasets which concatenate the
    # datasets from all files.
    #
    output_filename = f"{output_basename}.hdf5"
    outfile = h5py.File(output_filename,"w")

    # Copy metadata
    sourcefile.copy("Header",    outfile)
    sourcefile.copy("Cosmology", outfile)
    sourcefile.copy("Units",     outfile)
    sourcefile.copy("Cells",     outfile)
    # Record dataset attributes
    pos_attrs = dict(sourcefile["PartType1/Coordinates"].attrs)
    vel_attrs = dict(sourcefile["PartType1/Velocities"].attrs)
    ids_attrs = dict(sourcefile["PartType1/ParticleIDs"].attrs)
    sourcefile.close()

    # Make updates where necessary
    outfile["Header"].attrs["NumFilesPerSnapshot"] = 1
    numpart_thisfile = np.zeros(6, dtype=np.int64)
    numpart_thisfile[1] = total_nr_parts
    outfile["Header"].attrs["NumPart_ThisFile"] = numpart_thisfile

    # Need to update cells:
    # All cells are now in file 0 and offsets need to be recomputed
    outfile["Cells"]["Files"]["PartType1"][...] = 0
    count = outfile["Cells"]["Counts"]["PartType1"][...]
    offset = np.cumsum(count) - count
    outfile["Cells"]["OffsetsInFile"]["PartType1"][...] = offset

    # Create dataset layouts
    pos_layout = h5py.VirtualLayout(shape=(total_nr_parts,3), dtype=pos_dtype)
    vel_layout = h5py.VirtualLayout(shape=(total_nr_parts,3), dtype=vel_dtype)
    ids_layout = h5py.VirtualLayout(shape=(total_nr_parts,),  dtype=ids_dtype)
    offset = 0
    for file_nr in range(nr_files):
        nr_parts = hashtable[file_nr][3]
        source_filename = f"{output_basename}.{file_nr}.hdf5"
        pos_layout[offset:offset+nr_parts,:] = h5py.VirtualSource(source_filename, "PartType1/Coordinates", shape=(nr_parts,3))
        vel_layout[offset:offset+nr_parts,:] = h5py.VirtualSource(source_filename, "PartType1/Velocities", shape=(nr_parts,3))
        ids_layout[offset:offset+nr_parts]   = h5py.VirtualSource(source_filename, "PartType1/ParticleIDs", shape=(nr_parts,))
        offset += nr_parts

    # Create datasets
    group = outfile.create_group("PartType1")
    pos_dset = group.create_virtual_dataset("Coordinates", pos_layout)
    for name in pos_attrs:
        pos_dset.attrs[name] = pos_attrs[name]
    vel_dset = group.create_virtual_dataset("Velocities",  vel_layout)
    for name in vel_attrs:
        vel_dset.attrs[name] = vel_attrs[name]
    ids_dset = group.create_virtual_dataset("ParticleIDs", ids_layout)
    for name in ids_attrs:
        ids_dset.attrs[name] = ids_attrs[name]
    outfile.close()


if __name__ == "__main__":
    
    # Full Millennium simulation
    input_basename = "/cosma6/data/dp004/Millennium/Full/snapdir_063/snap_millennium_063"
    output_basename = "/cosma6/data/dp004/jch/wrap_millennium/snap_millennium_063"

    wrap_millennium_snapshot(input_basename, output_basename)
    
