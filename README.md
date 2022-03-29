# Metadata Wrapper utility

This package contains a python module which can be used to generate HDF5
"wrapper" files which contain metadata describing external files. The target
files may be in HDF5 or binary format.

This can be used to add metadata to simulation outputs which don't contain
any metadata without modifying the original outputs. It also allows access
to binary file formats using HDF5. This is particularly helpful for older
simulation outputs which were written as Fortran unformatted files, for
example: once wrapped, the simulation files can easily be read in python using
the h5py module.

## Installation

To install the module to your home directory:
```
python3 ./setup.py install --user
```
It should then be possible to import the module with
```
import metadata_wrapper
```

## Wrapping HDF5 files

Wrapping HDF5 files is straightforward. First, create an empty wrapper file:
```
from metadata_wrapper.hdf5 import HDF5MetadataWrapper
wrapper = HDF5MetadataWrapper("wrapped_test_hdf5.hdf5")
```
Then, open the target HDF5 file to be wrapped:
```
input_hdf5_file = h5py.File("test.hdf5", "r")
```
Then, for each dataset in the target file call the add_dataset method:
```
wrapper.add_dataset(name, dataset, unit, description, extra_attributes=None)
```
This will create a virtual dataset in the wrapper file which refers to the data
in the target file. Attributes are added with the units and a description of
the dataset. Here, `name` is the name the dataset will be given in the wrapper
file, `dataset` is the h5py.Dataset object from the target file and
`description` is a description to add to the dataset in the wrapper file.

The units of the dataset are specifed by passing a unyt unit object as the unit
parameter. See the [unyt documentation](https://unyt.readthedocs.io/) for details.
The extra_attributes parameter can be used to attach additional attributes to
the dataset in the wrapper file.

Any other data which should be added to the wrapper can be written using the
wrapper.file attribute, which is a reference to underlying h5py.File object.

Once all datasets have been wrapped, the wrapper file should be closed.
```
wrapper.close()
```
It should then be possible to read the target file by opening the wrapper file
with any tool that can read HDF5.

See examples/wrap_hdf5_file.py for a simple example of wrapping a HDF5 file.

## Wrapping binary files

Wrapping binary files is more complicated because it is first necessary to
establish how to read the file. This is done using the class
`metadata_wrapper.read_binary.BinaryFile`, which provides a h5py-like interface
to binary files.

The file is assumed to consist of a series of data blocks which can be
described by specifying their name, data type and shape. To do this, create a
new class which inherits from BinaryFile and provide an __init__ method which
defines the data blocks in the file by calling the add_dataset method for each
block in the order the blocks appear in the file.

For example, to wrap a binary file that just contains single array of integers
with an integer header specifying the array size:
```
class TestFile(BinaryFile):
    def __init__(self, fname, *args):
        BinaryFile.__init__(self, fname, *args)
        # Define the block with the array length
        self.add_dataset(name="array_length", dtype=np.int32, shape=(,))
        # Now that the block has been defined, we can read it
        N = self["array_length"][()]
        # Define the block with the array data
        self.add_dataset(name="array_data",   dtype=np.int32, shape=(N,))
```
The parameters to add_dataset are:
  * name: name of the dataset
  * dtype: numpy data type corresponding to the data type in the file
  * shape: dimensions of the array in the file, or (,) for scalars

Once the new class has been defined it can be used to read the binary file in
the same way as a h5py.File object. For example, to read the 100 element array
above:
```
testfile = TestFile("test.hdf5")
test_data = testfile["test_data"][...]
```

We can then create a wrapper file as before:
```
from metadata_wrapper.hdf5 import HDF5MetadataWrapper
wrapper = HDF5MetadataWrapper("wrapped_binary.hdf5")
input_binary_file = TestFile("test.dat", "r")
```
and call the add_dataset method to wrap individual datasets:
```
wrapper.add_dataset(name, dataset, unit, description, extra_attributes=None)
```
where in this case the `dataset` parameter is an instance of the class derived
from BinaryFile.

See examples/wrap_binary_file.py for a simple example of wrapping a binary file.

## Examples

### Wrapping a binary Gadget-2 snapshot

The file examples/wrap_gadget_snapshot.py contains a more complicated example
which takes a binary simulation snapshot from the Gadget-2 n-body simulation
code and wraps it so that it can be read using HDF5 tools and adds attributes
specifying the units in the file.
