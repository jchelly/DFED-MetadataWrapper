#!/bin/env python
#
# Read in and plot part of the Millennium simulation
# from the HDF5 wrapped snapshot using swiftsimio.
#

import matplotlib.pyplot as plt
import swiftsimio as sw

mask = sw.mask("snap_millennium_063.hdf5")

# Choose region to read (in terms of the boxsize to get the units right)
b = mask.metadata.boxsize[0]
load_region=((0*b,.1*b),(0*b,.1*b),(0*b,.1*b))

# Seelct the region
mask.constrain_spatial(load_region)
data = sw.load("snap_millennium_063.hdf5", mask=mask)

# Read particle positions
pos = data.dark_matter.coordinates

# Make a dotplot
plt.plot(pos[:,0], pos[:,1], "k,", rasterized=True, alpha=0.05)
plt.gca().set_aspect("equal")
plt.xlabel("x [Mpc/h]")
plt.ylabel("y [Mpc/h]")
plt.title("Millennium simulation z=0")
