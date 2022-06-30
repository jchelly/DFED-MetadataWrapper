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
mass = data.dark_matter.masses

# Generate smoothing lengths for the dark matter
from swiftsimio.visualisation.smoothing_length_generation import generate_smoothing_lengths
data.dark_matter.smoothing_length = generate_smoothing_lengths(
    data.dark_matter.coordinates,
    data.metadata.boxsize,
    kernel_gamma=1.8,
    neighbours=57,
    speedup_fac=2,
    dimension=3,
)

# Project the dark matter mass
from swiftsimio.visualisation.projection import project_pixel_grid
dm_mass = project_pixel_grid(
    # Note here that we pass in the dark matter dataset not the whole
    # data object, to specify what particle type we wish to visualise
    data=data.dark_matter,
    boxsize=0.1*data.metadata.boxsize,
    resolution=1024,
    project="masses",
    parallel=True,
    region=None
)

from matplotlib.pyplot import imshow
from matplotlib.colors import LogNorm
imshow(LogNorm()(dm_mass).T, cmap="inferno", origin="lower",
       extent=(load_region[0][0], load_region[0][1], load_region[1][0], load_region[1][1]))
plt.xlabel("x [Mpc/h]")
plt.ylabel("y [Mpc/h]")
plt.gca().set_aspect("equal")
plt.title("Millennium simulation, z=0")
