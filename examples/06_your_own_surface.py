"""
Simulate your own surface
=========================

Use a measured height map, for example from a profilometer or confocal
microscope, as the sample. Pass the file and its pixel spacing on the command
line::

    python 06_your_own_surface.py heights.npy 0.05      # 0.05 mm pixels
    python 06_your_own_surface.py heights.csv 0.05      # comma-separated text

Heights must be in the same unit as the spacing (mm here). Without arguments,
a synthetic rough surface is used.
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

import bse_recon as br
from bse_recon.plotting import plot_images, plot_reconstruction

# %%
# Load or generate the height map
# -------------------------------
# A measured map only needs centred coordinates, built from its pixel spacing.
# The beam axis passes through the centre of the scan.

if len(sys.argv) >= 3:
    path, spacing = Path(sys.argv[1]), float(sys.argv[2])
    heights = np.load(path) if path.suffix == ".npy" else np.loadtxt(path, delimiter=",")
    ny, nx = heights.shape
    x = (np.arange(nx) - (nx - 1) / 2) * spacing
    y = (np.arange(ny) - (ny - 1) / 2) * spacing
    surface = br.SurfaceGrid(x, y, heights)
else:
    surface = br.SurfaceGrid.from_extent(200, 200, (-10, 10), (-10, 10))
    h = br.random_surface(surface, rms=0.1, correlation_length=0.4, seed=42)
    surface = surface.with_heights(h * br.edge_taper(surface.shape, 0.1))
print(surface)

# %%
# Simulate the images
# -------------------
# The image stack can be saved and processed later like real data. The order
# is left, top, right, bottom.

offset, height = 132.0, 87.0
detectors = br.four_quadrant_detectors(diameter=80, tilt=75, offset=offset, height=height)
images = br.simulate_images(surface, detectors, br.ElectronSource([0, 0, 1600]))

out = Path("output")
out.mkdir(exist_ok=True)
np.save(out / "06_images.npy", images)
plot_images(images, surface)

# %%
# Reconstruct and compare
# -----------------------

gx, gy = br.estimate_slopes(images, surface.X, surface.Y, offset, height)
z = br.integrate_gradients(gx, gy, surface.dx, surface.dy)
print("height error [mm]:", br.height_errors(surface.heights, z))
plot_reconstruction(surface, z)
plt.show()
