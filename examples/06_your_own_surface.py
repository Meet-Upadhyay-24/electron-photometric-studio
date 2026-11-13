"""
06 — Simulate BSE images of your own measured surface.

Supply a height map (e.g. from a profilometer or confocal microscope) as a 2-D
array plus its pixel spacing. Without an input file this script makes a
synthetic rough surface so it runs out of the box:

    python 06_your_own_surface.py                       # synthetic
    python 06_your_own_surface.py heights.npy 0.05      # your data, 0.05 mm pixels
    python 06_your_own_surface.py heights.csv 0.05      # comma-separated also works

Heights must be in the same unit as the spacing (mm here).
"""
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

import bse_recon as br
from bse_recon.plotting import plot_images, plot_reconstruction

OUT = Path(__file__).parent / "output"
OUT.mkdir(exist_ok=True)

if len(sys.argv) >= 3:
    path, spacing = Path(sys.argv[1]), float(sys.argv[2])
    heights = np.load(path) if path.suffix == ".npy" else np.loadtxt(path, delimiter=",")
    ny, nx = heights.shape
    x = (np.arange(nx) - (nx - 1) / 2) * spacing      # centre the scan on the beam axis
    y = (np.arange(ny) - (ny - 1) / 2) * spacing
    surface = br.SurfaceGrid(x, y, heights)
else:
    surface = br.SurfaceGrid.from_extent(200, 200, (-10, 10), (-10, 10))
    h = br.random_surface(surface, rms=0.1, correlation_length=0.4, seed=42)
    surface = surface.with_heights(h * br.edge_taper(surface.shape, 0.1))
print(surface)

offset, height = 132.0, 87.0
detectors = br.four_quadrant_detectors(diameter=80, tilt=75, offset=offset, height=height)
images = br.simulate_images(surface, detectors, br.ElectronSource([0, 0, 1600]))
np.save(OUT / "06_images.npy", images)                # (4, ny, nx): left, top, right, bottom

gx, gy = br.estimate_slopes(images, surface.X, surface.Y, offset, height)
z = br.integrate_gradients(gx, gy, surface.dx, surface.dy)
print("height error [mm]:", br.height_errors(surface.heights, z))

plot_images(images, surface).savefig(OUT / "06_images.png", dpi=150)
plot_reconstruction(surface, z).savefig(OUT / "06_reconstruction.png", dpi=150)
plt.show()
