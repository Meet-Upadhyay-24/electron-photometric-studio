"""
01 — Build and inspect detector geometry.

Shows the three ways to get a detector (shipped mesh, generated shape, your own
Abaqus .inp file), how to place it, and the standard four-detector layout.
"""
from pathlib import Path

import matplotlib.pyplot as plt

import bse_recon as br
from bse_recon.plotting import plot_detector_2d, plot_setup

OUT = Path(__file__).parent / "output"
OUT.mkdir(exist_ok=True)

# 1. The circular mesh shipped with the package (90 triangles).
disk = br.Detector.default_disk()
print(disk)

# 2. Generated shapes — any size, any resolution.
ring = br.Detector.disk(diameter=40, n_rings=6, n_sectors=32)
plate = br.Detector.rectangle(width=30, height=20, nx=10, ny=6)
print(ring, plate, sep="\n")

# 3. Your own mesh: br.Detector.from_inp("my_detector.inp")
#    or br.Detector.from_mesh(nodes_xyz, triangles_0_based)

# Detectors are immutable; transformations return new objects and can be chained.
placed = plate.rotated(30, "y").translated(-120, 0, 90)
print("placed:", placed)

# The standard symmetric arrangement used throughout the examples:
detectors = br.four_quadrant_detectors(diameter=80, tilt=75, offset=132, height=87)
for d in detectors:
    print(d)

surface = br.SurfaceGrid.from_extent(50, 50, (-20, 20), (-20, 20))
source = br.ElectronSource(position=[0, 0, 400])

plot_detector_2d(ring, show_centroids=True).savefig(OUT / "01_ring_mesh.png", dpi=150)
plot_setup(detectors, surface, source).savefig(OUT / "01_setup.png", dpi=150)
plt.show()
