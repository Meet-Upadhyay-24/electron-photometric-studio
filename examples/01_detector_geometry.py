"""
Detector geometry
=================

Detectors in ``bse_recon`` are triangle meshes. This example shows the three
ways to obtain one, how to position it, and the standard symmetric
four-detector layout used throughout the other examples.
"""

import matplotlib.pyplot as plt

import bse_recon as br
from bse_recon.plotting import plot_detector_2d, plot_setup

# %%
# Getting a detector mesh
# -----------------------
# The package ships a circular mesh of 90 triangles (about 5.6 mm across).
# It is meant to be scaled to the size you need.

disk = br.Detector.default_disk()
print(disk)

# %%
# Simple shapes can be generated at any size and resolution:

ring = br.Detector.disk(diameter=40, n_rings=6, n_sectors=32)
plate = br.Detector.rectangle(width=30, height=20, nx=10, ny=6)
print(ring)
print(plate)

plot_detector_2d(ring, show_centroids=True)

# %%
# To use your own mesh, load an Abaqus ``.inp`` file or pass node and triangle
# arrays directly::
#
#     det = br.Detector.from_inp("my_detector.inp")
#     det = br.Detector.from_mesh(nodes_xyz, triangles)   # 0-based connectivity
#
# Only the triangle vertices matter. The detector is assumed to collect every
# electron that reaches its surface.

# %%
# Placing a detector
# ------------------
# Detectors are immutable: each transformation returns a new detector, so
# calls can be chained. Rotations are in degrees about an axis through the
# origin, unless a ``center`` is given.

placed = plate.rotated(30, "y").translated(-120, 0, 90)
print(placed)

# %%
# The four-detector layout
# ------------------------
# :func:`~bse_recon.four_quadrant_detectors` scales a mesh to ``diameter``,
# tilts each copy by ``tilt`` degrees so that it faces the beam axis, and
# places the copies at distance ``offset`` from the axis and ``height`` above
# the sample. The detectors are returned in the order the inversion expects:
# ``[left, top, right, bottom]``.

detectors = br.four_quadrant_detectors(diameter=80, tilt=75, offset=132, height=87)
for d in detectors:
    print(d)

surface = br.SurfaceGrid.from_extent(50, 50, (-20, 20), (-20, 20))
source = br.ElectronSource(position=[0, 0, 400])
plot_setup(detectors, surface, source)
plt.show()
