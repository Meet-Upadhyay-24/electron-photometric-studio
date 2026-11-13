"""
Simulate and reconstruct a surface
==================================

The complete photometric-stereo loop on a synthetic surface:

    surface → detector images → contrast → slopes → height map → error

All lengths are in millimetres.
"""

import matplotlib.pyplot as plt

import bse_recon as br
from bse_recon.plotting import plot_images, plot_reconstruction, plot_slopes

# %%
# Set up the geometry
# -------------------
# Four 80 mm detectors are tilted 75°, placed 132 mm from the beam axis and
# 87 mm above the sample. The beam pivots 1600 mm above the sample.

offset, height = 132.0, 87.0
detectors = br.four_quadrant_detectors(diameter=80, tilt=75, offset=offset, height=height)
source = br.ElectronSource(position=[0, 0, 1600])

# %%
# The sample is a 40 × 40 mm "egg-crate" surface, ``0.25·(cos 2x + cos 2y)``.

surface = br.SurfaceGrid.from_extent(150, 150, (-20, 20), (-20, 20))
surface = surface.with_heights(br.sine_surface(surface, amplitude=0.25, wavenumber=2.0))
print(surface)

# %%
# Simulate the detector images
# ----------------------------
# :func:`~bse_recon.simulate_images` returns one image per detector, shape
# ``(4, ny, nx)``. By default each image is divided by the image of a flat
# surface, which removes the shading caused by the beam moving relative to
# the detectors.

images = br.simulate_images(surface, detectors, source)
plot_images(images, surface)

# %%
# Estimate the slopes
# -------------------
# Opposed detectors give the contrasts ``Mx = (L − R)/(L + R)`` and
# ``My = (B − T)/(B + T)``. :func:`~bse_recon.estimate_slopes` converts them
# into surface slopes, correcting for the beam position ``(x, y)``.

gx, gy = br.estimate_slopes(images, surface.X, surface.Y, offset, height)
plot_slopes(surface, gx, gy)

# %%
# Integrate to a height map
# -------------------------
# :func:`~bse_recon.integrate_gradients` solves the least-squares integration
# exactly through a Sylvester equation. The absolute height offset cannot be
# observed from slopes, so the result has zero mean.

z = br.integrate_gradients(gx, gy, surface.dx, surface.dy)
print("height error [mm]:", br.height_errors(surface.heights, z))
plot_reconstruction(surface, z)
plt.show()
