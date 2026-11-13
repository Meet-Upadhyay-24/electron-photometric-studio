"""
02 — The complete photometric-stereo loop on a synthetic surface.

surface → simulated detector images → contrast → slopes → height map → error
"""
from pathlib import Path

import matplotlib.pyplot as plt

import bse_recon as br
from bse_recon.plotting import plot_images, plot_reconstruction, plot_slopes

OUT = Path(__file__).parent / "output"
OUT.mkdir(exist_ok=True)

# --- 1. geometry (all lengths in mm) ------------------------------------------
offset, height = 132.0, 87.0                       # detector distance from beam axis, height above sample
detectors = br.four_quadrant_detectors(diameter=80, tilt=75, offset=offset, height=height)
source = br.ElectronSource(position=[0, 0, 1600])  # beam pivot point

surface = br.SurfaceGrid.from_extent(150, 150, (-20, 20), (-20, 20))
surface = surface.with_heights(br.sine_surface(surface, amplitude=0.25, wavenumber=2.0))
print(surface)

# --- 2. forward model ----------------------------------------------------------
images = br.simulate_images(surface, detectors, source)   # (4, ny, nx), flat-field normalised

# --- 3. inversion ------------------------------------------------------------
gx, gy = br.estimate_slopes(images, surface.X, surface.Y, offset, height)
z = br.integrate_gradients(gx, gy, surface.dx, surface.dy)

# --- 4. assess ---------------------------------------------------------------
print("height error [mm]:", br.height_errors(surface.heights, z))

plot_images(images, surface).savefig(OUT / "02_images.png", dpi=150)
plot_slopes(surface, gx, gy).savefig(OUT / "02_slopes.png", dpi=150)
plot_reconstruction(surface, z).savefig(OUT / "02_reconstruction.png", dpi=150)
plt.show()
