"""
Detector design sweep
=====================

A template for choosing a detector layout. It sweeps tilt, height and diameter,
reconstructs a random rough surface for each design, and tabulates the
error. Edit the parameter lists to match the constraints of your machine.
"""

import itertools
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

import bse_recon as br

tilts = [0, 40, 75]           # degrees
heights = [87, 150, 272]      # mm
diameters = [20, 80]          # mm
offset = 132.0

# %%
# Test surface
# ------------
# A Gaussian random surface (RMS 0.2 mm, correlation length 2 mm). The edges
# are tapered to zero so that the boundaries do not dominate the error.

source = br.ElectronSource([0, 0, 1600])
s = br.SurfaceGrid.from_extent(80, 80, (-20, 20), (-20, 20))
s = s.with_heights(br.random_surface(s, rms=0.2, correlation_length=2.0, seed=0) * br.edge_taper(s.shape))

# %%
# Run the sweep
# -------------

rows = []
for tilt, h, dia in itertools.product(tilts, heights, diameters):
    dets = br.four_quadrant_detectors(diameter=dia, tilt=tilt, offset=offset, height=h)
    gx, gy = br.estimate_slopes(br.simulate_images(s, dets, source), s.X, s.Y, offset, h)
    err = br.height_errors(s.heights, br.integrate_gradients(gx, gy, s.dx, s.dy))
    rows.append((tilt, h, dia, err["rmse"]))
    print(f"tilt {tilt:3d} deg  height {h:4d}  diameter {dia:3d}  ->  RMSE {err['rmse']:.4f} mm")

# %%
# Compare the designs
# -------------------
# Large detectors sitting close to the sample give the largest errors.

rows = np.array(rows)
fig, axes = plt.subplots(1, len(diameters), figsize=(10, 3.8), sharey=True)
for ax, dia in zip(axes, diameters):
    for tilt in tilts:
        sel = (rows[:, 2] == dia) & (rows[:, 0] == tilt)
        ax.plot(rows[sel, 1], rows[sel, 3] * 1000, "o-", label=f"tilt {tilt}°")
    ax.set(title=f"diameter {dia} mm", xlabel="detector height [mm]")
axes[0].set_ylabel("height RMSE [µm]")
axes[0].legend()
fig.tight_layout()

# %%
# Save the table to CSV for further analysis:

out = Path("output")
out.mkdir(exist_ok=True)
np.savetxt(out / "05_sweep.csv", rows, delimiter=",", header="tilt,height,diameter,rmse", comments="")
plt.show()
