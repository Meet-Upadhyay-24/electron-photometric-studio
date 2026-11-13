"""
03 — How reconstruction error grows as the scan window moves off-axis.

This reproduces the "origins" study of the accompanying publication: a
40 × 40 mm egg-crate surface is scanned with its centre at (0, 0), (40, 40) and
(80, 80) mm, using four 80 mm detectors tilted 75° at d_w = 132 mm, d_h = 87 mm.

Pass a grid size on the command line (default 150; the paper used 300):

    python 03_beam_position_study.py 300
"""
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

import bse_recon as br

OUT = Path(__file__).parent / "output"
OUT.mkdir(exist_ok=True)
n = int(sys.argv[1]) if len(sys.argv) > 1 else 150

dw, dh = 132.0, 87.0
detectors = br.four_quadrant_detectors(diameter=80, tilt=75, offset=dw, height=dh)
source = br.ElectronSource([0, 0, 1600])

fig, axes = plt.subplots(1, 3, figsize=(12, 4), sharey=False)
for ax, (cx, cy) in zip(axes, [(0, 0), (40, 40), (80, 80)]):
    s = br.SurfaceGrid.from_extent(n, n, (cx - 20, cx + 20), (cy - 20, cy + 20))
    s = s.with_heights(br.sine_surface(s, 0.25, 2.0))
    gx, gy = br.estimate_slopes(br.simulate_images(s, detectors, source), s.X, s.Y, dw, dh)
    z = br.integrate_gradients(gx, gy, s.dx, s.dy)
    err = br.height_errors(s.heights, z)
    print(f"origin ({cx:2d},{cy:2d}) mm:  MAE {err['mae']:.4f}  RMSE {err['rmse']:.4f}  max {err['max']:.4f} mm")

    diff = (z - z.mean()) - (s.heights - s.heights.mean())
    im = ax.imshow(diff, origin="lower", cmap="viridis",
                   extent=np.array([s.x[0], s.x[-1], s.y[0], s.y[-1]]) / dw)
    ax.set(title=f"origin ({cx}, {cy}) — MAE {err['mae']:.3f} mm", xlabel="x / d_w", ylabel="y / d_w")
    fig.colorbar(im, ax=ax, shrink=0.8, label="height error [mm]")

fig.tight_layout()
fig.savefig(OUT / "03_beam_position.png", dpi=150)
plt.show()
