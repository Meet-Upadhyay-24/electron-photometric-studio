"""
05 — Sweep detector design parameters and tabulate reconstruction error.

Use this as a template for designing a detector layout for your own machine:
edit the parameter lists, run, and inspect ``output/05_sweep.csv``.
"""
import csv
import itertools
from pathlib import Path

import bse_recon as br

OUT = Path(__file__).parent / "output"
OUT.mkdir(exist_ok=True)

tilts = [0, 40, 75]           # degrees
heights = [87, 150, 272]      # mm
diameters = [20, 80]          # mm
offset = 132.0

source = br.ElectronSource([0, 0, 1600])
s = br.SurfaceGrid.from_extent(100, 100, (-20, 20), (-20, 20))
s = s.with_heights(br.random_surface(s, rms=0.2, correlation_length=2.0, seed=0) * br.edge_taper(s.shape))

rows = []
for tilt, h, dia in itertools.product(tilts, heights, diameters):
    dets = br.four_quadrant_detectors(diameter=dia, tilt=tilt, offset=offset, height=h)
    gx, gy = br.estimate_slopes(br.simulate_images(s, dets, source), s.X, s.Y, offset, h)
    err = br.height_errors(s.heights, br.integrate_gradients(gx, gy, s.dx, s.dy))
    rows.append({"tilt": tilt, "height": h, "diameter": dia, **err})
    print(f"tilt {tilt:3d}°  height {h:4d}  diameter {dia:3d}  →  RMSE {err['rmse']:.4f} mm")

with open(OUT / "05_sweep.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=rows[0].keys())
    w.writeheader()
    w.writerows(rows)
print("wrote", OUT / "05_sweep.csv")
