"""
Record reference outputs for the scenarios behind the publication figures.

Provenance script — NOT part of the package. It reproduces the set-ups of
``Detector variation/randomsurface/paperfigures.ipynb`` with *exactly* the
notebook's imports (``main_arrays_parallel``, ``det_funcs``,
``surf_reconstruction``) at a reduced grid size, and stores the original code's
outputs in ``reference_paper.npz`` for ``tests/test_paper.py``.

    python tests/data/make_paper_reference.py          # write reference_paper.npz
    python tests/data/make_paper_reference.py --full   # also compare full-resolution figure metrics

Scenarios
---------
origins_*   : sine surface 0.25(cos2x+cos2y), ±20 mm windows centred at (0,0), (40,40), (80,80);
              4 detectors d=80, tilt 75°, offset 132, height 87          (notebook cells 4-10)
asym_*      : sine surface 0.5(cos2x+cos2y) at the origin; d=80, tilt 80°; right detector
              perturbed in rotation / translation / height              (cells 24-28)
diam_*      : 1-D sine 0.5cos(2x); untilted detectors of diameter 1 or 80 (cell 11)
"""
from __future__ import annotations

import contextlib
import copy
import io
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
RESEARCH = HERE.parents[2]          # .../BSD
sys.path.insert(0, str(RESEARCH))

with contextlib.redirect_stdout(io.StringIO()):
    from det_funcs import detector_from_inp  # noqa: E402
    from main_arrays_parallel import ElectronSource, SurfacePlateGrid, difference_image  # noqa: E402
    from surf_reconstruction import sylvester_reconstruction  # noqa: E402

INP = RESEARCH / "detector_geometry.inp"
SOURCE = ElectronSource(position=[0, 0, 1600], voltage=30e3, current=10e-3, dwell_time=0.4e-6, transimpedance=1e5)
N = 80


def quad(dia, vra, translate, height):
    """Notebook detector construction (cell 4 / 24 / 29)."""
    base = detector_from_inp(str(INP))
    base.scale_to_diameter(dia)
    plates = [copy.deepcopy(base) for _ in range(4)]
    lt, t, r, b = plates
    lt.rotate_nodes(-vra, "y"); r.rotate_nodes(vra, "y")      # noqa: E702
    t.rotate_nodes(-vra, "x"); b.rotate_nodes(vra, "x")      # noqa: E702
    lt.translate_nodes(-translate, 0, height); r.translate_nodes(translate, 0, height)  # noqa: E702
    t.translate_nodes(0, translate, height); b.translate_nodes(0, -translate, height)  # noqa: E702
    return [lt, t, r, b], base


def diameter_quad(dia, translate=132, height=87):
    """Cell 11: flat detectors, translated, scaled, then re-centred."""
    plates = [detector_from_inp(str(INP)) for _ in range(4)]
    lt, t, r, b = plates
    lt.translate_nodes(-translate, 0, height); r.translate_nodes(translate, 0, height)  # noqa: E702
    t.translate_nodes(0, translate, height); b.translate_nodes(0, -translate, height)  # noqa: E702
    for p in plates:
        p.scale_to_diameter(dia)
    c = [np.mean([n.position for n in p.nodes], axis=0) for p in plates]
    lt.translate_nodes(-(c[0][0] + translate), -c[0][1], c[0][2] - height)
    r.translate_nodes(-(c[2][0] - translate), -c[2][1], c[2][2] - height)
    t.translate_nodes(-c[1][0], -(c[1][1] - translate), c[1][2] - height)
    b.translate_nodes(-c[3][0], -(c[3][1] + translate), c[3][2] - height)
    return plates


def verts(plate):
    return np.array([[e.node1.position, e.node2.position, e.node3.position] for e in plate.elements])


def run(name, detectors, heights_fn, origin, translate, height, n=N):
    x0, y0 = origin
    sp = SurfacePlateGrid(n, n, x0 - 20, x0 + 20, y0 - 20, y0 + 20)
    sp.heights = heights_fn(sp.x_grid, sp.y_grid)
    sp.calc_surface_normal()
    cp = SurfacePlateGrid(n, n, x0 - 20, x0 + 20, y0 - 20, y0 + 20)
    cp.set_heights(0)
    cp.calc_surface_normal()
    with contextlib.redirect_stdout(io.StringIO()):
        imgs = difference_image(SOURCE, sp, cp, detectors, dw=translate, dh=height, return_images=True)
    zcalc, z = sylvester_reconstruction(sp)
    out = {
        "intensities": np.array(imgs, dtype=np.float32),
        "nsx": sp.nsx.astype(np.float32), "nsy": sp.nsy.astype(np.float32),
        "zcalc": zcalc.astype(np.float32),
        "vertices": np.array([verts(d) for d in detectors]),
        "setup": np.array([n, x0, y0, translate, height]),
    }
    mae = np.mean(np.abs(zcalc - z))
    print(f"{name:22s} MAE(orig)={mae:.5f} mm")
    return {f"{name}/{k}": v for k, v in out.items()}


def egg(a):
    return lambda X, Y: a * (np.sin(2 * X + np.pi / 2) + np.sin(2 * Y + np.pi / 2))


def main(full=False):
    data = {}
    dets, _ = quad(80, 75, 132, 87)
    for o in [(0, 0), (40, 40), (80, 80)]:
        data.update(run(f"origins_{o[0]}_{o[1]}", dets, egg(0.25), o, 132, 87))

    dets80, base = quad(80, 80, 132, 87)
    for param, factor in [("rotation", 0.05), ("translation", 0.05), ("height", -0.10)]:
        vra, tr, h = 80.0, 132.0, 87.0
        if param == "rotation":
            vra *= 1 + factor
        elif param == "translation":
            tr *= 1 + factor
        else:
            h *= 1 + factor
        r = copy.deepcopy(base)
        r.rotate_nodes(vra, "y")
        r.translate_nodes(tr, 0, h)
        dets_a = [dets80[0], dets80[1], r, dets80[3]]
        # the notebook passes the *perturbed* translate/height as dw/dh
        data.update(run(f"asym_{param}_{factor:+.2f}", dets_a, egg(0.5), (0, 0), tr, h))

    for dia in [1, 80]:
        data.update(run(f"diam_{dia}", diameter_quad(dia), lambda X, Y: 0.5 * np.sin(2 * X + np.pi / 2),
                        (0, 0), 132, 87))

    np.savez_compressed(HERE / "reference_paper.npz", **data)
    print("wrote", HERE / "reference_paper.npz")

    if full:
        full_resolution_check()


def full_resolution_check():
    """Notebook cell 6-9 at the published 300×300 resolution: original vs package MAE per origin."""
    sys.path.insert(0, str(HERE.parents[1] / "src"))
    import bse_recon as br

    dets, _ = quad(80, 75, 132, 87)
    pkg_dets = br.four_quadrant_detectors(diameter=80, tilt=75, offset=132, height=87)
    src = br.ElectronSource([0, 0, 1600])
    print("\nfull resolution (300×300) — reconstruction MAE [mm]")
    for o in [(0, 0), (40, 40), (80, 80)]:
        sp = SurfacePlateGrid(300, 300, o[0] - 20, o[0] + 20, o[1] - 20, o[1] + 20)
        sp.set_heights(0.25, 2)
        sp.calc_surface_normal()
        cp = SurfacePlateGrid(300, 300, o[0] - 20, o[0] + 20, o[1] - 20, o[1] + 20)
        cp.set_heights(0)
        cp.calc_surface_normal()
        with contextlib.redirect_stdout(io.StringIO()):
            difference_image(SOURCE, sp, cp, dets, dw=132, dh=87)
        zc, z = sylvester_reconstruction(sp)
        mae_nb = np.mean(np.abs(zc - z))                               # as computed in the notebook
        mae_orig = np.mean(np.abs((zc - zc.mean()) - (z - z.mean())))   # offset removed

        s = br.SurfaceGrid.from_extent(300, 300, (o[0] - 20, o[0] + 20), (o[1] - 20, o[1] + 20))
        s = s.with_heights(br.sine_surface(s, 0.25, 2.0))
        gx, gy = br.estimate_slopes(br.simulate_images(s, pkg_dets, src), s.X, s.Y, 132, 87)
        zp = br.integrate_gradients(gx, gy, s.dx, s.dy)
        mae_pkg = np.mean(np.abs(zp - (s.heights - s.heights.mean())))
        print(f"  origin {o}: notebook-style {mae_nb:.4f} (offset {zc.mean() - z.mean():+.4f})   "
              f"original {mae_orig:.4f}   package {mae_pkg:.4f}")


if __name__ == "__main__":
    main(full="--full" in sys.argv)
