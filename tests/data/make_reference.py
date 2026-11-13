"""
Generate regression reference data from the ORIGINAL research code.

This script is kept for provenance only. It imports the original
``main_arrays_parallel.py`` from the parent research folder (not shipped with
the package) and records its outputs for a small, fixed scenario. The package's
test-suite compares against the resulting ``reference_sine.npz``.

Run from the repository root (needs the original research folder next to it):

    python tests/data/make_reference.py
"""
from __future__ import annotations

import sys
from copy import deepcopy
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
RESEARCH = REPO.parent  # folder that holds main_arrays_parallel.py

sys.path.insert(0, str(RESEARCH))
import main_arrays_parallel as orig  # noqa: E402


# --- minimal detector geometry, mirroring the research workflow --------------
class _Node:
    def __init__(self, num, xyz):
        self.node_number = num
        self.position = np.asarray(xyz, dtype=float)


class _Element:
    def __init__(self, n1, n2, n3):
        self.node1, self.node2, self.node3 = n1, n2, n3
        self.centroid = (n1.position + n2.position + n3.position) / 3.0


class _Plate:
    def __init__(self, elements):
        self.elements = elements


def _read_inp(path):
    nodes, tris, section = {}, [], None
    for raw in Path(path).read_text().splitlines():
        line = raw.strip()
        if line.startswith("*"):
            low = line.lower()
            section = "node" if low.startswith("*node") else "elem" if low.startswith("*element") else None
            continue
        parts = [p for p in (s.strip() for s in line.split(",")) if p]
        if section == "node" and len(parts) >= 3:
            z = float(parts[3]) if len(parts) > 3 else 0.0
            nodes[int(parts[0])] = np.array([float(parts[1]), float(parts[2]), z])
        elif section == "elem" and len(parts) >= 4:
            tris.append([int(p) for p in parts[1:4]])
    return nodes, tris


def _rot(axis, deg):
    t = np.deg2rad(deg)
    c, s = np.cos(t), np.sin(t)
    return {
        "x": np.array([[1, 0, 0], [0, c, -s], [0, s, c]]),
        "y": np.array([[c, 0, s], [0, 1, 0], [-s, 0, c]]),
    }[axis]


def four_detectors(inp, diameter, tilt, offset, height):
    nodes, tris = _read_inp(inp)
    pts = np.array(list(nodes.values()))
    cx, cy = (pts[:, 0].max() + pts[:, 0].min()) / 2, (pts[:, 1].max() + pts[:, 1].min()) / 2
    span = max(np.ptp(pts[:, 0]), np.ptp(pts[:, 1]))
    f = diameter / span
    base = {k: np.array([cx + (p[0] - cx) * f, cy + (p[1] - cy) * f, p[2]]) for k, p in nodes.items()}
    specs = [("y", -tilt, (-offset, 0, height)),   # left
             ("x", -tilt, (0, offset, height)),    # top
             ("y", tilt, (offset, 0, height)),     # right
             ("x", tilt, (0, -offset, height))]    # bottom
    plates, verts = [], []
    for axis, ang, tr in specs:
        R = _rot(axis, ang)
        moved = {k: R @ p + np.asarray(tr, float) for k, p in deepcopy(base).items()}
        nd = {k: _Node(k, p) for k, p in moved.items()}
        plates.append(_Plate([_Element(nd[a], nd[b], nd[c]) for a, b, c in tris]))
        verts.append(np.array([[moved[a], moved[b], moved[c]] for a, b, c in tris]))
    return plates, np.array(verts)


def main():
    n = 40
    amp, freq = 0.25, 1.0
    diameter, tilt, offset, height = 40.0, 25.0, 100.0, 120.0

    surface = orig.SurfacePlateGrid(n, n, -10, 10, -10, 10)
    surface.heights = amp * (np.sin(surface.x_grid * freq + np.pi / 2) + np.sin(surface.y_grid * freq + np.pi / 2))
    surface.calc_surface_normal()
    flat = orig.SurfacePlateGrid(n, n, -10, 10, -10, 10)
    flat.set_heights(0)
    flat.calc_surface_normal()

    source = orig.ElectronSource([0, 0, 1600], 60e3, 10e-3, 0.4e-6, 1e5)
    detectors, verts = four_detectors(REPO / "src/bse_recon/data/disk_detector.inp", diameter, tilt, offset, height)

    images = orig.difference_image(source, surface, flat, detectors, dw=offset, dh=height, return_images=True)
    raw = [orig.element_intensity_no_correction(d, flat, surface, source) for d in detectors]

    np.savez_compressed(
        HERE / "reference_sine.npz",
        params=np.array([n, amp, freq, diameter, tilt, offset, height]),
        heights=surface.heights,
        detector_vertices=verts,
        intensities=np.array(images),
        raw_intensities=np.array(raw),
        Mx=surface.Mx, My=surface.My, nsx=surface.nsx, nsy=surface.nsy,
    )
    print("wrote", HERE / "reference_sine.npz")


if __name__ == "__main__":
    main()
