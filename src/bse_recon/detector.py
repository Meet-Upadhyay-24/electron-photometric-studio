"""Triangulated detector geometry and standard detector arrangements."""
from __future__ import annotations

from importlib import resources
from pathlib import Path
from typing import List, Optional, Union

import numpy as np

PathLike = Union[str, Path]


def _rotation_matrix(angle_deg: float, axis: str) -> np.ndarray:
    t = np.deg2rad(angle_deg)
    c, s = np.cos(t), np.sin(t)
    axis = axis.lower()
    if axis == "x":
        return np.array([[1, 0, 0], [0, c, -s], [0, s, c]])
    if axis == "y":
        return np.array([[c, 0, s], [0, 1, 0], [-s, 0, c]])
    if axis == "z":
        return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])
    raise ValueError("axis must be 'x', 'y' or 'z'.")


class Detector:
    """
    A detector surface made of triangles.

    Detectors are immutable: transformation methods return a new ``Detector``,
    so they can be chained::

        det = Detector.default_disk().scaled_to_diameter(40).rotated(25, "y").translated(100, 0, 120)

    Parameters
    ----------
    vertices:
        Array of shape ``(n_triangles, 3, 3)`` — three ``[x, y, z]`` corners per triangle.
    name:
        Optional label used in plots.
    """

    def __init__(self, vertices: np.ndarray, name: str = ""):
        v = np.array(vertices, dtype=float)
        if v.ndim != 3 or v.shape[1:] != (3, 3):
            raise ValueError(f"vertices must have shape (n, 3, 3); got {v.shape}.")
        v.setflags(write=False)
        self.vertices = v
        self.name = name

    # -- constructors -------------------------------------------------------
    @classmethod
    def from_mesh(cls, nodes: np.ndarray, triangles: np.ndarray, name: str = "") -> "Detector":
        """From a node array ``(n_nodes, 3)`` and 0-based triangle connectivity ``(n_tri, 3)``."""
        nodes = np.asarray(nodes, dtype=float)
        if nodes.shape[1] == 2:
            nodes = np.column_stack([nodes, np.zeros(len(nodes))])
        return cls(nodes[np.asarray(triangles, dtype=int)], name)

    @classmethod
    def from_inp(cls, path: PathLike, name: str = "") -> "Detector":
        """
        Read the ``*Node`` and ``*Element`` blocks of an Abaqus ``.inp`` mesh.

        Only 3-node (triangular) elements are used; any other element types are skipped.
        """
        nodes, tris, section = {}, [], None
        for raw in Path(path).read_text(encoding="utf-8", errors="ignore").splitlines():
            line = raw.strip()
            if not line or line.startswith("**"):
                continue
            if line.startswith("*"):
                low = line.lower()
                section = "node" if low.startswith("*node") else "elem" if low.startswith("*element") else None
                continue
            parts = [p for p in (s.strip() for s in line.split(",")) if p]
            if section == "node" and len(parts) >= 3:
                z = float(parts[3]) if len(parts) > 3 else 0.0
                nodes[int(parts[0])] = (float(parts[1]), float(parts[2]), z)
            elif section == "elem" and len(parts) == 4:
                tris.append(tuple(int(p) for p in parts[1:4]))
        if not tris:
            raise ValueError(f"No triangular elements found in {path}.")
        try:
            verts = np.array([[nodes[a], nodes[b], nodes[c]] for a, b, c in tris])
        except KeyError as e:
            raise ValueError(f"Element refers to undefined node {e.args[0]} in {path}.") from None
        return cls(verts, name or Path(path).stem)

    @classmethod
    def default_disk(cls) -> "Detector":
        """The circular detector mesh (90 triangles, ~5.7 mm across) shipped with the package."""
        ref = resources.files("bse_recon").joinpath("data/disk_detector.inp")
        with resources.as_file(ref) as p:
            return cls.from_inp(p, name="disk")

    @classmethod
    def disk(cls, diameter: float, n_rings: int = 4, n_sectors: int = 24, name: str = "disk") -> "Detector":
        """A flat circular detector of the given diameter in the ``z = 0`` plane, centred on the origin."""
        r = np.linspace(0, diameter / 2, n_rings + 1)
        a = np.linspace(0, 2 * np.pi, n_sectors, endpoint=False)
        pts = [(0.0, 0.0)] + [(ri * np.cos(ai), ri * np.sin(ai)) for ri in r[1:] for ai in a]
        idx = lambda ring, k: 1 + (ring - 1) * n_sectors + (k % n_sectors)  # noqa: E731
        tris = [(0, idx(1, k), idx(1, k + 1)) for k in range(n_sectors)]
        for ring in range(1, n_rings):
            for k in range(n_sectors):
                a0, a1 = idx(ring, k), idx(ring, k + 1)
                b0, b1 = idx(ring + 1, k), idx(ring + 1, k + 1)
                tris += [(a0, b0, b1), (a0, b1, a1)]
        return cls.from_mesh(np.array(pts), np.array(tris), name)

    @classmethod
    def rectangle(cls, width: float, height: float, nx: int = 8, ny: int = 8, name: str = "rectangle") -> "Detector":
        """A flat rectangular detector in the ``z = 0`` plane, centred on the origin."""
        xs, ys = np.linspace(-width / 2, width / 2, nx + 1), np.linspace(-height / 2, height / 2, ny + 1)
        X, Y = np.meshgrid(xs, ys)
        nodes = np.column_stack([X.ravel(), Y.ravel()])
        i = lambda r, c: r * (nx + 1) + c  # noqa: E731
        tris = [t for r in range(ny) for c in range(nx)
                for t in ((i(r, c), i(r, c + 1), i(r + 1, c + 1)), (i(r, c), i(r + 1, c + 1), i(r + 1, c)))]
        return cls.from_mesh(nodes, np.array(tris), name)

    # -- derived geometry ---------------------------------------------------
    def __len__(self) -> int:
        return len(self.vertices)

    @property
    def centroids(self) -> np.ndarray:
        """Triangle centroids, shape ``(n, 3)``."""
        return self.vertices.mean(axis=1)

    @property
    def center(self) -> np.ndarray:
        """Area-weighted centre of the detector."""
        return np.average(self.centroids, axis=0, weights=self.areas)

    @property
    def areas(self) -> np.ndarray:
        v = self.vertices
        return 0.5 * np.linalg.norm(np.cross(v[:, 1] - v[:, 0], v[:, 2] - v[:, 0]), axis=-1)

    @property
    def area(self) -> float:
        return float(self.areas.sum())

    # -- transformations (return new detectors) -----------------------------
    def translated(self, dx: float, dy: float, dz: float) -> "Detector":
        return Detector(self.vertices + np.array([dx, dy, dz]), self.name)

    def rotated(self, angle_deg: float, axis: str = "z", center: Optional[np.ndarray] = None) -> "Detector":
        """Rotate by ``angle_deg`` (right-hand rule) about an axis through ``center`` (default: origin)."""
        c = np.zeros(3) if center is None else np.asarray(center, dtype=float)
        R = _rotation_matrix(angle_deg, axis)
        return Detector((self.vertices - c) @ R.T + c, self.name)

    def scaled_to_diameter(self, diameter: float) -> "Detector":
        """Scale in x and y (about the bounding-box centre) so the larger in-plane span equals ``diameter``."""
        pts = self.vertices.reshape(-1, 3)
        lo, hi = pts[:, :2].min(axis=0), pts[:, :2].max(axis=0)
        span = float(np.max(hi - lo))
        if span <= 0:
            raise ValueError("Cannot scale a detector with zero in-plane extent.")
        c = np.r_[(lo + hi) / 2, 0.0]
        f = np.array([diameter / span, diameter / span, 1.0])
        return Detector((self.vertices - c) * f + c, self.name)

    def renamed(self, name: str) -> "Detector":
        return Detector(self.vertices, name)

    def __repr__(self) -> str:
        c = self.center
        return (f"Detector({self.name!r}, {len(self)} triangles, area={self.area:.4g}, "
                f"center=[{c[0]:.4g}, {c[1]:.4g}, {c[2]:.4g}])")


def four_quadrant_detectors(
    diameter: float = 40.0,
    tilt: float = 25.0,
    offset: float = 100.0,
    height: float = 120.0,
    azimuth: float = 0.0,
    mesh: Optional[Detector] = None,
) -> List[Detector]:
    """
    Four identical detectors arranged symmetrically around the beam axis.

    Each detector is scaled to ``diameter``, tilted by ``tilt`` degrees so that it
    faces the beam axis, placed ``offset`` from the axis at ``height`` above the
    ``z = 0`` plane, and finally the whole set is rotated by ``azimuth`` degrees
    about the z-axis.

    Returns
    -------
    list of Detector
        In the order ``[left (-x), top (+y), right (+x), bottom (-y)]`` (before
        the azimuthal rotation). This order is what :func:`bse_recon.contrast`
        expects.
    """
    base = (mesh if mesh is not None else Detector.default_disk()).scaled_to_diameter(diameter)
    placements = [
        ("left", "y", -tilt, (-offset, 0.0, height)),
        ("top", "x", -tilt, (0.0, offset, height)),
        ("right", "y", tilt, (offset, 0.0, height)),
        ("bottom", "x", tilt, (0.0, -offset, height)),
    ]
    return [
        base.rotated(angle, axis).translated(*shift).rotated(azimuth, "z").renamed(name)
        for name, axis, angle, shift in placements
    ]
