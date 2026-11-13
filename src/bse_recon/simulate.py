"""Forward model: backscattered-electron signal collected by triangulated detectors."""
from __future__ import annotations

from typing import Callable, Optional, Sequence, Union

import numpy as np

from .detector import Detector
from .materials import TI64, Material
from .scattering import lambertian
from .source import ElectronSource
from .surface import SurfaceGrid

ScatteringModel = Callable[..., np.ndarray]


def _unit(v: np.ndarray) -> np.ndarray:
    return v / np.linalg.norm(v, axis=-1, keepdims=True)


def triangle_solid_angles(points: np.ndarray, vertices: np.ndarray) -> np.ndarray:
    """
    Solid angle subtended by each triangle at each point (Van Oosterom & Strackee, 1983).

    Parameters
    ----------
    points:
        ``(n_points, 3)``.
    vertices:
        ``(n_triangles, 3, 3)``.

    Returns
    -------
    ndarray, shape ``(n_points, n_triangles)``, in steradians (always ≥ 0).
    """
    r1 = vertices[None, :, 0, :] - points[:, None, :]
    r2 = vertices[None, :, 1, :] - points[:, None, :]
    r3 = vertices[None, :, 2, :] - points[:, None, :]
    n1, n2, n3 = (np.linalg.norm(r, axis=-1) for r in (r1, r2, r3))
    triple = np.einsum("...i,...i->...", r1, np.cross(r2, r3))
    denom = (n1 * n2 * n3
             + np.einsum("...i,...i->...", r1, r2) * n3
             + np.einsum("...i,...i->...", r1, r3) * n2
             + np.einsum("...i,...i->...", r2, r3) * n1)
    return np.abs(2.0 * np.arctan2(triple, denom))


def detector_signal_arrays(
    points: np.ndarray,
    normals: np.ndarray,
    source_position: np.ndarray,
    vertices: np.ndarray,
    *,
    centroids: Optional[np.ndarray] = None,
    material: Material = TI64,
    model: ScatteringModel = lambertian,
    max_memory_mb: float = 256.0,
) -> np.ndarray:
    """
    Array-level forward model (no grid assumptions).

    For every surface point ``p`` with normal ``n`` the signal is::

        S(p) = Σ_triangles  g(cos_in, cos_out, cos_spec) · Ω_triangle(p)

    where ``g`` is the scattering model, the angles are evaluated towards the
    triangle centroid and ``Ω`` is the exact solid angle of the triangle.
    Shadowing and multiple scattering are not modelled.

    Parameters
    ----------
    points, normals:
        Arrays of shape ``(..., 3)``. Normals need not be unit length.
    source_position:
        ``(3,)`` beam pivot point.
    vertices:
        ``(n_triangles, 3, 3)`` detector triangles.
    centroids:
        Optional ``(n_triangles, 3)`` points used for the emission angles
        (default: triangle centroids).
    max_memory_mb:
        Approximate working-memory budget; triangles are processed in chunks to stay below it.

    Returns
    -------
    ndarray of shape ``points.shape[:-1]``.
    """
    out_shape = points.shape[:-1]
    P = np.asarray(points, dtype=float).reshape(-1, 3)
    N = _unit(np.asarray(normals, dtype=float).reshape(-1, 3))
    V = np.asarray(vertices, dtype=float)
    C = V.mean(axis=1) if centroids is None else np.asarray(centroids, dtype=float)

    to_source = _unit(np.asarray(source_position, dtype=float) - P)           # (P, 3)
    cos_in = np.einsum("pi,pi->p", to_source, N)[:, None]                     # (P, 1)
    specular = 2.0 * cos_in * N - to_source                                   # mirror of the beam about n

    # ~15 float64 temporaries of shape (P, chunk, 3) are alive at peak
    chunk = max(1, int(max_memory_mb * 2**20 / (len(P) * 3 * 8 * 15)))
    signal = np.zeros(len(P))
    for start in range(0, len(V), chunk):
        v, c = V[start:start + chunk], C[start:start + chunk]
        to_det = _unit(c[None, :, :] - P[:, None, :])                         # (P, e, 3)
        cos_out = np.einsum("pei,pi->pe", to_det, N)
        cos_spec = np.einsum("pei,pi->pe", to_det, specular)
        g = np.maximum(model(cos_in, cos_out, cos_spec, material), 0.0)
        signal += np.sum(g * triangle_solid_angles(P, v), axis=1)
    return signal.reshape(out_shape)


def simulate_signal(
    surface: SurfaceGrid,
    detector: Detector,
    source: ElectronSource,
    *,
    material: Material = TI64,
    model: ScatteringModel = lambertian,
    normalize: bool = True,
    max_memory_mb: float = 256.0,
) -> np.ndarray:
    """
    Signal image ``(ny, nx)`` collected by one detector while the beam scans ``surface``.

    With ``normalize=True`` (default) the image is divided by the signal from a
    flat surface on the same grid. This removes the purely geometric shading
    (detector solid angle varying with beam position) — the simulated equivalent
    of flat-field correcting with a reference plate. With ``normalize=False``
    the result is in backscattered electrons per primary electron.
    """
    kw = dict(material=material, model=model, max_memory_mb=max_memory_mb)
    pos = np.asarray(source.position)
    raw = detector_signal_arrays(surface.points, surface.normals(), pos, detector.vertices, **kw)
    if not normalize:
        return raw
    flat = surface.flat()
    ref = detector_signal_arrays(flat.points, flat.normals(), pos, detector.vertices, **kw)
    return raw / ref


def simulate_images(
    surface: SurfaceGrid,
    detectors: Union[Detector, Sequence[Detector]],
    source: ElectronSource,
    **kwargs,
) -> np.ndarray:
    """
    Simulate one image per detector; returns an array of shape ``(n_detectors, ny, nx)``.

    Keyword arguments are passed to :func:`simulate_signal`.
    """
    if isinstance(detectors, Detector):
        detectors = [detectors]
    return np.stack([simulate_signal(surface, d, source, **kwargs) for d in detectors])
