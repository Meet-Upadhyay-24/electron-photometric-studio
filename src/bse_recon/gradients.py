"""From four detector images to surface slopes (photometric stereo inversion)."""
from __future__ import annotations

from typing import Sequence, Tuple

import numpy as np

_EPS = 1e-12


def contrast(images: Sequence[np.ndarray]) -> Tuple[np.ndarray, np.ndarray]:
    """
    Normalised difference signals from four opposed detectors.

    Parameters
    ----------
    images:
        ``[left, top, right, bottom]`` — the order returned by
        :func:`bse_recon.four_quadrant_detectors`.

    Returns
    -------
    Mx, My:
        ``Mx = (L − R)/(L + R)`` and ``My = (B − T)/(B + T)``. A surface rising
        towards +x faces the left detector, so ``Mx > 0`` where ``∂h/∂x > 0``.
    """
    left, top, right, bottom = (np.asarray(im, dtype=float) for im in images)
    Mx = (left - right) / (left + right + _EPS)
    My = (bottom - top) / (bottom + top + _EPS)
    return Mx, My


def slopes_from_contrast(
    Mx: np.ndarray,
    My: np.ndarray,
    x: np.ndarray,
    y: np.ndarray,
    offset: float,
    height: float,
    azimuth: float = 0.0,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Convert contrast signals into surface slopes ``(∂h/∂x, ∂h/∂y)``.

    Treating each detector as a point at ``(±offset, 0, height)`` /
    ``(0, ±offset, height)`` and the emission as Lambertian, the contrasts are::

        Mx = gx·w / (gx·x + gy·y + h),   My = gy·w / (gx·x + gy·y + h)

    with ``w = offset`` and ``h = height``, which inverts exactly to::

        gx = Mx·h / (w − Mx·x − My·y),   gy = My·h / (w − Mx·x − My·y)

    The ``x, y`` terms correct for the beam position moving off-axis. The
    derivation assumes small detectors and surface relief much smaller than
    ``height``; a relief of ``Δz`` causes a relative slope error of roughly ``Δz / height``.

    Parameters
    ----------
    x, y:
        Beam position of every pixel (2-D arrays, e.g. ``surface.X, surface.Y``).
    offset, height:
        Detector distance from the beam axis and height above the sample.
    azimuth:
        Rotation (degrees, about z) applied to the detector set, as in
        :func:`bse_recon.four_quadrant_detectors`. The slopes are computed in the
        detector frame and rotated back to the x/y frame.
    """
    a = np.deg2rad(azimuth)
    c, s = np.cos(a), np.sin(a)
    xd, yd = c * x + s * y, -s * x + c * y          # pixel positions in the detector frame
    denom = offset - Mx * xd - My * yd
    gxd, gyd = Mx * height / denom, My * height / denom
    return c * gxd - s * gyd, s * gxd + c * gyd


def estimate_slopes(images, x, y, offset: float, height: float, azimuth: float = 0.0):
    """:func:`contrast` followed by :func:`slopes_from_contrast`."""
    Mx, My = contrast(images)
    return slopes_from_contrast(Mx, My, x, y, offset, height, azimuth)
