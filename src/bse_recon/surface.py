"""Sample surfaces: a regular height-map grid plus a few synthetic generators."""
from __future__ import annotations

from typing import Optional, Tuple

import numpy as np


class SurfaceGrid:
    """
    A height map ``z = h(x, y)`` sampled on a regular grid.

    All lengths (``x``, ``y``, heights, detector and source positions) must use
    the same unit — the examples use millimetres.

    Parameters
    ----------
    x, y:
        1-D, strictly increasing, evenly spaced coordinates of the columns and rows.
    heights:
        Array of shape ``(len(y), len(x))``; defaults to a flat surface.
    """

    def __init__(self, x, y, heights: Optional[np.ndarray] = None):
        self.x = np.asarray(x, dtype=float)
        self.y = np.asarray(y, dtype=float)
        if self.x.ndim != 1 or self.y.ndim != 1 or len(self.x) < 2 or len(self.y) < 2:
            raise ValueError("x and y must be 1-D with at least two samples each.")
        self.X, self.Y = np.meshgrid(self.x, self.y)
        if heights is None:
            heights = np.zeros(self.shape)
        self.heights = np.asarray(heights, dtype=float)
        if self.heights.shape != self.shape:
            raise ValueError(f"heights has shape {self.heights.shape}, expected {self.shape}.")

    @classmethod
    def from_extent(cls, nx: int, ny: int, x_range: Tuple[float, float], y_range: Tuple[float, float]) -> "SurfaceGrid":
        """Flat grid with ``nx × ny`` nodes spanning ``x_range`` and ``y_range`` (inclusive)."""
        return cls(np.linspace(*x_range, nx), np.linspace(*y_range, ny))

    # -- geometry -----------------------------------------------------------
    @property
    def shape(self) -> Tuple[int, int]:
        return (len(self.y), len(self.x))

    @property
    def dx(self) -> float:
        return float(self.x[1] - self.x[0])

    @property
    def dy(self) -> float:
        return float(self.y[1] - self.y[0])

    @property
    def points(self) -> np.ndarray:
        """Node positions, shape ``(ny, nx, 3)``."""
        return np.stack([self.X, self.Y, self.heights], axis=-1)

    def gradients(self) -> Tuple[np.ndarray, np.ndarray]:
        """Surface slopes ``(∂h/∂x, ∂h/∂y)`` by central differences (one-sided at edges)."""
        gy, gx = np.gradient(self.heights, self.dy, self.dx)
        return gx, gy

    def normals(self) -> np.ndarray:
        """Unit outward (upward) surface normals ``(-∂h/∂x, -∂h/∂y, 1)/‖·‖``, shape ``(ny, nx, 3)``."""
        gx, gy = self.gradients()
        n = np.stack([-gx, -gy, np.ones_like(gx)], axis=-1)
        return n / np.linalg.norm(n, axis=-1, keepdims=True)

    # -- convenience --------------------------------------------------------
    def with_heights(self, heights: np.ndarray) -> "SurfaceGrid":
        """A new grid on the same coordinates with different heights."""
        return SurfaceGrid(self.x, self.y, heights)

    def flat(self) -> "SurfaceGrid":
        """A flat (``h = 0``) grid on the same coordinates."""
        return SurfaceGrid(self.x, self.y)

    def __repr__(self) -> str:
        return (f"SurfaceGrid({self.shape[1]}x{self.shape[0]}, x=[{self.x[0]:g}, {self.x[-1]:g}], "
                f"y=[{self.y[0]:g}, {self.y[-1]:g}], h=[{self.heights.min():.3g}, {self.heights.max():.3g}])")


# -- synthetic height maps ------------------------------------------------------
def sine_surface(grid: SurfaceGrid, amplitude: float, wavenumber: float = 1.0) -> np.ndarray:
    """``amplitude · (cos(k x) + cos(k y))`` — the egg-crate test surface (peak-to-valley 4·amplitude)."""
    return amplitude * (np.cos(wavenumber * grid.X) + np.cos(wavenumber * grid.Y))


def gaussian_bump(grid: SurfaceGrid, amplitude: float, sigma: float,
                  center: Tuple[float, float] = (0.0, 0.0)) -> np.ndarray:
    """A single Gaussian hill (``amplitude > 0``) or pit (``amplitude < 0``)."""
    r2 = (grid.X - center[0]) ** 2 + (grid.Y - center[1]) ** 2
    return amplitude * np.exp(-r2 / (2.0 * sigma**2))


def random_surface(grid: SurfaceGrid, rms: float, correlation_length: float, seed: Optional[int] = None) -> np.ndarray:
    """Gaussian random rough surface with the given RMS height and correlation length."""
    from scipy.ndimage import gaussian_filter

    rng = np.random.default_rng(seed)
    h = gaussian_filter(rng.standard_normal(grid.shape),
                        sigma=(correlation_length / grid.dy, correlation_length / grid.dx), mode="wrap")
    h -= h.mean()
    return h * (rms / h.std())


def edge_taper(shape: Tuple[int, int], fraction: float = 0.15) -> np.ndarray:
    """
    Separable cosine window that rolls off to zero over ``fraction`` of each edge.

    Multiply a height map by it to avoid slope discontinuities at the borders.
    """
    def window(n):
        w = np.ones(n)
        m = max(int(round(fraction * n)), 1)
        ramp = 0.5 * (1 - np.cos(np.pi * np.arange(m) / m))
        w[:m] = ramp
        w[-m:] = ramp[::-1]
        return w

    return np.outer(window(shape[0]), window(shape[1]))
