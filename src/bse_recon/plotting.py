"""Matplotlib helpers. All functions return the figure so you can save or tweak it."""
from __future__ import annotations

from typing import Optional, Sequence

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

from .detector import Detector
from .source import ElectronSource
from .surface import SurfaceGrid

_DET_COLORS = ["tab:blue", "tab:orange", "tab:green", "tab:red", "tab:purple", "tab:brown"]


def apply_science_style() -> None:
    """Use the SciencePlots ``science`` style if installed (``pip install bse-recon[science]``); otherwise no-op."""
    try:
        import scienceplots  # noqa: F401
        plt.style.use(["science", "no-latex"])
    except Exception:
        pass


def plot_setup(
    detectors: Sequence[Detector],
    surface: Optional[SurfaceGrid] = None,
    source: Optional[ElectronSource] = None,
    *,
    show_beam: bool = True,
    ax=None,
):
    """3-D view of the detectors, the scanned area and the beam."""
    if ax is None:
        fig = plt.figure(figsize=(7, 6))
        ax = fig.add_subplot(projection="3d")
    fig = ax.figure
    for i, det in enumerate(detectors):
        col = _DET_COLORS[i % len(_DET_COLORS)]
        ax.add_collection3d(Poly3DCollection(det.vertices, facecolor=col, edgecolor="k", linewidth=0.2,
                                             alpha=0.8, label=det.name or f"detector {i}"))
    pts = [d.vertices.reshape(-1, 3) for d in detectors]
    if surface is not None:
        x0, x1, y0, y1 = surface.x[0], surface.x[-1], surface.y[0], surface.y[-1]
        ax.plot([x0, x1, x1, x0, x0], [y0, y0, y1, y1, y0], 0, color="k")
        pts.append(np.array([[x0, y0, 0], [x1, y1, 0]]))
    if source is not None and show_beam:
        p = source.position
        ax.scatter(*p, color="k", marker="v")
        ax.plot([p[0], 0], [p[1], 0], [p[2], 0], "k--", linewidth=0.8)
    pts = np.vstack(pts)
    lo, hi = pts.min(axis=0), pts.max(axis=0)
    ax.set(xlim=(lo[0], hi[0]), ylim=(lo[1], hi[1]), zlim=(min(lo[2], 0), hi[2]),
           xlabel="x", ylabel="y", zlabel="z")
    ax.set_box_aspect(np.maximum(hi - lo, 1e-9))
    ax.legend(loc="upper left")
    return fig


def plot_detector_2d(detector: Detector, *, ax=None, show_centroids: bool = False):
    """Top-down (x-y) view of a detector mesh."""
    if ax is None:
        _, ax = plt.subplots(figsize=(5, 5))
    for tri in detector.vertices:
        closed = np.vstack([tri, tri[:1]])
        ax.plot(closed[:, 0], closed[:, 1], "k-", linewidth=0.5)
    if show_centroids:
        c = detector.centroids
        ax.plot(c[:, 0], c[:, 1], "r.", markersize=3)
    ax.set(aspect="equal", xlabel="x", ylabel="y", title=detector.name)
    return ax.figure


def _imshow(ax, surface, data, title, cmap="viridis", vmin=None, vmax=None, label=None):
    ext = [surface.x[0], surface.x[-1], surface.y[0], surface.y[-1]]
    im = ax.imshow(data, origin="lower", extent=ext, cmap=cmap, vmin=vmin, vmax=vmax)
    ax.set(title=title, xlabel="x", ylabel="y")
    ax.figure.colorbar(im, ax=ax, shrink=0.8, label=label)
    return im


def plot_images(images: Sequence[np.ndarray], surface: SurfaceGrid,
                names: Sequence[str] = ("left", "top", "right", "bottom")):
    """One panel per detector image."""
    n = len(images)
    fig, axes = plt.subplots(1, n, figsize=(3.6 * n, 3.4), squeeze=False)
    for ax, im, name in zip(axes[0], images, list(names) + [""] * n):
        _imshow(ax, surface, im, name)
    fig.tight_layout()
    return fig


def plot_slopes(surface: SurfaceGrid, gx_est: np.ndarray, gy_est: np.ndarray):
    """True vs estimated slopes, their difference, and a centre-line profile."""
    gx, gy = surface.gradients()
    fig, axes = plt.subplots(2, 4, figsize=(15, 6.5))
    for row, (true, est, lab) in enumerate([(gx, gx_est, "∂h/∂x"), (gy, gy_est, "∂h/∂y")]):
        lim = np.abs(true).max() or None
        _imshow(axes[row, 0], surface, true, f"true {lab}", "RdBu_r", -lim if lim else None, lim)
        _imshow(axes[row, 1], surface, est, f"estimated {lab}", "RdBu_r", -lim if lim else None, lim)
        _imshow(axes[row, 2], surface, est - true, "estimated − true", "RdBu_r")
        ax = axes[row, 3]
        if row == 0:
            k = surface.shape[0] // 2
            ax.plot(surface.x, true[k], "k-", label="true")
            ax.plot(surface.x, est[k], "r--", label="estimated")
            ax.set(xlabel="x", title="centre row")
        else:
            k = surface.shape[1] // 2
            ax.plot(surface.y, true[:, k], "k-", label="true")
            ax.plot(surface.y, est[:, k], "r--", label="estimated")
            ax.set(xlabel="y", title="centre column")
        ax.legend()
    fig.tight_layout()
    return fig


def plot_reconstruction(surface: SurfaceGrid, z_est: np.ndarray):
    """True heights, reconstructed heights, their difference and a centre-line profile."""
    z_true = surface.heights - surface.heights.mean()
    z_est = z_est - z_est.mean()
    vmin, vmax = z_true.min(), z_true.max()
    fig, axes = plt.subplots(1, 4, figsize=(16, 3.6))
    _imshow(axes[0], surface, z_true, "true height", vmin=vmin, vmax=vmax)
    _imshow(axes[1], surface, z_est, "reconstructed", vmin=vmin, vmax=vmax)
    _imshow(axes[2], surface, z_est - z_true, "reconstructed − true", "RdBu_r")
    k = surface.shape[0] // 2
    axes[3].plot(surface.x, z_true[k], "k-", label="true")
    axes[3].plot(surface.x, z_est[k], "r--", label="reconstructed")
    axes[3].set(xlabel="x", ylabel="height", title="centre row")
    axes[3].legend()
    fig.tight_layout()
    return fig
