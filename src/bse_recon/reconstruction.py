"""Integrate slope fields into height maps, and compare heights."""
from __future__ import annotations

from typing import Dict

import numpy as np
from scipy.linalg import solve_sylvester


def _difference_matrix(n: int, h: float) -> np.ndarray:
    """Second-order accurate first-derivative matrix (central inside, one-sided at the ends)."""
    if n < 3:
        raise ValueError("Need at least 3 samples per axis.")
    D = np.zeros((n, n))
    i = np.arange(1, n - 1)
    D[i, i - 1], D[i, i + 1] = -1.0, 1.0
    D[0, :3] = [-3.0, 4.0, -1.0]
    D[-1, -3:] = [1.0, -4.0, 3.0]
    return D / (2.0 * h)


def integrate_gradients(
    gx: np.ndarray,
    gy: np.ndarray,
    dx: float,
    dy: float,
    *,
    regularization: float = 0.0,
    boundary: str = "free",
) -> np.ndarray:
    """
    Least-squares surface from slope fields (Harker & O'Leary style Sylvester solver).

    Solves ``min ‖Dy Z − gy‖² + ‖Z Dxᵀ − gx‖² (+ λ²‖Z‖²)`` exactly via the
    Sylvester equation ``A Z + Z B = C``.

    Parameters
    ----------
    gx, gy:
        Slopes ``∂h/∂x`` and ``∂h/∂y``, shape ``(ny, nx)``.
    dx, dy:
        Grid spacing.
    regularization:
        Tikhonov weight λ pulling heights towards zero (useful for noisy data).
    boundary:
        ``"free"`` (default) or ``"dirichlet"`` (height fixed to zero just outside the grid).

    Returns
    -------
    ndarray ``(ny, nx)`` with zero mean (the absolute offset is not observable from slopes).
    """
    gx, gy = np.asarray(gx, dtype=float), np.asarray(gy, dtype=float)
    ny, nx = gx.shape
    lam2 = float(regularization) ** 2

    if boundary == "dirichlet":
        Dx, Dy = _difference_matrix(nx + 2, dx), _difference_matrix(ny + 2, dy)
        Pn = np.zeros((nx + 2, nx)); Pn[1:-1] = np.eye(nx)  # noqa: E702
        Pm = np.zeros((ny + 2, ny)); Pm[1:-1] = np.eye(ny)  # noqa: E702
        gxp, gyp = np.pad(gx, 1), np.pad(gy, 1)
        DyP, DxP = Dy @ Pm, Dx @ Pn
        A = DyP.T @ DyP + lam2 * np.eye(ny)
        B = DxP.T @ DxP + lam2 * np.eye(nx)
        C = DyP.T @ gyp @ Pn + Pm.T @ gxp @ DxP
    elif boundary == "free":
        Dx, Dy = _difference_matrix(nx, dx), _difference_matrix(ny, dy)
        # The free problem is singular (constant offset); a vanishing ridge keeps the solver well posed.
        lam2 = max(lam2, 1e-12 * (np.trace(Dx.T @ Dx) / nx + np.trace(Dy.T @ Dy) / ny))
        A = Dy.T @ Dy + lam2 * np.eye(ny)
        B = Dx.T @ Dx + lam2 * np.eye(nx)
        C = Dy.T @ gy + gx @ Dx
    else:
        raise ValueError("boundary must be 'free' or 'dirichlet'.")

    Z = solve_sylvester(A, B, C)
    return Z - Z.mean()


def remove_background(x: np.ndarray, y: np.ndarray, z: np.ndarray, degree: int = 1) -> np.ndarray:
    """Subtract the least-squares 2-D polynomial of total degree ``degree`` (1 = plane)."""
    xs, ys, zs = (np.ravel(a) for a in (x, y, z))
    cols = [xs**i * ys**j for i in range(degree + 1) for j in range(degree + 1 - i)]
    A = np.column_stack(cols)
    coef, *_ = np.linalg.lstsq(A, zs, rcond=None)
    return z - (A @ coef).reshape(np.shape(z))


def height_errors(z_true: np.ndarray, z_est: np.ndarray) -> Dict[str, float]:
    """RMS, mean-absolute and maximum-absolute error after removing the mean offset of each map."""
    d = (np.asarray(z_est) - np.mean(z_est)) - (np.asarray(z_true) - np.mean(z_true))
    return {"rmse": float(np.sqrt(np.mean(d**2))),
            "mae": float(np.mean(np.abs(d))),
            "max": float(np.max(np.abs(d)))}
