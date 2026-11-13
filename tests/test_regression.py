"""
Regression against the original research code.

The reference files were produced by ``tests/data/make_reference.py`` and
``tests/data/make_paper_reference.py`` (the latter reproduces the scenarios of
the publication notebook). They are compared at two levels:

1. **Exact** — the package's forward model, fed the original code's inputs
   (its surface-normal convention, grid spacing and float16 detector vertices),
   must reproduce the original images; the Sylvester port must reproduce the
   original reconstruction.
2. **End-to-end** — the package with its own (corrected) conventions must reach
   the same reconstruction accuracy as the original pipeline.
"""
from pathlib import Path

import numpy as np
import pytest

import bse_recon as br

DATA = Path(__file__).parent / "data"
PAPER = np.load(DATA / "reference_paper.npz")
SCENARIOS = sorted({k.split("/")[0] for k in PAPER.files})
SOURCE = np.array([0.0, 0.0, 1600.0])


def _scenario(name):
    g = {k.split("/")[1]: PAPER[k] for k in PAPER.files if k.startswith(name + "/")}
    n, x0, y0, offset, height = g["setup"]
    n = int(n)
    x, y = np.linspace(x0 - 20, x0 + 20, n), np.linspace(y0 - 20, y0 + 20, n)
    X, Y = np.meshgrid(x, y)
    if name.startswith("diam"):
        H = 0.5 * np.cos(2 * X)
    else:
        amp = 0.5 if name.startswith("asym") else 0.25
        H = amp * (np.cos(2 * X) + np.cos(2 * Y))
    return g, x, y, X, Y, H, float(offset), float(height)


def _legacy_inputs(X, Y, H, spacing):
    """Normals as the original code built them: (+gx, +gy, 1) with spacing = extent / n."""
    gy, gx = np.gradient(H, spacing, spacing)
    normals = np.stack([gx, gy, np.ones_like(gx)], axis=-1)
    flat_n = np.zeros_like(normals)
    flat_n[..., 2] = 1.0
    return np.stack([X, Y, H], -1), normals, np.stack([X, Y, 0 * H], -1), flat_n


@pytest.mark.parametrize("name", SCENARIOS)
def test_forward_model_reproduces_original_images(name):
    g, x, y, X, Y, H, *_ = _scenario(name)
    P, N, Pf, Nf = _legacy_inputs(X, Y, H, 40.0 / len(x))
    for i in range(4):
        v16 = g["vertices"][i].astype(np.float16).astype(float)  # original stored vertices as float16
        cen = g["vertices"][i].mean(axis=1)
        img = (br.detector_signal_arrays(P, N, SOURCE, v16, centroids=cen)
               / br.detector_signal_arrays(Pf, Nf, SOURCE, v16, centroids=cen))
        np.testing.assert_allclose(img, g["intensities"][i], rtol=1e-6)


@pytest.mark.parametrize("name", SCENARIOS)
def test_sylvester_port_reproduces_original(name):
    g, x, *_ = _scenario(name)
    d = 40.0 / len(x)
    z = br.integrate_gradients(g["nsx"], g["nsy"], d, d)
    ref = g["zcalc"] - g["zcalc"].mean()
    np.testing.assert_allclose(z, ref, atol=1e-6)


@pytest.mark.parametrize("name", SCENARIOS)
def test_end_to_end_accuracy_matches_original(name):
    g, x, y, X, Y, H, offset, height = _scenario(name)
    surface = br.SurfaceGrid(x, y, H)
    detectors = [br.Detector(v) for v in g["vertices"]]
    images = br.simulate_images(surface, detectors, br.ElectronSource(SOURCE))
    gx, gy = br.estimate_slopes(images, surface.X, surface.Y, offset, height)
    z = br.integrate_gradients(gx, gy, surface.dx, surface.dy)
    mae_pkg = br.height_errors(H, z)["mae"]
    mae_orig = br.height_errors(H, g["zcalc"])["mae"]
    assert mae_pkg == pytest.approx(mae_orig, rel=0.05)


def test_package_detectors_match_research_geometry():
    """four_quadrant_detectors reproduces the notebook geometry up to its 0.04 mm y-centring quirk."""
    g, *_ = _scenario("origins_0_0")
    dets = br.four_quadrant_detectors(diameter=80, tilt=75, offset=132, height=87)
    for det, ref in zip(dets, g["vertices"]):
        assert np.abs(det.vertices - ref).max() < 0.05


def test_simple_sine_reference():
    r = np.load(DATA / "reference_sine.npz")
    n, *_ = r["params"]
    n = int(n)
    x = np.linspace(-10, 10, n)
    X, Y = np.meshgrid(x, x)
    P, N, Pf, Nf = _legacy_inputs(X, Y, r["heights"], 20.0 / n)
    for i, V in enumerate(r["detector_vertices"]):
        v16 = V.astype(np.float16).astype(float)
        img = (br.detector_signal_arrays(P, N, SOURCE, v16, centroids=V.mean(1))
               / br.detector_signal_arrays(Pf, Nf, SOURCE, v16, centroids=V.mean(1)))
        np.testing.assert_allclose(img, r["intensities"][i], rtol=1e-7)
