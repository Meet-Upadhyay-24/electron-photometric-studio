"""Checks against analytic results and symmetry arguments."""
import numpy as np
import pytest

import bse_recon as br

SOURCE = br.ElectronSource([0, 0, 1600])
# The shipped mesh is very slightly asymmetric (node-mean y = 0.003 mesh units), so
# exact-symmetry tests use a generated, perfectly symmetric disk instead.
SYM_DISK = br.Detector.disk(1.0, n_rings=3, n_sectors=16)


def sym_detectors(**kw):
    return br.four_quadrant_detectors(mesh=SYM_DISK, **kw)


# -- solid angles ---------------------------------------------------------------
def test_closed_surface_subtends_4pi():
    # regular tetrahedron around an interior point
    v = np.array([[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]], dtype=float)
    faces = np.array([[v[a], v[b], v[c]] for a, b, c in [(0, 1, 2), (0, 3, 1), (0, 2, 3), (1, 3, 2)]])
    omega = br.triangle_solid_angles(np.array([[0.1, -0.2, 0.05]]), faces).sum()
    assert omega == pytest.approx(4 * np.pi, rel=1e-12)


def test_disk_on_axis_solid_angle():
    r, h = 20.0, 87.0
    det = br.Detector.disk(2 * r, n_rings=40, n_sectors=256).translated(0, 0, h)
    omega = br.triangle_solid_angles(np.zeros((1, 3)), det.vertices).sum()
    exact = 2 * np.pi * (1 - h / np.hypot(h, r))
    assert omega == pytest.approx(exact, rel=2e-4)  # polygonal disk is slightly smaller


# -- forward model -----------------------------------------------------------
def test_tilted_plane_faces_the_left_detector():
    s = br.SurfaceGrid.from_extent(5, 5, (-1, 1), (-1, 1))
    s = s.with_heights(0.1 * s.X)                                   # rises towards +x
    left, top, right, bottom = br.simulate_images(s, sym_detectors(), SOURCE)
    assert np.all(left > right)
    np.testing.assert_allclose(top, bottom[::-1, :], rtol=1e-6)  # mirror pair in y


def test_flat_surface_contrast_is_antisymmetric():
    s = br.SurfaceGrid.from_extent(21, 21, (-10, 10), (-10, 10))
    images = br.simulate_images(s, sym_detectors(), SOURCE, normalize=False)
    Mx, My = br.contrast(images)
    # vertices are mirror-symmetric but the quad-split diagonals are not → ~1e-7 residual
    np.testing.assert_allclose(Mx, -Mx[:, ::-1], atol=1e-6)
    np.testing.assert_allclose(My, -My[::-1, :], atol=1e-6)
    assert abs(Mx[10, 10]) < 1e-6 and abs(My[10, 10]) < 1e-6


def test_normalized_flat_surface_is_unity():
    s = br.SurfaceGrid.from_extent(9, 9, (-10, 10), (-10, 10))
    images = br.simulate_images(s, br.four_quadrant_detectors(), SOURCE)
    np.testing.assert_allclose(images, 1.0, rtol=1e-12)


def test_chunking_does_not_change_result():
    s = br.SurfaceGrid.from_extent(15, 12, (-5, 5), (-4, 4))
    s = s.with_heights(br.gaussian_bump(s, 0.5, 2.0))
    det = br.four_quadrant_detectors()[2]
    a = br.simulate_signal(s, det, SOURCE, max_memory_mb=1e-6)    # one triangle per chunk
    b = br.simulate_signal(s, det, SOURCE)
    np.testing.assert_allclose(a, b, rtol=1e-12)


def test_custom_scattering_model_is_used():
    s = br.SurfaceGrid.from_extent(5, 5, (-1, 1), (-1, 1))
    det = br.four_quadrant_detectors()[0]
    double = lambda ci, co, cs, m: 2 * br.lambertian(ci, co, cs, m)  # noqa: E731
    a = br.simulate_signal(s, det, SOURCE, normalize=False)
    b = br.simulate_signal(s, det, SOURCE, normalize=False, model=double)
    np.testing.assert_allclose(b, 2 * a)


# -- inversion ---------------------------------------------------------------
@pytest.mark.parametrize("azimuth", [0.0, 30.0, 45.0])
def test_plane_slopes_are_recovered(azimuth):
    # small window: the closed-form inversion neglects the sample height relative to the detector height
    s = br.SurfaceGrid.from_extent(11, 11, (-3, 3), (-3, 3))
    s = s.with_heights(0.05 * s.X - 0.03 * s.Y)
    dets = sym_detectors(diameter=1.0, tilt=0, offset=132, height=87, azimuth=azimuth)
    gx, gy = br.estimate_slopes(br.simulate_images(s, dets, SOURCE), s.X, s.Y, 132, 87, azimuth=azimuth)
    # relief up to 0.24 mm vs 87 mm detector height → expected relative error ~ Δz/height ≈ 0.3 %
    np.testing.assert_allclose(gx, 0.05, rtol=0.01)
    np.testing.assert_allclose(gy, -0.03, rtol=0.01)


def test_point_detector_inversion_is_exact():
    """The closed-form inversion is exact for the point-detector Lambertian model it is derived from."""
    rng = np.random.default_rng(0)
    gx, gy = rng.uniform(-0.3, 0.3, 50), rng.uniform(-0.3, 0.3, 50)
    x, y = rng.uniform(-20, 20, 50), rng.uniform(-20, 20, 50)
    w, h = 132.0, 87.0

    def lam(dx, dy):  # n · (d − p), n = (−gx, −gy, 1)
        return -gx * (dx - x) - gy * (dy - y) + h

    images = [lam(-w, 0), lam(0, w), lam(w, 0), lam(0, -w)]
    ex, ey = br.estimate_slopes(images, x, y, w, h)
    np.testing.assert_allclose(ex, gx, atol=1e-12)
    np.testing.assert_allclose(ey, gy, atol=1e-12)


@pytest.mark.parametrize("boundary", ["free", "dirichlet"])
def test_integration_of_exact_gradients(boundary):
    s = br.SurfaceGrid.from_extent(60, 50, (-3, 3), (-2, 2))
    if boundary == "dirichlet":
        h = br.gaussian_bump(s, 1.0, 0.4)           # ~0 at the border, as Dirichlet assumes
    else:
        h = 0.3 * s.X**2 - 0.2 * s.X * s.Y + 0.1 * s.Y
    s = s.with_heights(h)
    z = br.integrate_gradients(*s.gradients(), s.dx, s.dy, boundary=boundary)
    assert br.height_errors(h, z)["max"] < 0.02 * np.ptp(h)


def test_tikhonov_shrinks_heights():
    s = br.SurfaceGrid.from_extent(30, 30, (-3, 3), (-3, 3))
    s = s.with_heights(br.sine_surface(s, 0.2, 1.5))
    z0 = br.integrate_gradients(*s.gradients(), s.dx, s.dy)
    z1 = br.integrate_gradients(*s.gradients(), s.dx, s.dy, regularization=1.0)
    assert np.std(z1) < np.std(z0)


def test_remove_background_plane():
    s = br.SurfaceGrid.from_extent(20, 20, (0, 1), (0, 1))
    np.testing.assert_allclose(br.remove_background(s.X, s.Y, 3 + 2 * s.X - s.Y), 0, atol=1e-12)
