import numpy as np
import pytest

import bse_recon as br


def test_default_disk_loads():
    d = br.Detector.default_disk()
    assert len(d) == 90
    lo, hi = d.vertices.reshape(-1, 3).min(0), d.vertices.reshape(-1, 3).max(0)
    assert hi[0] - lo[0] == pytest.approx(5.64244128)


def test_generated_disk_area():
    d = br.Detector.disk(10.0, n_rings=10, n_sectors=200)
    assert d.area == pytest.approx(np.pi * 25, rel=1e-3)


def test_rectangle_area_and_scaling():
    d = br.Detector.rectangle(4, 2)
    assert d.area == pytest.approx(8.0)
    assert d.scaled_to_diameter(8).area == pytest.approx(32.0)


def test_transformations_are_immutable_and_rigid():
    d = br.Detector.default_disk()
    moved = d.rotated(30, "y").translated(1, 2, 3).rotated(45, "z")
    assert d.vertices is not moved.vertices
    assert moved.area == pytest.approx(d.area)
    with pytest.raises(ValueError):
        d.vertices[0, 0, 0] = 1.0


def test_four_quadrant_layout():
    dets = br.four_quadrant_detectors(diameter=40, tilt=25, offset=100, height=120)
    centers = np.array([d.center for d in dets])
    np.testing.assert_allclose(centers[:, 2], 120, atol=0.5)
    assert [d.name for d in dets] == ["left", "top", "right", "bottom"]
    assert centers[0, 0] < 0 < centers[2, 0] and centers[3, 1] < 0 < centers[1, 1]


def test_from_inp_roundtrip(tmp_path):
    p = tmp_path / "tri.inp"
    p.write_text("*Node\n1, 0., 0.\n2, 1., 0.\n3, 0., 1., 2.\n*Element, type=S3\n1, 1, 2, 3\n*End Part\n")
    d = br.Detector.from_inp(p)
    np.testing.assert_allclose(d.vertices[0], [[0, 0, 0], [1, 0, 0], [0, 1, 2]])


def test_surface_spacing_and_normals():
    s = br.SurfaceGrid.from_extent(41, 21, (-10, 10), (0, 5))
    assert s.dx == pytest.approx(0.5) and s.dy == pytest.approx(0.25)
    s = s.with_heights(0.2 * s.X)
    n = s.normals()
    np.testing.assert_allclose(np.linalg.norm(n, axis=-1), 1.0)
    assert np.all(n[..., 0] < 0)  # surface rising towards +x faces -x


def test_material():
    assert br.TI64.eta(1.0) == pytest.approx(br.TI64.eta0)
    assert br.TI64.eta(0.0) == pytest.approx(0.89)
    with pytest.raises(ValueError):
        br.Material("bad", ((26, 0.5),))


def test_surface_generators():
    s = br.SurfaceGrid.from_extent(64, 64, (0, 10), (0, 10))
    h = br.random_surface(s, rms=0.1, correlation_length=0.5, seed=1)
    assert h.std() == pytest.approx(0.1)
    t = br.edge_taper(s.shape, 0.2)
    assert t[0, 0] == 0 and t[32, 32] == 1
