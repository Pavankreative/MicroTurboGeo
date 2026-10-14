"""Tests for the parametric section generator."""
import numpy as np
import pytest

from microturbogeo import sections, geom


def test_ps_ss_share_endpoints():
    ps, ss = sections.make_blade_sections([0.03, 0.035], 0.01, [0, 0], [60, 55], [0.15, 0.12], 0.0002)
    for i in range(2):
        for key in ('x', 'y', 'z'):
            assert ps[i][0][key] == pytest.approx(ss[i][0][key], abs=1e-15)
            assert ps[i][-1][key] == pytest.approx(ss[i][-1][key], abs=1e-12)


def test_sections_lie_on_cylinders():
    radii = [0.03, 0.036]
    ps, ss = sections.make_blade_sections(radii, 0.01, [0, 0], [60, 60], [0.15, 0.15], 0.0002)
    for i, r in enumerate(radii):
        np.testing.assert_allclose(geom.radial_of(ps[i]), r, rtol=1e-12)
        np.testing.assert_allclose(geom.radial_of(ss[i]), r, rtol=1e-12)


def test_pressure_side_is_on_positive_theta():
    ps, ss = sections.make_blade_sections([0.03], 0.01, [0], [60], [0.15], 0.0002)
    th_ps = np.arctan2(ps[0]['y'], ps[0]['x'])
    th_ss = np.arctan2(ss[0]['y'], ss[0]['x'])
    mid = len(th_ps) // 2
    assert th_ps[mid] > th_ss[mid]


def test_max_thickness_matches_request():
    chord, t = 0.02, 0.12
    ps2, ss2 = sections.build_section_2d(chord, 0, 0, t, 0.0, n_pts=400)   # zero camber, sharp TE
    thickness = np.max(ps2[:, 1] - ss2[:, 1])
    assert thickness == pytest.approx(t * chord, rel=0.02)


def test_blunt_trailing_edge_thickness():
    ps2, ss2 = sections.build_section_2d(0.01, 0, 0, 0.1, 0.0003, n_pts=100)
    # the SS curve ends on the PS corner; the corners are t_te apart in r*theta
    te_corner_ss = ss2[-1]
    assert np.allclose(te_corner_ss, ps2[-1])
    assert abs(ps2[-1][1]) == pytest.approx(0.0003 / 2, rel=1e-3)


def test_fpd_roundtrip(tmp_path):
    ps, _ = sections.make_blade_sections([0.03, 0.035], 0.01, [0, 0], [60, 55], [0.15, 0.12], 0.0002)
    sections.write_fpd(tmp_path / 'ps.fpd', ps)
    back = geom.load_curves_from_fpd(tmp_path / 'ps.fpd')
    np.testing.assert_allclose(back['x'], ps['x'], rtol=1e-13)
    assert back.shape == ps.shape


def test_radii_must_ascend():
    with pytest.raises(ValueError):
        sections.make_blade_sections([0.035, 0.03], 0.01, [0, 0], [1, 1], [0.1, 0.1], 0.0)


def test_free_vortex_angle_keeps_rtan_constant():
    r = np.array([0.03, 0.035, 0.04])
    a = sections.free_vortex_angle(r, 0.035, 65.0)
    k = r * np.tan(np.radians(a))
    np.testing.assert_allclose(k, k[0], rtol=1e-12)


def test_mirror_theta():
    ps, _ = sections.make_blade_sections([0.03], 0.01, [0], [60], [0.15], 0.0002)
    m = sections.mirror_theta(ps)
    np.testing.assert_allclose(m['y'], -ps['y'])
    np.testing.assert_allclose(m['z'], ps['z'])
