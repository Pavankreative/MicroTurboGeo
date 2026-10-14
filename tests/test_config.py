"""Tests for blade-count/pitch control, tip-clearance options, parametric endwalls and axis validation."""
import numpy as np
import pytest

from microturbogeo import blade, geom, machine, stage


def test_n_blade_from_pitch():
    assert blade._resolve_blade_count({'pitch_deg': 15.0}) == 24
    assert blade._resolve_blade_count({'n_blade': 24, 'pitch_deg': 15.0}) == 24


def test_inconsistent_pitch_rejected():
    with pytest.raises(ValueError):
        blade._resolve_blade_count({'pitch_deg': 14.0})            # 360/14 is not an integer
    with pytest.raises(ValueError):
        blade._resolve_blade_count({'n_blade': 20, 'pitch_deg': 15.0})
    with pytest.raises(ValueError):
        blade._resolve_blade_count({})


def test_clearance_modes():
    b = blade.Blade()
    b.tip_clearance = 0.0003
    assert b.clearance_at(0.01) == (0.0003, 0.0003)
    b.tip_clearance_te = 0.0002
    assert b.clearance_at(0.01) == (0.0003, 0.0002)
    b.tip_clearance_mode = 'percent_span'
    b.tip_clearance, b.tip_clearance_te = 2.0, -1.0
    c_le, c_te = b.clearance_at(0.01)
    assert c_le == pytest.approx(0.0002) and c_te == pytest.approx(0.0002)


def test_parametric_endwall_curve():
    c = geom.make_parametric_meridional(0.03, 0.04, -0.01, 0.02, 11, 'cosine')
    r = geom.radial_of(c)
    assert r[0] == pytest.approx(0.03) and r[-1] == pytest.approx(0.04)
    assert c['z'][0] == -0.01 and c['z'][-1] == 0.02
    assert np.all(np.diff(r) >= 0)
    lin = geom.make_parametric_meridional(0.03, 0.04, -0.01, 0.02, 11, 'linear')
    assert geom.radial_of(lin)[5] == pytest.approx(0.035)


def test_parametric_endwall_validation():
    with pytest.raises(ValueError):
        geom.make_parametric_meridional(0.03, 0.04, 0.02, 0.0)
    with pytest.raises(ValueError):
        stage.Endwalls.from_config({'type': 'parametric', 'hub_radius': [0.03, 0.03]})
    with pytest.raises(ValueError):    # shroud inside hub
        stage.Endwalls.from_config({'type': 'parametric', 'hub_radius': [0.04, 0.04],
                                    'shroud_radius': [0.03, 0.03], 'z_inlet': -0.01, 'z_outlet': 0.02})


def test_radial_offset_curve():
    c = geom.make_parametric_meridional(0.03, 0.04, -0.01, 0.02, 5, 'linear')
    o = geom.radial_offset_curve(c, -0.0003)
    np.testing.assert_allclose(geom.radial_of(o), geom.radial_of(c) - 0.0003)
    with pytest.raises(ValueError):
        geom.radial_offset_curve(c, -1.0)


def test_interp_radius():
    c = geom.make_parametric_meridional(0.03, 0.04, 0.0, 0.01, 11, 'linear')
    assert geom.interp_radius(c, 0.005) == pytest.approx(0.035)


def test_axis_must_be_z():
    machine._check_axis(((0., 0., 0.), (0., 0., 1.)))
    machine._check_axis(((0., 0., 0.), (0., 0., 5.)))
    with pytest.raises(ValueError):
        machine._check_axis(((0., 0., 0.), (1., 0., 0.)))
    with pytest.raises(ValueError):
        machine._check_axis(((0., 0., 0.), (0., 0., -1.)))
    with pytest.raises(ValueError):
        machine._check_axis(((0.1, 0., 0.), (0.1, 0., 1.)))


def test_extrude_radially_handles_negative_x():
    pts = np.array([(-1.0, 1.0, 0.0)], dtype=geom.cartesian_type)
    out = geom.extrude_radially(pts, 1.0)
    assert np.hypot(out['x'][0], out['y'][0]) == pytest.approx(np.sqrt(2) + 1.0)
    assert out['x'][0] < 0
