"""End-to-end CAD test: build the micro turbine rotor domain and verify topology, volume, tip gap and export."""
import json
import pathlib
import warnings

import numpy as np
import pytest

pytestmark = pytest.mark.slow
cq = pytest.importorskip('cadquery')

from microturbogeo import export, geom                      # noqa: E402
from microturbogeo.cad import DomainCreator                 # noqa: E402
from microturbogeo.machine import Machine                   # noqa: E402

CASE = pathlib.Path(__file__).parent.parent / 'examples' / 'micro_turbojet_turbine_rotor' / 'micro_turbojet_turbine_rotor.toml'


@pytest.fixture(scope='module')
def built():
    warnings.simplefilter('ignore')
    m = Machine.from_config_file(CASE)
    st = m.stages[0]
    c = DomainCreator(st.blades[0], st.endwalls, m.units, m.axis)
    c.create_domain()
    return m, st, c


def test_domain_is_one_valid_solid(built):
    _, _, c = built
    assert c.domain.val().isValid()
    assert len(c.domain.solids().vals()) == 1


def test_domain_volume_is_sector_minus_blade(built):
    _, _, c = built
    sector = cq.Workplane('XY').add(c.per).intersect(c.cad_endwalls)
    blade_inside = cq.Workplane('XY').add(c.blade).intersect(sector)
    expected = sector.val().Volume() - blade_inside.val().Volume()
    assert c.domain.val().Volume() == pytest.approx(expected, rel=2e-3)


def test_sector_volume_matches_annulus_geometry(built):
    """Sector volume ~ (2*pi/N) * integral of the annulus cross-section about the axis."""
    m, st, c = built
    z = np.linspace(-0.004, 0.016, 2001)
    rh, rs = geom.interp_radius(st.endwalls.hub, z), geom.interp_radius(st.endwalls.shroud, z)
    # cosine-blended spline vs piecewise-linear interpolation of the same points: allow 1 %
    ann = np.trapezoid(0.5 * (rs ** 2 - rh ** 2), z)
    expected = ann * (2 * np.pi / st.blades[0].n_blade)
    sector = cq.Workplane('XY').add(c.per).intersect(c.cad_endwalls).val().Volume()
    assert sector == pytest.approx(expected, rel=0.01)


def test_patches(built):
    m, st, c = built
    patches = export.classify_faces(c, m.units)
    assert set(patches) == {'inlet', 'outlet', 'hub', 'shroud', 'blade', 'blade_tip', 'periodic_a', 'periodic_b'}
    assert sum(len(v) for v in patches.values()) == len(c.domain.faces().vals())
    pa, pb = patches['periodic_a'][0], patches['periodic_b'][0]
    assert pa.Area() == pytest.approx(pb.Area(), rel=1e-6)
    pitch = st.blades[0].pitch_angle_rad
    ca, cb = np.array(pa.Center().toTuple()), np.array(pb.Center().toTuple())
    rot = np.array([[np.cos(-pitch), -np.sin(-pitch)], [np.sin(-pitch), np.cos(-pitch)]]) @ ca[:2]
    np.testing.assert_allclose(rot, cb[:2], atol=1e-8)


def test_tip_gap_matches_requested_taper(built):
    """Every point of the tip face must sit between c_te and c_le below the shroud."""
    m, st, c = built
    tip = export.classify_faces(c, m.units)['blade_tip'][0]
    pts = export._sample_points(tip, 2e-6, n=60)
    r, z = np.hypot(pts[:, 0], pts[:, 1]), pts[:, 2]
    gap = geom.interp_radius(st.endwalls.shroud, z) - r
    c_le, c_te = 0.0003, 0.00025
    assert gap.min() > c_te * 0.97 and gap.max() < c_le * 1.03
    # and it tapers: gap near the LE is larger than near the TE
    assert gap[np.argmin(z)] > gap[np.argmax(z)]


def test_trim_removed_blade_material(built):
    _, _, c = built
    assert c.blade.Volume() < c.blade_untrimmed.Volume()
    assert c.tip_clearance_applied == (0.0003, 0.00025)


def test_export_all(built, tmp_path):
    m, st, c = built
    out = export.export_all(c, m, str(tmp_path / 'rotor.step'), formats=('step', 'stl', 'patch-step'))
    names = {pathlib.Path(p).name for p in out}
    assert {'rotor.step', 'rotor_patches.json', 'rotor_patches.stl', 'rotor_blade_tip.step'} <= names
    man = json.loads((tmp_path / 'rotor_patches.json').read_text())
    assert man['n_blade'] == 25 and 'periodic' in man
    stl = (tmp_path / 'rotor_patches.stl').read_text()
    for patch in man['patches']:
        assert f'solid {patch}\n' in stl
    reread = cq.importers.importStep(str(tmp_path / 'rotor.step'))
    assert reread.val().Volume() == pytest.approx(c.domain.val().Volume(), rel=1e-3)
