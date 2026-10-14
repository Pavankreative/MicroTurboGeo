"""Pre-flight checks on the shipped example cases (no CAD kernel needed)."""
import pathlib

import pytest

from microturbogeo import preflight
from microturbogeo.machine import Machine

EX = pathlib.Path(__file__).parent.parent / 'examples'
CASES = ['micro_turbojet_ngv', 'micro_turbojet_turbine_rotor', 'micro_turbojet_axial_compressor_rotor']


def _report(case):
    m = Machine.from_config_file(EX / case / f'{case}.toml')
    st = m.stages[0]
    return m, st, st.blades[0], preflight.run(m, st, st.blades[0])


@pytest.mark.parametrize('case', CASES)
def test_examples_pass_preflight(case):
    _, _, _, rep = _report(case)
    assert rep['errors'] == []
    assert 0.8 < rep['solidity'] < 2.0


def test_rotor_tip_speed_and_clearance_report():
    _, _, _, rep = _report('micro_turbojet_turbine_rotor')
    assert rep['tip_speed_m_s'] == pytest.approx(110000 * 2 * 3.14159265 / 60 * 0.0408, rel=0.02)
    assert rep['tip_clearance_le'] == pytest.approx(0.0003)
    assert rep['tip_clearance_te'] == pytest.approx(0.00025)


def test_percent_span_clearance_resolved():
    _, _, b, rep = _report('micro_turbojet_axial_compressor_rotor')
    assert rep['clearance_pct_span'] == pytest.approx(1.2, rel=1e-6)
    assert rep['tip_clearance_le'] == pytest.approx(0.012 * rep['span'])


def test_sections_not_reaching_shroud_is_an_error():
    m, st, b, _ = _report('micro_turbojet_ngv')
    st.endwalls.shroud = geom_shift(st.endwalls.shroud, 0.003)
    rep = preflight.run(m, st, b)
    assert any('casing' in e for e in rep['errors'])


def geom_shift(curve, dr):
    from microturbogeo import geom
    return geom.radial_offset_curve(curve, dr)


def test_clearance_with_step_endwalls_is_an_error():
    m, st, b, _ = _report('micro_turbojet_turbine_rotor')
    st.endwalls.type = 'step'
    rep = preflight.run(m, st, b)
    assert any('tip clearance' in e for e in rep['errors'])


def test_report_formats():
    m, st, b, rep = _report('micro_turbojet_turbine_rotor')
    text = preflight.format_report(rep, st.name, b.name)
    assert 'tip speed' in text and 'meshing hint' in text
