"""Regenerate the three micro-turbojet example cases (section files and endwall files).

Run from the repository root:

    python examples/generate_examples.py

All dimensions are in metres. The cases are *illustrative reference geometries* sized for a
~100-150 N class micro-turbojet (core diameter ~ 80 mm, 90-110 krpm). They exercise the geometry tool;
they are not aerodynamically optimised designs and no performance claim is made for them.
"""
import pathlib

import numpy as np

from microturbogeo import geom, sections

HERE = pathlib.Path(__file__).parent


def _toml(path: pathlib.Path, text: str) -> None:
    path.write_text(text)


def ngv() -> None:
    """Axial nozzle guide vane row, free-vortex outlet swirl, constant-radius annulus."""
    out = HERE / 'micro_turbojet_ngv'
    out.mkdir(exist_ok=True)
    radii = np.array([0.0295, 0.0345, 0.0395])            # 0.5 mm beyond hub (0.030) and shroud (0.039)
    beta_out = sections.free_vortex_angle(radii, 0.0345, 67.0)
    ps, ss = sections.make_blade_sections(radii, chord=0.012, beta_in_deg=[0, 0, 0], beta_out_deg=beta_out,
                                          t_ratio=[0.20, 0.16, 0.13], t_te=0.00025, n_pts=60)
    sections.write_fpd(out / 'ngv_ps.fpd', ps)
    sections.write_fpd(out / 'ngv_ss.fpd', ss)
    _toml(out / 'micro_turbojet_ngv.toml', '''# Axial NGV row, micro-turbojet (~100-150 N class). Units: metres. Sign convention: pressure side on +theta.
[machine]
name = "Micro-turbojet NGV"
units = 'metres'
axis = [ [0.0,0.0,0.0] , [0.0,0.0,1.0] ]
rpm = 0.0                       # stator

[[stage]]
name = 'stage_1'

[[stage.endwall]]
type = 'parametric'
z_inlet = -0.006                # 0.5 axial chord upstream of the leading edge
z_outlet = 0.024                # 1.0 axial chord downstream of the trailing edge
hub_radius = [0.0300, 0.0300]
shroud_radius = [0.0390, 0.0390]
n_points = 11
blend = 'linear'

[[stage.blade_section]]
name = 'ngv'
n_blade = 17                    # pitch at mid-span = 12.7 mm; alternatively: pitch_deg = 21.17647...
ps_section_fname = 'ngv_ps.fpd'
ss_section_fname = 'ngv_ss.fpd'
''')


def turbine_rotor() -> None:
    """Axial turbine rotor with tapered tip clearance and a gently diverging annulus."""
    out = HERE / 'micro_turbojet_turbine_rotor'
    out.mkdir(exist_ok=True)
    radii = np.array([0.0285, 0.0352, 0.0420])             # hub min 0.0295, shroud max 0.0412
    k = np.tan(np.radians(57.0)) / 0.0352                   # tan(beta_out) ~ r (blade speed ~ r)
    beta_out = np.degrees(np.arctan(k * radii))
    ps, ss = sections.make_blade_sections(radii, chord=0.008, beta_in_deg=[-20, -28, -32], beta_out_deg=beta_out,
                                          t_ratio=[0.22, 0.17, 0.10], t_te=0.0002, n_pts=60)
    sections.write_fpd(out / 'rotor_ps.fpd', ps)
    sections.write_fpd(out / 'rotor_ss.fpd', ss)
    _toml(out / 'micro_turbojet_turbine_rotor.toml', '''# Axial turbine rotor, micro-turbojet. Units: metres.
# Tip speed at 110 krpm and r = 41.2 mm is ~475 m/s. Tip clearance tapers 0.30 mm (LE) -> 0.25 mm (TE).
[machine]
name = "Micro-turbojet turbine rotor"
units = 'metres'
axis = [ [0.0,0.0,0.0] , [0.0,0.0,1.0] ]
rpm = 110000.0

[[stage]]
name = 'stage_1'

[[stage.endwall]]
type = 'parametric'
z_inlet = -0.004
z_outlet = 0.016
hub_radius = [0.0300, 0.0295]
shroud_radius = [0.0405, 0.0412]
n_points = 21
blend = 'cosine'

[[stage.blade_section]]
name = 'rotor'
n_blade = 25
tip_clearance = 0.00030         # absolute, metres (at the leading edge)
tip_clearance_te = 0.00025      # optional, at the trailing edge
tip_clearance_mode = 'absolute'
ps_section_fname = 'rotor_ps.fpd'
ss_section_fname = 'rotor_ss.fpd'
''')


def compressor_rotor() -> None:
    """Axial compressor inducer rotor, endwalls supplied as fpd curves, clearance as % span."""
    out = HERE / 'micro_turbojet_axial_compressor_rotor'
    out.mkdir(exist_ok=True)
    radii = np.array([0.0115, 0.0210, 0.0315])             # hub min 0.0125, shroud 0.0300
    k = np.tan(np.radians(70.0)) / 0.0300                   # tan(beta_in) ~ r
    beta_in = -np.degrees(np.arctan(k * radii))
    turning = np.array([16.0, 13.0, 11.0])
    ps, ss = sections.make_blade_sections(radii, chord=0.010, beta_in_deg=beta_in, beta_out_deg=beta_in + turning,
                                          t_ratio=[0.10, 0.07, 0.04], t_te=0.00015, n_pts=60)
    sections.write_fpd(out / 'inducer_ps.fpd', ps)
    sections.write_fpd(out / 'inducer_ss.fpd', ss)
    hub = geom.make_parametric_meridional(0.0125, 0.0140, -0.006, 0.020, 21, 'cosine')
    shroud = geom.make_parametric_meridional(0.0300, 0.0300, -0.006, 0.020, 21, 'linear')
    sections.write_fpd(out / 'hub.fpd', hub)
    sections.write_fpd(out / 'shroud.fpd', shroud)
    _toml(out / 'micro_turbojet_axial_compressor_rotor.toml', '''# Axial compressor inducer rotor, micro-turbojet, hub/shroud supplied as .fpd meridional curves.
# Tip speed at 90 krpm and r = 30 mm is ~283 m/s. Tip clearance 1.2 % of span (~0.2 mm).
[machine]
name = "Micro-turbojet axial compressor rotor"
units = 'metres'
axis = [ [0.0,0.0,0.0] , [0.0,0.0,1.0] ]
rpm = 90000.0

[[stage]]
name = 'stage_1'

[[stage.endwall]]
type = 'fpd'
hub_fname = 'hub.fpd'
shroud_fname = 'shroud.fpd'

[[stage.blade_section]]
name = 'inducer'
n_blade = 12
tip_clearance = 1.2
tip_clearance_mode = 'percent_span'
ps_section_fname = 'inducer_ps.fpd'
ss_section_fname = 'inducer_ss.fpd'
''')


if __name__ == '__main__':
    ngv()
    turbine_rotor()
    compressor_rotor()
    print('examples written to', HERE)
