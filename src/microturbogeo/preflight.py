"""Pre-flight geometry checks for a blade row. New in MicroTurboGeo.

Cheap numerical checks (no CAD kernel) that catch the usual reasons a CFD domain build fails or
produces a poor tip region: sections that do not span the annulus, clearance that is large or tiny
relative to span/chord, and unusual pitch-to-chord ratios.
"""
from __future__ import annotations

import numpy as np

from microturbogeo import geom

_UNIT_TO_M = {'metres': 1.0, 'millimetres': 1e-3}


def run(machine, stage, blade) -> dict:
    """Return a dictionary of derived quantities, warnings and errors for one blade row."""
    errors, warnings = [], []
    r_min, r_max = blade.radial_extent()
    z_le, z_te = blade.axial_extent()
    ax_chord = z_te - z_le
    pitch_deg = 360.0 / blade.n_blade

    rep = dict(units=machine.units, n_blade=blade.n_blade, pitch_deg=pitch_deg,
               section_r_min=r_min, section_r_max=r_max, axial_chord=ax_chord,
               n_sections=int(blade.ps_sections.shape[0]))

    ps = np.concatenate((blade.ps_sections, blade.ss_sections), axis=1)
    mid_i = ps.shape[0] // 2
    r_mid = float(np.mean(geom.radial_of(ps[mid_i])))
    pitch = r_mid * np.radians(pitch_deg)
    # true chord: distance from LE to TE point of the mid section along the cylindrical section
    le = blade.ps_sections[mid_i][0]
    te = blade.ps_sections[mid_i][-1]
    theta_le, theta_te = np.arctan2(le['y'], le['x']), np.arctan2(te['y'], te['x'])
    true_chord = float(np.hypot(ax_chord, r_mid * (theta_te - theta_le)))
    rep.update(mid_radius=r_mid, pitch_mid=pitch, true_chord_mid=true_chord,
               pitch_to_chord=pitch / true_chord, solidity=true_chord / pitch)
    if not 0.6 <= true_chord / pitch <= 2.5:
        warnings.append(f'solidity {true_chord / pitch:.2f} at mid-span is outside the usual 0.6-2.5 range; '
                        f'check blade count ({blade.n_blade}).')

    if np.any(np.diff(geom.radial_of(blade.ps_sections)[:, 0]) <= 0):
        errors.append('sections must be ordered hub -> tip (ascending radius)')

    ew = stage.endwalls
    if ew.type in ('fpd', 'parametric'):
        zz = np.linspace(z_le, z_te, 25)
        hub_r = geom.interp_radius(ew.hub, zz)
        sh_r = geom.interp_radius(ew.shroud, zz)
        span = float(np.mean(sh_r - hub_r))
        rep.update(hub_radius_mid=float(np.mean(hub_r)), shroud_radius_mid=float(np.mean(sh_r)), span=span,
                   hub_tip_ratio=float(np.mean(hub_r / sh_r)), aspect_ratio=span / true_chord)
        if r_min > hub_r.min():
            errors.append(f'sections start at r={r_min:g} but the hub reaches r={hub_r.min():g}: add a section at '
                          f'or below the hub radius so the blade penetrates the hub')
        if r_max < sh_r.max():
            errors.append(f'sections end at r={r_max:g} but the shroud reaches r={sh_r.max():g}: add a section at '
                          f'or above the shroud radius so the blade penetrates the casing')
        z_ew = (min(ew.hub['z'].min(), ew.shroud['z'].min()), max(ew.hub['z'].max(), ew.shroud['z'].max()))
        rep.update(z_inlet=float(z_ew[0]), z_outlet=float(z_ew[1]),
                   upstream_chords=float((z_le - z_ew[0]) / ax_chord), downstream_chords=float((z_ew[1] - z_te) / ax_chord))
        if z_ew[0] > z_le - 0.25 * ax_chord:
            warnings.append('domain inlet is closer than 0.25 axial chords to the leading edge')
        if z_ew[1] < z_te + 0.5 * ax_chord:
            warnings.append('domain outlet is closer than 0.5 axial chords to the trailing edge')

        if blade.tip_clearance > 0:
            c_le, c_te = blade.clearance_at(span)
            rep.update(tip_clearance_le=c_le, tip_clearance_te=c_te, clearance_pct_span=100 * max(c_le, c_te) / span,
                       clearance_pct_chord=100 * max(c_le, c_te) / true_chord)
            if max(c_le, c_te) > 0.05 * span:
                warnings.append('tip clearance exceeds 5 % of span - unusually large for a small turbomachine')
            if min(c_le, c_te) < 0.002 * span:
                warnings.append('tip clearance below 0.2 % of span - difficult to manufacture and mesh')
            rep['gap_cells_hint'] = ('resolve the gap with at least 8-10 cells in the radial direction; '
                                     f'target first-cell size <= {min(c_le, c_te) / 10:g} {machine.units}')
    else:
        if blade.tip_clearance > 0:
            errors.append("tip clearance requires 'fpd' or 'parametric' endwalls")

    if machine.rpm > 0:
        omega = machine.rpm * 2 * np.pi / 60.0
        scale = _UNIT_TO_M[machine.units]
        rep['rpm'] = machine.rpm
        rep['tip_speed_m_s'] = float(omega * r_max * scale) if 'shroud_radius_mid' not in rep else \
            float(omega * rep['shroud_radius_mid'] * scale)

    rep['errors'], rep['warnings'] = errors, warnings
    return rep


def format_report(rep: dict, stage_name: str, blade_name: str) -> str:
    """Human-readable text report."""
    u = rep['units']
    lines = [f'--- Pre-flight: {stage_name} / {blade_name} ---',
             f"blades: {rep['n_blade']}  (angular pitch {rep['pitch_deg']:.4f} deg), sections: {rep['n_sections']}",
             f"mid-span radius {rep['mid_radius']:g} {u}, pitch {rep['pitch_mid']:g} {u}, "
             f"axial chord {rep['axial_chord']:g} {u}, true chord {rep['true_chord_mid']:g} {u}",
             f"pitch/chord {rep['pitch_to_chord']:.3f}  (solidity {rep['solidity']:.3f})"]
    if 'span' in rep:
        lines.append(f"span {rep['span']:g} {u}, hub/tip ratio {rep['hub_tip_ratio']:.3f}, "
                     f"aspect ratio (span/true chord) {rep['aspect_ratio']:.2f}")
        lines.append(f"domain: {rep['upstream_chords']:.2f} axial chords upstream, "
                     f"{rep['downstream_chords']:.2f} downstream of the blade")
    if 'tip_clearance_le' in rep:
        lines.append(f"tip clearance {rep['tip_clearance_le']:g} (LE) -> {rep['tip_clearance_te']:g} (TE) {u} = "
                     f"{rep['clearance_pct_span']:.2f} % span, {rep['clearance_pct_chord']:.2f} % true chord")
        lines.append(f"meshing hint: {rep['gap_cells_hint']}")
    if 'tip_speed_m_s' in rep:
        lines.append(f"rotational speed {rep['rpm']:g} rpm -> tip speed {rep['tip_speed_m_s']:.1f} m/s")
    for w in rep['warnings']:
        lines.append(f'WARNING: {w}')
    for e in rep['errors']:
        lines.append(f'ERROR: {e}')
    return '\n'.join(lines)
