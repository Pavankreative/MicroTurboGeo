"""Parametric blade-section generator for compact axial turbine and compressor rows.

New in MicroTurboGeo. Produces pressure-surface (PS) and suction-surface (SS) curves on cylindrical
sections at user-chosen radii, in the formatted point data (.fpd) layout read by :mod:`microturbogeo.geom`.

Section model (per radius, unwrapped to the z - r*theta plane):

* mean line: linearly varying tangent of metal angle (parabolic camber) between ``beta_in`` and ``beta_out``
* thickness: NACA 4-digit style distribution (round leading edge) plus a linear term that gives a finite
  trailing-edge thickness ``t_te``; the blunt trailing edge is closed by a semicircular arc
* the SS curve runs LE -> suction corner -> TE arc -> pressure corner so that PS and SS share exactly the
  same first and last point (required by the solid construction in :mod:`microturbogeo.cad`)

Sign convention (matches the loader and the periodic-domain code): the pressure side lies on the +theta
side of the mean line, i.e. the camber line is concave towards +theta and ``tan(beta_out) > tan(beta_in)``.
Use :func:`mirror_theta` for the opposite hand.

These sections are a geometry-generation aid for building CFD domains quickly. They are NOT an
aerodynamically optimised design; supply your own section data for real design work.
"""
from __future__ import annotations

import pathlib
from typing import Sequence, Tuple

import numpy as np

from microturbogeo.geom import cartesian_type

_NACA = (0.2969, -0.1260, -0.3516, 0.2843, -0.1036)


def _half_thickness(s: np.ndarray, chord: float, t_ratio: float, t_te: float) -> np.ndarray:
    a0, a1, a2, a3, a4 = _NACA
    yt = 5.0 * t_ratio * (a0 * np.sqrt(s) + a1 * s + a2 * s ** 2 + a3 * s ** 3 + a4 * s ** 4)
    return chord * yt + s * t_te / 2.0


def build_section_2d(chord: float, beta_in_deg: float, beta_out_deg: float, t_ratio: float, t_te: float,
                     n_pts: int = 60, n_te_arc: int = 9) -> Tuple[np.ndarray, np.ndarray]:
    """Build one unwrapped section.

    Returns:
        (ps, ss): each an (N, 2) array of (z, r*theta) points running LE -> TE.
    """
    if chord <= 0 or not 0 < t_ratio < 0.5 or t_te < 0:
        raise ValueError('chord must be >0, 0 < t_ratio < 0.5 and t_te >= 0')
    u = np.linspace(0.0, 1.0, n_pts)
    s = 0.5 * (1.0 - np.cos(np.pi * u))                    # cosine clustering at LE and TE
    tb1, tb2 = np.tan(np.radians(beta_in_deg)), np.tan(np.radians(beta_out_deg))
    z = chord * s
    m = chord * (s * tb1 + 0.5 * s ** 2 * (tb2 - tb1))      # mean line, r*theta
    phi = np.arctan(tb1 + (tb2 - tb1) * s)                  # local tangent angle
    h = _half_thickness(s, chord, t_ratio, t_te)
    nz, nt = -np.sin(phi), np.cos(phi)                      # unit normal towards +theta (pressure side)

    ps = np.column_stack((z + h * nz, m + h * nt))
    ss = np.column_stack((z - h * nz, m - h * nt))
    ps[0] = ss[0] = (0.0, 0.0)                              # exact shared leading-edge point

    if t_te > 0:
        r_te = t_te / 2.0
        c = np.array([z[-1], m[-1]])
        T = np.array([np.cos(phi[-1]), np.sin(phi[-1])])
        N = np.array([-np.sin(phi[-1]), np.cos(phi[-1])])
        t = np.linspace(0.0, np.pi, n_te_arc + 2)[1:-1]
        arc = c[None, :] - r_te * np.outer(np.cos(t), N) + r_te * np.outer(np.sin(t), T)
        ss = np.vstack((ss, arc, ps[-1]))                   # SS ends on the pressure corner
    else:
        ps[-1] = ss[-1]
    return ps, ss


def make_blade_sections(radii: Sequence[float], chord: float, beta_in_deg: Sequence[float],
                        beta_out_deg: Sequence[float], t_ratio: Sequence[float], t_te: float,
                        n_pts: int = 60, z_le: float = 0.0) -> Tuple[np.ndarray, np.ndarray]:
    """Stack cylindrical sections (leading edge on the theta = 0 ray) from hub to tip.

    Args:
        radii: section radii, ascending (hub first). Extend them slightly past the hub and shroud
               radii so the blade fully penetrates both endwalls.
        chord: axial chord (constant over span)
        beta_in_deg, beta_out_deg: metal angles from axial at each radius
        t_ratio: maximum thickness / chord at each radius
        t_te: trailing-edge thickness (absolute, model units)
        z_le: axial position of the leading edge

    Returns:
        ps, ss: arrays of shape (n_sections, n_points) with dtype geom.cartesian_type
    """
    radii = np.asarray(radii, dtype=float)
    if np.any(np.diff(radii) <= 0):
        raise ValueError('radii must be strictly ascending (hub first)')
    n = len(radii)
    if not (len(beta_in_deg) == len(beta_out_deg) == len(t_ratio) == n):
        raise ValueError('radii, beta_in_deg, beta_out_deg and t_ratio must have the same length')

    ps_all, ss_all = [], []
    for i, r in enumerate(radii):
        ps2, ss2 = build_section_2d(chord, beta_in_deg[i], beta_out_deg[i], t_ratio[i], t_te, n_pts)
        ps_all.append(_wrap(ps2, r, z_le))
        ss_all.append(_wrap(ss2, r, z_le))
    return np.array(ps_all), np.array(ss_all)


def _wrap(curve2d: np.ndarray, r: float, z_le: float) -> np.ndarray:
    out = np.empty(curve2d.shape[0], dtype=cartesian_type)
    theta = curve2d[:, 1] / r
    out['x'] = r * np.cos(theta)
    out['y'] = r * np.sin(theta)
    out['z'] = curve2d[:, 0] + z_le
    return out


def mirror_theta(sections: np.ndarray) -> np.ndarray:
    """Mirror sections about the y=0 plane (reverses the hand of the blade row)."""
    out = sections.copy()
    out['y'] = -out['y']
    return out


def free_vortex_angle(r: np.ndarray, r_ref: float, angle_ref_deg: float) -> np.ndarray:
    """Flow/metal angle distribution (deg from axial) for a free-vortex swirl, r*tan(alpha) = const."""
    k = r_ref * np.tan(np.radians(angle_ref_deg))
    return np.degrees(np.arctan(k / np.asarray(r)))


def write_fpd(fname, curves: np.ndarray) -> None:
    """Write curves to the formatted point data layout: header 'n_points n_curves', then x y z rows."""
    curves = np.atleast_2d(curves)
    n_curve, n_pts = curves.shape
    with open(pathlib.Path(fname), 'w') as f:
        f.write(f'{n_pts} {n_curve}\n')
        for c in curves:
            for p in c:
                f.write(f"{p['x']:.15e} {p['y']:.15e} {p['z']:.15e}\n")
