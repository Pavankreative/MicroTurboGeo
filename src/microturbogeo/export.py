"""Boundary-patch classification and solver-oriented export. New in MicroTurboGeo.

The fluid domain is a single periodic sector. Its faces are classified geometrically into named
boundary patches::

    inlet, outlet, hub, shroud, blade_ps_ss (blade wall), blade_tip (if a tip gap exists),
    periodic_a, periodic_b

``periodic_b`` is ``periodic_a`` rotated by -pitch about the machine axis (+Z).

Outputs (``export_all``):

* ``<name>.step``                 full domain solid
* ``<name>_patches.stl``          one ASCII STL, one named ``solid`` per patch -> OpenFOAM snappyHexMesh /
                                  surface-mesh tools (region names are preserved)
* ``<name>_<patch>.step``         one STEP per patch (for ANSYS Named Selections / SpaceClaim)
* ``<name>_patches.json``         manifest: patch -> face count, area, centroid; periodic rotation data; and
                                  ready-to-edit OpenFOAM / SU2 periodic boundary snippets

Always check the snippets against the solver version you use; they are starting points, not a guarantee.
"""
from __future__ import annotations

import json
import os
import re
from typing import Dict, List

import numpy as np
import cadquery as cq

from microturbogeo import geom

PATCH_ORDER = ['inlet', 'outlet', 'hub', 'shroud', 'blade', 'blade_tip', 'periodic_a', 'periodic_b']


def _rot_z(p: np.ndarray, ang: float) -> np.ndarray:
    c, s = np.cos(ang), np.sin(ang)
    return np.array([c * p[0] - s * p[1], s * p[0] + c * p[1], p[2]])


def _sample_points(face: cq.Face, tol: float, n: int = 15) -> np.ndarray:
    """Centroids of n evenly spread tessellation triangles (interior points of the face)."""
    verts, tris = face.tessellate(tol, 0.2)
    v = np.array([[p.x, p.y, p.z] for p in verts])
    idx = np.unique(np.linspace(0, len(tris) - 1, min(n, len(tris))).astype(int))
    return np.array([v[list(tris[i])].mean(axis=0) for i in idx])


def classify_faces(creator, units: str = 'metres') -> Dict[str, List[cq.Face]]:
    """Assign every face of ``creator.domain`` to a named patch."""
    faces = creator.domain.faces().vals()
    length_scale = 1.0 if units == 'metres' else 1e3
    tol = 1e-6 * length_scale
    pitch = creator.blade_def.pitch_angle_rad
    blade = creator.blade
    ew = creator.endwalls

    centers = [np.array(f.Center().toTuple()) for f in faces]
    areas = np.array([f.Area() for f in faces])
    patches: Dict[str, List[cq.Face]] = {k: [] for k in PATCH_ORDER}
    assigned = [None] * len(faces)

    # periodic pair: face B is face A rotated by -pitch (same area, rotated centroid)
    for i in range(len(faces)):
        if assigned[i]:
            continue
        for j in range(len(faces)):
            if i == j or assigned[j]:
                continue
            if abs(areas[i] - areas[j]) > 1e-3 * areas[i]:
                continue
            if np.linalg.norm(_rot_z(centers[i], -pitch) - centers[j]) < 1e-5 * length_scale:
                assigned[i], assigned[j] = 'periodic_a', 'periodic_b'
                break

    bbs = [f.BoundingBox() for f in faces]
    z_lo, z_hi = float(min(bb.zmin for bb in bbs)), float(max(bb.zmax for bb in bbs))
    prof = creator.tip_gap_profile()
    tess_tol = 2e-6 * length_scale
    near = 10 * tess_tol
    for i, f in enumerate(faces):
        if assigned[i]:
            continue
        bb = bbs[i]
        # inlet / outlet: faces lying in a plane perpendicular to the axis at the ends of the domain
        if bb.zmax - bb.zmin < 10 * tol:
            if abs(bb.zmin - z_lo) < 10 * tol:
                assigned[i] = 'inlet'
                continue
            if abs(bb.zmax - z_hi) < 10 * tol:
                assigned[i] = 'outlet'
                continue
        pts = _sample_points(f, tess_tol)
        r = np.hypot(pts[:, 0], pts[:, 1])
        z = pts[:, 2]
        # solid classifier with a tolerance: points on (or within `near` of) the blade boundary count as on-blade.
        # (Vertex.distance() against B-spline loft faces proved unreliable, so it is not used here.)
        on_blade = np.mean([blade.isInside(cq.Vector(*p), near) for p in pts])
        if on_blade > 0.6:
            assigned[i] = 'blade'
            if prof is not None:
                z_le, z_te, c_le, c_te = prof
                frac = np.clip((z - z_le) / (z_te - z_le), 0.0, 1.0)
                r_off = geom.interp_radius(ew.shroud, z) - (c_le + (c_te - c_le) * frac)
                if np.mean(np.abs(r - r_off) < near) > 0.6:
                    assigned[i] = 'blade_tip'
            continue
        if ew.type in ('fpd', 'parametric'):
            dh = np.abs(r - geom.interp_radius(ew.hub, z))
            ds = np.abs(r - geom.interp_radius(ew.shroud, z))
            assigned[i] = 'hub' if np.mean(dh < ds) > 0.5 else 'shroud'
        else:  # CAD endwalls: split about the mid radius of the domain
            r_all = [np.hypot(c[0], c[1]) for c in centers]
            assigned[i] = 'hub' if np.mean(r) < 0.5 * (min(r_all) + max(r_all)) else 'shroud'

    for f, name in zip(faces, assigned):
        patches[name].append(f)
    return {k: v for k, v in patches.items() if v}


def _ascii_stl_for(faces: List[cq.Face], name: str, tolerance: float, angular_tol: float = 0.1) -> str:
    """Tessellate faces and return an ASCII STL solid block named ``name``."""
    out = [f'solid {name}']
    for f in faces:
        verts, tris = f.tessellate(tolerance, angular_tol)
        v = np.array([[p.x, p.y, p.z] for p in verts])
        for a, b, c in tris:
            n = np.cross(v[b] - v[a], v[c] - v[a])
            norm = np.linalg.norm(n)
            n = n / norm if norm > 0 else n
            out.append(f' facet normal {n[0]:.6e} {n[1]:.6e} {n[2]:.6e}\n  outer loop')
            for k in (a, b, c):
                out.append(f'   vertex {v[k][0]:.9e} {v[k][1]:.9e} {v[k][2]:.9e}')
            out.append('  endloop\n endfacet')
    out.append(f'endsolid {name}')
    return '\n'.join(out) + '\n'


def manifest(creator, machine, patches: Dict[str, List[cq.Face]]) -> dict:
    """Build the JSON-serialisable patch manifest."""
    pitch = float(creator.blade_def.pitch_angle_rad)
    info = {}
    for name, fl in patches.items():
        c = np.array([f.Center().toTuple() for f in fl])
        a = np.array([f.Area() for f in fl])
        info[name] = dict(n_faces=len(fl), area=float(a.sum()),
                          centroid=[float(x) for x in (c * a[:, None]).sum(0) / a.sum()])
    m = dict(units=machine.units, n_blade=int(creator.blade_def.n_blade), axis=[[0, 0, 0], [0, 0, 1]],
             sector_angle_rad=pitch, sector_angle_deg=float(np.degrees(pitch)),
             tip_clearance_le_te=list(getattr(creator, 'tip_clearance_applied', (0.0, 0.0))),
             patches=info)
    if 'periodic_a' in patches and 'periodic_b' in patches:
        m['periodic'] = dict(
            note='periodic_b = periodic_a rotated by -sector_angle about +Z through the origin',
            openfoam_boundary_snippet=(
                'periodic_a { type cyclic; neighbourPatch periodic_b; transform rotational; '
                'rotationAxis (0 0 1); rotationCentre (0 0 0); }\n'
                'periodic_b { type cyclic; neighbourPatch periodic_a; transform rotational; '
                'rotationAxis (0 0 1); rotationCentre (0 0 0); }'),
            su2_snippet=(f'MARKER_PERIODIC= ( periodic_a, periodic_b, 0.0, 0.0, 0.0, 0.0, 0.0, '
                         f'{pitch:.12f}, 0.0, 0.0, 0.0 )   % check sign and angle units for your SU2 version'),
            ansys_note='Import the per-patch STEP files (or select faces by the centroids above) and create Named '
                       'Selections; define periodic_a/periodic_b as a Rotational periodic pair about Z.')
    return m


def _safe(name: str) -> str:
    return re.sub(r'[^A-Za-z0-9_.-]+', '_', name)


def export_all(creator, machine, step_fname: str, formats=('step',), stl_tolerance=None) -> List[str]:
    """Write the requested export formats next to ``step_fname``. Returns the list of written files."""
    base = re.sub(r'\.step$', '', str(step_fname))
    written = []
    creator.export('domain', f'{base}.step')
    written.append(f'{base}.step')

    patches = classify_faces(creator, machine.units)
    man = manifest(creator, machine, patches)
    with open(f'{base}_patches.json', 'w') as fp:
        json.dump(man, fp, indent=2)
    written.append(f'{base}_patches.json')

    if 'stl' in formats:
        if stl_tolerance is None:
            stl_tolerance = creator.annulus_span() / 2000.0 if creator.endwalls.type != 'step' else \
                (1e-5 if machine.units == 'metres' else 1e-2)
        names = {'blade': 'blade', 'blade_tip': 'blade_tip'}
        with open(f'{base}_patches.stl', 'w') as fp:
            for pname, fl in patches.items():
                fp.write(_ascii_stl_for(fl, names.get(pname, pname), stl_tolerance))
        written.append(f'{base}_patches.stl')

    if 'patch-step' in formats:
        for pname, fl in patches.items():
            fn = f'{base}_{_safe(pname)}.step'
            cq.exporters.export(cq.Workplane('XY').add(cq.Compound.makeCompound(fl)), fn)
            written.append(fn)

    for w in written:
        print('wrote', w)
    return written
