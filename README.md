# MicroTurboGeo

**CFD-ready fluid-domain geometry for micro-turbojet and UAV axial turbomachinery blade rows.**

MicroTurboGeo turns blade section definitions (pressure- and suction-surface curves at a few radii) plus a hub and
shroud definition into a single periodic fluid sector: one blade passage, bounded by hub, shroud, inlet, outlet,
a pair of rotationally periodic faces, and the blade wall, with an optional, explicitly modelled **tip-clearance gap**.
The result is exported as STEP, plus named-patch STL/STEP and a JSON manifest, ready for meshing in ANSYS,
OpenFOAM or SU2 workflows.

It targets small, high-speed machines (core diameters of roughly 60-100 mm, 90-150 krpm), where the tip gap is a
tenth of a millimetre and a poorly defined endwall or a coarse geometric tolerance is a large fraction of the
passage height.

## Contents
1. [Why a micro-scale tool](#why-a-micro-scale-tool)
2. [Installation](#installation)
3. [Quick start](#quick-start)
4. [Case file reference](#case-file-reference)
5. [Geometry conventions](#geometry-conventions)
6. [Tip clearance](#tip-clearance)
7. [Hub and shroud definition](#hub-and-shroud-definition)
8. [Periodic faces and boundary patches](#periodic-faces-and-boundary-patches)
9. [Pre-flight report](#pre-flight-report)
10. [Example cases](#example-cases)
11. [Meshing and solver notes](#meshing-and-solver-notes)
12. [Limitations](#limitations)
13. [Testing](#testing)
14. [What I contributed](#what-i-contributed)
15. [Validation](#validation)
16. [Attribution and license](#attribution-and-license)

## Why a micro-scale tool

| Micro-turbojet feature | How MicroTurboGeo handles it |
|---|---|
| Tip gap of 0.1-0.3 mm is 1-3 % of span | Gap is part of the geometry: blade is trimmed against an offset shroud; constant or LE-to-TE tapered; absolute or % span |
| Short blades, small chord (5-15 mm) | Geometric tolerances scale with the domain size, not fixed at 1 mm |
| Annulus changes radius through the row | Parametric hub/shroud (linear or smooth cosine blend), or your own meridional curves |
| Few blades, large pitch angle | `n_blade` or `pitch_deg`; pre-flight reports pitch, pitch/chord and solidity |
| Very high rotational speed | Optional `rpm` gives tip speed in the report |
| Mesh needs to resolve the gap | Pre-flight prints a first-cell size hint; named `blade_tip` patch |

## Installation

Requires Python 3.9+ and CadQuery (OpenCASCADE).

```bash
pip install -e .
```

## Quick start

```bash
# 1. Check a case without building CAD (instant)
microturbogeo examples/micro_turbojet_turbine_rotor/micro_turbojet_turbine_rotor.toml --check

# 2. Build the domain and export STEP + multi-region STL + patch manifest
microturbogeo examples/micro_turbojet_turbine_rotor/micro_turbojet_turbine_rotor.toml --formats step stl
```

Outputs are written next to the case file, named `<case>-<stage>-<blade>`:

| File | Contents |
|---|---|
| `*.step` | Full fluid domain solid |
| `*_patches.stl` | One ASCII STL, one named `solid` per patch (`--formats stl`) |
| `*_<patch>.step` | One STEP per patch (`--formats patch-step`) |
| `*_patches.json` | Patch areas and centroids, sector angle, periodic rotation data, OpenFOAM/SU2 snippets |

Build time depends on the case; the three shipped examples took about 0.5 to 2 minutes each (28 s, 98 s and 123 s) on a cloud
container, measured 2026-10-15. Run `--check` first, since it catches most input mistakes in well under a second.

## Case file reference

A case is a TOML file. File names are relative to the case file. Units are metres or millimetres
(`units`); the machine axis must be +Z through the origin.

```toml
[machine]
name = "Micro-turbojet turbine rotor"
units = 'metres'
axis = [ [0.0,0.0,0.0] , [0.0,0.0,1.0] ]
rpm = 110000.0                  # optional, only used for the tip-speed report

[[stage]]
name = 'stage_1'

[[stage.endwall]]               # see "Hub and shroud definition"
type = 'parametric'
z_inlet = -0.004
z_outlet = 0.016
hub_radius = [0.0300, 0.0295]   # at inlet, at outlet
shroud_radius = [0.0405, 0.0412]
n_points = 21                   # optional, default 21
blend = 'cosine'                # 'cosine' (default) or 'linear'

[[stage.blade_section]]
name = 'rotor'
n_blade = 25                    # or pitch_deg = 14.4 (must divide 360 exactly); both may be given if consistent
tip_clearance = 0.00030         # at the leading edge
tip_clearance_te = 0.00025      # optional; default = same as leading edge
tip_clearance_mode = 'absolute' # 'absolute' (model units) or 'percent_span'
ps_section_fname = 'rotor_ps.fpd'
ss_section_fname = 'rotor_ss.fpd'
```

### Section files (`.fpd`)

Plain text: a header `n_points n_curves`, then `x y z` rows for each curve in turn. Sections run from hub to tip
(ascending radius), each curve from leading edge to trailing edge. Requirements:

* PS and SS curves of a section must share **exactly** the same first point (leading edge) and last point
  (trailing edge). A blunt trailing edge is handled by running the SS curve around the trailing-edge arc so it
  ends on the pressure-side corner.
* Sections must extend past the hub and shroud (and past the tip-clearance offset surface) so the blade fully
  penetrates both endwalls. Pre-flight reports an error if they do not.

`microturbogeo.sections` can generate sections (camber, NACA-style thickness, blunt trailing edge) and write
`.fpd` files; see `examples/generate_examples.py`. These generated sections are a convenience for building CFD
domains quickly and are not aerodynamically optimised.

## Geometry conventions

* Axial direction is +Z; the angular coordinate is theta about Z.
* **The pressure side lies on the +theta side of the mean line** (the camber line is concave towards +theta).
  For the opposite hand, use `sections.mirror_theta`.
* The periodic sector is generated from a midline surface between the blade and its +theta neighbour. Face
  `periodic_b` is `periodic_a` rotated by **-pitch** about Z.
* Blade, pitch and endwall radii must be consistent: sections of neighbouring blades must not overlap.

## Tip clearance

Set `tip_clearance` on the blade row. The blade is trimmed against a copy of the shroud offset radially towards
the axis, and the resulting gap stays in the fluid domain, so tip leakage is part of the CFD geometry.

| Option | Meaning |
|---|---|
| `tip_clearance_mode = 'absolute'` | value in model units |
| `tip_clearance_mode = 'percent_span'` | percentage of the hub-to-shroud height at mid-chord |
| `tip_clearance_te` | clearance at the trailing edge; the gap varies linearly between the blade leading and trailing edges, and is constant upstream and downstream |

Requires `parametric` or `fpd` endwalls (the shroud curve is offset radially). The exported `blade_tip` patch is
the trimmed tip surface. The gap is verified in the test suite: every point of the tip face lies between the
requested trailing-edge and leading-edge clearances below the shroud.

## Hub and shroud definition

Three endwall types:

* **`parametric`**: end radii and axial limits (above). `blend = 'linear'` gives a cone, `'cosine'` a smooth
  S-shaped transition with zero slope at inlet and outlet. Validated: the shroud must exceed the hub everywhere.
* **`fpd`**: your own meridional curves for `hub_fname` and `shroud_fname` (points in the y = 0 plane).
* **`step`**: a CAD solid of the endwall region (supports cavities). Tip clearance is not available with this type.

The inlet and outlet planes are the axial limits of the endwall curves. Pre-flight reports how many axial chords
the domain extends upstream and downstream of the blade and warns when they are short (less than 0.25 chord
upstream or 0.5 chord downstream).

## Periodic faces and boundary patches

The domain is classified into named patches:

`inlet`, `outlet`, `hub`, `shroud`, `blade`, `blade_tip` (when a gap exists), `periodic_a`, `periodic_b`

Classification is geometric (see `src/microturbogeo/export.py`): the periodic pair is found as two equal-area faces
related by a rotation of one pitch; the blade faces are those lying on the blade boundary; the rest are split by
radius against the hub and shroud curves. The manifest contains the rotation axis and sector angle, plus starting
snippets for OpenFOAM cyclic patches and an SU2 `MARKER_PERIODIC` line.

**The snippets are starting points, not validated solver input.** In particular, check the sign and angle units of the
SU2 periodic rotation against your SU2 version, and confirm patch types and `transform` settings against your
OpenFOAM release.

## Pre-flight report

`--check` (and every build) prints derived quantities and warnings: blade count and angular pitch, mid-span
pitch, axial and true chord, pitch/chord and solidity, span, hub/tip ratio, aspect ratio, domain extent in chords,
tip clearance as absolute, % span and % chord, a first-cell size hint for the gap, and tip speed from `rpm`.
Errors (sections not spanning the annulus, tip clearance with `step` endwalls, sections not ordered hub to tip)
stop the build.

## Example cases

Three illustrative cases for a roughly 100-150 N class micro-turbojet are in `examples/`, regenerated by
`python examples/generate_examples.py`. They exercise the tool's features. They are **illustrative reference
geometries**, not aerodynamically optimised or validated designs, and no performance claim is made for them.
Values below are from the pre-flight report.

| | NGV | Turbine rotor | Axial compressor rotor |
|---|---|---|---|
| Directory | `micro_turbojet_ngv` | `micro_turbojet_turbine_rotor` | `micro_turbojet_axial_compressor_rotor` |
| Blades | 17 | 25 | 12 |
| Mid-span radius | 34.5 mm | 35.2 mm | 21.0 mm |
| Span | 9.0 mm | 10.9 mm | 16.9 mm |
| Hub/tip ratio | 0.77 | 0.73 | 0.44 |
| Axial chord | 12.1 mm | 8.1 mm | 10.1 mm |
| Solidity (mid-span) | 1.46 | 1.03 | 1.68 |
| Endwalls | parametric, constant annulus | parametric, diverging, cosine blend | `.fpd` curves, rising hub |
| Tip clearance | none (stator) | 0.30 mm LE to 0.25 mm TE (2.7 % span) | 1.2 % span (0.20 mm) |
| Speed / tip speed | n/a | 110 krpm, 469 m/s | 90 krpm, 283 m/s |
| Features shown | free-vortex outlet swirl | tapered clearance, diverging annulus | `.fpd` endwalls, %-span clearance |

Section data are generated by the in-repo generator with the stated metal angles and thickness ratios; the NGV
uses a free-vortex swirl distribution. Replace them with your own sections for real design work.

## Meshing and solver notes

* **Tip gap**: put at least 8-10 cells across the gap; the pre-flight hint gives a first-cell size (clearance / 10).
* **STL tolerance**: tessellation deflection defaults to 1/2000 of the span (`--stl-tolerance` overrides, in model
  units). A deflection chosen for millimetre-scale parts is far too coarse at this size.
* **ANSYS**: import the STEP (or the per-patch STEP files) and create Named Selections using the patch names or the
  centroids in the manifest; define `periodic_a` / `periodic_b` as a rotational periodic pair about Z.
* **OpenFOAM**: use the multi-region STL with snappyHexMesh, then the cyclic snippet from the manifest.
* **SU2**: mesh externally, keep the patch names as markers, then use the periodic line from the manifest after
  checking it as noted above.

None of these workflows has been run end to end with a solver as part of this repository's tests; the tests
verify the geometry, topology and export files (below).

## Limitations

* **Axial rows only.** Radial/centrifugal impellers and diffusers need a meridional channel definition and are not
  supported. This matters for micro-turbojets, whose compressors are very often centrifugal.
* Tip clearance needs `parametric` or `fpd` endwalls. No hub clearance, fillets or shrouded blades.
* Section hand: pressure side on +theta (mirror for the opposite hand).
* Build times of minutes per blade row; the blade cut dominates.
* Section definitions are used as supplied, with no aerodynamic checks beyond the pre-flight report.

## Testing

```bash
pytest -m "not slow"     # 33 fast tests, a few seconds
pytest -m slow           # 7 tests that build the turbine rotor domain: topology, volume, patches, tip gap, export (about 1 min)
```

The slow tests check that the domain is a single valid solid, that its volume equals the sector volume minus the
blade, that the sector volume agrees with the annulus geometry, that all eight patches are found, that the
periodic faces are a rotation of one pitch, and that the tip gap lies within the requested taper.

## What I contributed

MicroTurboGeo is derived from ProtoBlade (Tesseract Software Limited, Apache-2.0),
which supplies the core idea and code for building a periodic fluid sector from blade sections. The files and changes below
are mine; the complete, file-level record is in [`NOTICE`](NOTICE) and [`docs/CHANGES.md`](docs/CHANGES.md).

| Area | What I added | Where |
|---|---|---|
| Tip clearance | Gap modelled in the geometry: constant or LE-to-TE tapered, absolute or % span | `blade.py`, `cad.py` |
| Endwalls | Parametric hub and shroud (linear or cosine blend) | `stage.py`, `geom.py` |
| Pre-flight | `--check` report: solidity, pitch/chord, span, hub/tip ratio, clearance as % span and chord, first-cell hint, tip speed | `preflight.py` |
| CFD export | Named boundary patches, multi-region STL, per-patch STEP, JSON manifest with periodic-rotation data | `export.py` |
| Sections and examples | Parametric section generator, `.fpd` writer, three micro-turbojet example cases | `sections.py`, `examples/` |
| Robustness | Tolerances that scale with domain size, vectorised duplicate rejection, `arctan2`, validation of the machine axis, periodic surface extended past both endwalls | `geom.py`, `blade.py` |
| Tests | Section, config, pre-flight and full-domain tests | `tests/test_{sections,config,preflight,domain}.py` |

Modified from ProtoBlade: `blade`, `cad`, `geom`, `machine`, `stage`, `cli`, `__main__` and `tests/test_geom.py`.

## Validation

**Status: geometry and export are verified; solver workflows are not validated.**

| Level | What it shows | Status |
|---|---|---|
| Fast tests | Config parsing, section generation, geometry helpers, pre-flight numbers | 33 passed (2026-10-15) |
| Slow tests | Single valid solid; volume equals sector volume minus blade; all eight patches found; periodic faces are one pitch apart; tip gap within the requested taper | 7 passed (2026-10-15) |
| Example builds | All three shipped cases build and export | Passed, 28 s to 123 s |
| Export check | Turbine-rotor case: eight patches present, equal periodic face areas (2.3445e-4 m^2 each), sector angle 14.4 deg for 25 blades | Passed |
| Meshing | A mesh generated from the STEP or STL in ANSYS, OpenFOAM or SU2 with a quality report | Planned |
| Solver run | A periodic single-passage run with the OpenFOAM and SU2 snippets, checking mass-flow conservation across the periodic pair | Planned |
| Geometry vs reference | Compare a built passage with a published or hand-measured blade row | Planned |

The aerodynamic quality of the example sections is not assessed, and the OpenFOAM and SU2 snippets are unvalidated
starting points.

## Attribution and license

Licensed under the Apache License, Version 2.0 (see `LICENSE`). Changes relative to the upstream project are listed
in `docs/CHANGES.md`; copyright and attribution are recorded in `NOTICE`.

MicroTurboGeo is a modified derivative of ProtoBlade (© 2024 Tesseract Software Limited, Apache License 2.0).
