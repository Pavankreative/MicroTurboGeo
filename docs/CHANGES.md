# Changes relative to ProtoBlade

MicroTurboGeo is derived from ProtoBlade (Copyright 2024 Tesseract Software Limited, Apache-2.0).
This file satisfies the Apache-2.0 requirement to state significant changes.

## New features
| Area | Change | File |
|---|---|---|
| Tip clearance | `tip_clearance`, `tip_clearance_te`, `tip_clearance_mode` (`absolute` / `percent_span`). The blade is trimmed against the shroud offset radially inwards, constant or tapering LE to TE. | `blade.py`, `cad.py` |
| Endwalls | `type = 'parametric'`: hub and shroud from end radii and axial limits, linear or cosine blend. | `stage.py`, `geom.py` |
| Blade count | `n_blade` or `pitch_deg`, validated for consistency and integer count. | `blade.py` |
| Pre-flight | `--check`: solidity, pitch/chord, span, hub/tip ratio, clearance as % span and chord, tip speed from `rpm`, section coverage of the annulus. | `preflight.py` |
| Export | Face classification into named patches; multi-region ASCII STL, per-patch STEP, JSON manifest with periodic-rotation data and OpenFOAM/SU2 snippets. | `export.py` |
| Sections | Parametric section generator (camber, thickness, blunt trailing edge) and `.fpd` writer. | `sections.py` |
| Examples | Three micro-turbojet cases and a script that regenerates them. | `examples/` |

## Fixes and robustness
- Midline duplicate tolerance was a fixed `1e-3` (1 mm in metres); it now scales with the domain length. Identical results on the original test fixture.
- The O(N^2) Python duplicate-rejection loop is vectorised with NumPy.
- `arctan(y/x)` replaced by `arctan2(y, x)`.
- The machine axis is validated to be +Z through the origin (the code always assumed this).
- The periodic surface is extended past both endwalls and past the blade ends (previously hub side only, and only when needed). Coincident end faces made OpenCASCADE booleans extremely slow.
- Midline spline point count is configurable (`DomainCreator.midline_points`, default 40); 200-point splines made booleans several times slower.
- File names inside a case file are relative to the case file, not the working directory.
- Deprecated `&` / `-` Workplane operators replaced by `intersect` / `cut`.

## Not implemented
- Radial / centrifugal blade rows (meridional channel definition).
- Tip clearance with CAD (`step`) endwalls.
- Hub clearance, fillets, shrouded blades.
