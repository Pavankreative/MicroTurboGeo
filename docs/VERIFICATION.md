# Verification — 21 November 2026

Python 3.12, Linux, CadQuery/OpenCASCADE installed.

- `python -m pytest tests -q`: **48 passed** (40 existing + 8 regressions); one host pytest/AnyIO deprecation warning.
- All three bundled micro-turbojet TOML cases passed `python -m microturbogeo <case.toml> --check`.
- The full test suite includes CAD solid validity, sector/blade volume, named boundaries, tip clearance, trimming, STEP/STL export and STEP reimport checks.

Fixed shifted arc-coordinate interpolation, duplicate-parameter handling and invalid-input errors. The fixture comparison uses rtol=1e-13 and atol=1e-15 instead of bitwise equality: the vectorised calculation differed by at most 2.23e-16 in the original fixture. No fixture values were regenerated to hide a geometric difference.

`requirements.txt` delegates to pyproject.toml. These checks establish software behaviour on bundled cases, not aerodynamic or experimental validation. 3D RANS runs remain pending. Windows/8 GB operation was not benchmarked in this Linux execution.
