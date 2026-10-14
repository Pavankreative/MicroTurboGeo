"""Command line interface for MicroTurboGeo."""
import argparse


def create_parser():
    parser = argparse.ArgumentParser(
        prog='microturbogeo',
        description='Create a CFD-ready periodic fluid domain from blade section definitions.')
    parser.add_argument('filepath', help='Location of the TOML case file.')
    parser.add_argument('-o', '--output', default=None,
                        help='Output STEP file name (stage and blade names are appended). Default: <case>.step')
    parser.add_argument('--check', action='store_true',
                        help='Run the geometry pre-flight report only (no CAD is built).')
    parser.add_argument('--formats', nargs='+', default=['step'],
                        choices=['step', 'stl', 'patch-step'],
                        help='Export formats. step: full domain. stl: one multi-region ASCII STL with a named '
                             'solid per boundary patch (OpenFOAM snappyHexMesh, SU2 via mesher). patch-step: one '
                             'STEP per boundary patch (ANSYS named selections).')
    parser.add_argument('--stl-tolerance', type=float, default=None,
                        help='STL linear deflection in model units. Default: 1/2000 of the annulus span.')
    return parser
