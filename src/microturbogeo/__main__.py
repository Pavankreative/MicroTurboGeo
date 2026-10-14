"""Entry point: ``python -m microturbogeo case.toml``.

Modified for MicroTurboGeo (derived from ProtoBlade, (c) 2024 Tesseract Software Limited, Apache-2.0).
Changes: pre-flight report, --check mode, patch classification and multi-format export.
"""
import os
import sys

from microturbogeo.machine import Machine
from microturbogeo.cli import create_parser


def main(fname, output_filename=None, check=False, formats=('step',), stl_tolerance=None):
    from microturbogeo import preflight
    if not output_filename:
        output_filename = os.fspath(fname).replace('.toml', '.step')

    machine = Machine.from_config_file(fname)
    ok = True
    for stage in machine.stages:
        for blade_def in stage.blades:
            report = preflight.run(machine, stage, blade_def)
            print(preflight.format_report(report, stage.name, blade_def.name))
            ok = ok and not report['errors']
    if check:
        return 0 if ok else 1
    if not ok:
        print('Pre-flight errors found - aborting. Fix the case file or run with --check for details.')
        return 1

    from microturbogeo.cad import DomainCreator
    from microturbogeo import export
    for stage in machine.stages:
        for blade_def in stage.blades:
            temp_out_fname = str(output_filename).replace('.step', f'-{stage.name}-{blade_def.name}.step')
            creator = DomainCreator(blade_def, stage.endwalls, machine.units, machine.axis)
            creator.create_domain()
            export.export_all(creator, machine, temp_out_fname, formats=formats, stl_tolerance=stl_tolerance)
    return 0


def run():
    args = create_parser().parse_args()
    sys.exit(main(args.filepath, args.output, args.check, tuple(args.formats), args.stl_tolerance))


if __name__ == "__main__":
    run()
