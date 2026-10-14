"""A set of functions and classes to represent a turbomachinery machine.

Modified for MicroTurboGeo (derived from ProtoBlade, (c) 2024 Tesseract Software Limited, Apache-2.0).
Changes: optional rotational speed (rpm) for tip-speed reporting; the machine axis is validated to be the
global Z axis, which the geometry code assumes.
"""
from __future__ import annotations
from atom.api import Atom, Int, List, Enum, Str, Typed, Property, Float, Tuple
from .stage import Stage

import os

import numpy as np
import tomli

UNITS = ['metres', 'millimetres']


class Machine(Atom):
    name = Str()
    n_blade = Int()
    rpm = Float(0.0)
    units = Enum(*UNITS)
    axis = Tuple(Tuple(float))
    stages = List(Stage)

    @classmethod
    def from_config_file(cls, fname: str) -> Machine:
        """Create instance of class from a toml file."""
        config = _read_toml(fname)

        stages_config = config.pop('stage')
        _check_axis(config['machine']['axis'])
        machine = cls(**config['machine'])
        stages = []
        for stage in stages_config:
            stages.append(Stage.from_config(stage['name'], stage['endwall'][0], stage['blade_section']))
        machine.stages = stages
        return machine


def _read_toml(fname: str) -> dict:
    """Read a toml configuration file."""
    with open(fname, mode="rb") as fp:
        config = tomli.load(fp)
    # file names inside a case file are relative to the case file, not to the current directory
    base = os.path.dirname(os.path.abspath(fname))
    for stage in config.get('stage', []):
        for section in stage.get('endwall', []) + stage.get('blade_section', []):
            for key, value in list(section.items()):
                if key.endswith('_fname') and value and not os.path.isabs(value):
                    section[key] = os.path.join(base, value)
    config['machine']['axis'] = tuple([tuple(x) for x in config['machine']['axis']])
    return config


def _check_axis(axis) -> None:
    """The geometry code assumes the machine axis is the global Z axis through the origin."""
    (p0, p1) = axis
    d = np.array(p1, dtype=float) - np.array(p0, dtype=float)
    if np.linalg.norm(p0) > 1e-12 or np.linalg.norm(np.cross(d, [0, 0, 1])) > 1e-12 or d[2] <= 0:
        raise ValueError('machine.axis must run along +Z through the origin: [[0,0,0],[0,0,1]]')
