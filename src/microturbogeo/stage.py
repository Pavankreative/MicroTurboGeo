"""A set of functions and classes to represent a turbomachinery stage and its endwalls.

Modified for MicroTurboGeo (derived from ProtoBlade, (c) 2024 Tesseract Software Limited, Apache-2.0).
Changes: new 'parametric' endwall type (conical / smooth-blend hub and shroud from a handful of radii).
"""
from __future__ import annotations

import numpy as np
from atom.api import Atom, Int,List,Enum,Str,Typed,Property,Float,Tuple,Instance
from atom.atom import Atom
from atom.enum import Enum
from atom.scalars import Str
from atom.typed import Typed
from microturbogeo import geom
from .blade import Blade

ENDWALL_TYPES = ['fpd', 'step', 'parametric']

class Endwalls(Atom):
    """Hold objects required to define endwalls."""

    hub = Typed(np.ndarray)
    shroud = Typed(np.ndarray)
    type= Enum(*ENDWALL_TYPES)

    step_fname = Str()
    hub_fname = Str()
    shroud_fname = Str()

    @classmethod
    def from_config(cls,config:dict) -> Endwalls:
        if config['type'] not in ENDWALL_TYPES:
            raise ValueError('Invalid endwall type')

        endwall = Endwalls()
        endwall.type = config['type']

        if config['type'] == 'fpd':
            endwall.hub = geom.load_curves_from_fpd(config['hub_fname'])
            endwall.shroud = geom.load_curves_from_fpd(config['shroud_fname'])
        elif config['type'] == 'parametric':
            missing = [k for k in ('hub_radius', 'shroud_radius', 'z_inlet', 'z_outlet') if k not in config]
            if missing:
                raise ValueError(f'parametric endwalls need: {", ".join(missing)}')
            n_points = int(config.get('n_points', 21))
            blend = config.get('blend', 'cosine')
            z_in, z_out = float(config['z_inlet']), float(config['z_outlet'])
            endwall.hub = geom.make_parametric_meridional(*config['hub_radius'], z_in, z_out, n_points, blend)
            endwall.shroud = geom.make_parametric_meridional(*config['shroud_radius'], z_in, z_out, n_points, blend)
            if np.any(geom.radial_of(endwall.shroud) <= geom.radial_of(endwall.hub)):
                raise ValueError('shroud radius must exceed hub radius everywhere')

        endwall.step_fname = config.get('step_fname','')

        return endwall


class Stage(Atom):
    name = Str()
    blades = List(Blade)
    endwalls = Instance(Endwalls)

    @classmethod
    def from_config(
            cls,
            name,
            endwall_config,
            blade_config,
        ):
        stage = Stage()
        stage.name = name
        stage.endwalls = Endwalls.from_config(endwall_config)
        stage.blades = [Blade.from_config(config) for config in blade_config]
        return stage
