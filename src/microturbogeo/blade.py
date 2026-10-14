"""A set of functions and classes to represent a turbomachinery blade row.

Modified for MicroTurboGeo (derived from ProtoBlade, (c) 2024 Tesseract Software Limited, Apache-2.0).
Changes: pitch/blade-count control (n_blade or pitch_deg), tip clearance definition (absolute or percent of
span, optionally tapering from LE to TE), section extent helpers.
"""
from __future__ import annotations
from numpy.typing import NDArray
import numpy as np
from microturbogeo import geom
from atom.api import Atom, Int, Enum,Str,Typed,Property, Tuple,Float

TIP_CLEARANCE_MODES = ['absolute', 'percent_span']

class Blade(Atom):
    #TODO rename BladeSection?
    """Represents a turbomachinery blade."""
    name = Str()
    n_blade = Int()
    n_sections = Property(Int)
    ps_sections = Typed(np.ndarray)
    ss_sections = Typed(np.ndarray)
    interface_location = Float(0.0)

    #: tip clearance at the blade leading edge. Units depend on tip_clearance_mode
    tip_clearance = Float(0.0)
    #: optional tip clearance at the trailing edge (same mode). Negative means 'same as leading edge'
    tip_clearance_te = Float(-1.0)
    tip_clearance_mode = Enum(*TIP_CLEARANCE_MODES)


    pitch_angle_rad = Property()

    @classmethod
    def from_config(cls,config):
        blade = Blade()
        blade.name = config['name']
        blade.n_blade = _resolve_blade_count(config)

        for surf,fname in  [('ps_sections',config['ps_section_fname']),('ss_sections',config['ss_section_fname'])]:
            blade._load_sections(fname,surf)

        if 'interface_location' in config.keys():
            blade.interface_location = config['interface_location']

        blade.tip_clearance_mode = config.get('tip_clearance_mode', 'absolute')
        blade.tip_clearance = float(config.get('tip_clearance', 0.0))
        blade.tip_clearance_te = float(config.get('tip_clearance_te', -1.0))
        if blade.tip_clearance < 0:
            raise ValueError('tip_clearance must not be negative')
        return blade

    def _load_sections(self,fname:str,surf:str):
        try:
            setattr(self,surf,geom.load_curves_from_fpd(fname))
        except FileNotFoundError:
            raise

    def _get_n_sections(self) -> int:
        return len(self.ps_sections)


    def _get_pitch_angle_rad(self) ->float:
        return 2.0 * np.pi / self.n_blade

    def radial_extent(self) -> tuple:
        """Return (r_min, r_max) over all supplied PS and SS sections."""
        r = np.concatenate([geom.radial_of(self.ps_sections).ravel(), geom.radial_of(self.ss_sections).ravel()])
        return float(r.min()), float(r.max())

    def axial_extent(self) -> tuple:
        """Return (z_le, z_te): axial extent of the sections."""
        z = np.concatenate([self.ps_sections['z'].ravel(), self.ss_sections['z'].ravel()])
        return float(z.min()), float(z.max())

    def clearance_at(self, span: float) -> tuple:
        """Return absolute tip clearance (at LE, at TE) in model units for a given annulus span.

        Args:
            span: hub-to-shroud radial height (model units); only used in 'percent_span' mode.
        """
        scale = span / 100.0 if self.tip_clearance_mode == 'percent_span' else 1.0
        c_le = self.tip_clearance * scale
        c_te = c_le if self.tip_clearance_te < 0 else self.tip_clearance_te * scale
        return c_le, c_te


def _resolve_blade_count(config: dict) -> int:
    """Blade count from either n_blade or pitch_deg (angular pitch). Both may be given if consistent."""
    n = config.get('n_blade')
    pitch = config.get('pitch_deg')
    if n is None and pitch is None:
        raise ValueError("A blade row needs 'n_blade' or 'pitch_deg'")
    if pitch is not None:
        if pitch <= 0:
            raise ValueError('pitch_deg must be positive')
        n_from_pitch = 360.0 / pitch
        if abs(n_from_pitch - round(n_from_pitch)) > 1e-3:
            raise ValueError(f'pitch_deg={pitch} does not divide 360 degrees into an integer number of blades')
        if n is not None and int(n) != round(n_from_pitch):
            raise ValueError('n_blade and pitch_deg are inconsistent')
        n = round(n_from_pitch)
    if int(n) < 2:
        raise ValueError('A periodic sector needs at least 2 blades')
    return int(n)


#TODO: this should use the radial extrusion function
def create_sections_from_2D_profile(ps:NDArray,ss:NDArray,N_sections:int,r_extents:tuple,use_r_theta:bool=True,del_x:float=0,n_resample:int=0):
    """
    Create a set of sections from the R-Theta Z plane to the cartesian plane across a radius span.

    Args:
        ps: array of Z and r-Theta points for pressure surface
        ss: array of Z and r-Theta points for pressure surface
        N_sections: Number of sections to create across the radial span
        r_extents: radial extents of the span, the first entry of the tuple must be lower than the second
        use_r_theta: True will extrude the geometry with constant R-Theta co-ordinates across the span, False will
                        extrude the geometry with constant Theta values
        del_x: an offset to apply in the z/x-axis
        n_resample: set to a non-zero value to resample an array.

    Returns:
        ps_sections : NxM array of new pressure surface sections
        ss_sections : NxM array of new suction surface sections
    """
    x_ps = ps['x']
    rt_ps = ps['y']

    x_ss = ss['x']
    rt_ss = ss['y']
    r_min, r_max = r_extents

    if n_resample > 1:
        x_ps,rt_ps,_ = geom.reinterpolate_curve(x_ps,rt_ps,geom.calculate_curve_length(x_ps,rt_ps),base=1.5,N_new=n_resample)
        x_ss,rt_ss,_ = geom.reinterpolate_curve(x_ss,rt_ss,geom.calculate_curve_length(x_ss,rt_ss),base=1.5,N_new=n_resample)

        N_ps = x_ps.shape[0]
        N_ss = x_ss.shape[0]
    else:
        N_ps = len(x_ps)
        N_ss = len(x_ss)

    t_ps = rt_ps / r_min
    t_ss = rt_ss / r_min

    r = np.linspace(r_min, r_max, N_sections)


    ps_section=[]
    ss_section =[]
    for i in range(N_sections):
        ps_out = np.empty(shape=(N_ps,), dtype=geom.cartesian_type)
        ss_out = np.empty(shape=(N_ss,), dtype=geom.cartesian_type)
        if use_r_theta:
            theta_ps = rt_ps / r[i]
            theta_ss = rt_ss / r[i]
        else:
            theta_ps = t_ps
            theta_ss = t_ss

        ps_out['x'] = r[i] * np.cos(theta_ps)
        ps_out['y'] = r[i] * np.sin(theta_ps)
        ps_out['z'] = x_ps +del_x

        ss_out['x'] = r[i] * np.cos(theta_ss)
        ss_out['y'] = r[i] * np.sin(theta_ss)
        ss_out['z'] = x_ss +del_x

        ps_section.append(ps_out)
        ss_section.append(ss_out)

    return ps_section,ss_section