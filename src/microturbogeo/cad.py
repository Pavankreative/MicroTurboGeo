"""Module with classes and functions to create CAD models from microturbogeo classes.

Modified for MicroTurboGeo (derived from ProtoBlade, (c) 2024 Tesseract Software Limited, Apache-2.0).
Changes: tip-clearance trimming of the blade against an offset shroud envelope (constant or LE->TE tapered),
tip-side extension of the periodic surface, parametric endwall support, non-deprecated boolean calls.
"""
import cadquery
from numpy.typing import NDArray
from typing import Tuple,List, Literal
import numpy as np
from  microturbogeo import  geom, stage
from microturbogeo.blade import Blade
def _convert_array_to_list(pts:NDArray)-> List[Tuple]:
    return [tuple(pt) for pt in pts]

class DomainCreator:
    """Class to create a Fluid Domain that can be used in CFD analysis."""

    blade_def : Blade
    #: points per spanwise midline of the periodic surface. Lower = faster booleans.
    midline_points = 40

    def __init__(self,
                 blade_def:Blade,
                 endwalls:stage.Endwalls,
                 units:str,
                 axis:tuple,
                 cq=cadquery
                 ):
        """Create the object from a Stage instance."""
        #TODO : probaly want this to be a stage rather than blade - actually maybe not?
        self.blade_def = blade_def
        self.endwalls = endwalls
        self.units = units
        self.axis = axis
        #TODO: add function to cut domain at a given location
        self._cq = cq # add cadquery as an object to allow for test mock to be easily added

    def extrude_blade(self):
        """Extrude/loft the blade sections to create the main blade."""
        N_sections = self.blade_def.ps_sections.shape[0]

        ps_edges = []
        ss_edges = []
        for i in range(N_sections):
            ss_pts = _convert_array_to_list(self.blade_def.ss_sections[i])
            ss_edges.append(self._cq.Edge.makeSpline([self._cq.Vector(p) for p in ss_pts]))

            ps_pts = _convert_array_to_list(self.blade_def.ps_sections[i])
            ps_edges.append(self._cq.Edge.makeSpline([self._cq.Vector(p) for p in ps_pts]))

        blade_ss = self._cq.Solid.makeLoft(
            [self._cq.Wire.assembleEdges([edge]) for edge in ss_edges]
        )

        blade_ps = self._cq.Solid.makeLoft(
            [self._cq.Wire.assembleEdges([edge]) for edge in ps_edges]
        )

        top = self._cq.Solid.makeLoft(
            [self._cq.Wire.assembleEdges([edge]) for edge in [ps_edges[0], ss_edges[0]]]
        )

        bottom = self._cq.Solid.makeLoft(
            [self._cq.Wire.assembleEdges([edge]) for edge in [ps_edges[-1], ss_edges[-1]]]
        )

        shell = self._cq.Shell.makeShell([blade_ss.faces(), bottom.faces(), blade_ps.faces(), top.faces()])
        solid = self._cq.Solid.makeSolid(shell)


        self.blade = solid

    def create_endwalls(self):
        """Create CAD objects for the endwalls."""
        if self.endwalls.type in ('fpd', 'parametric'):
            hub_pts = _convert_array_to_list(self.endwalls.hub)
            shroud_pts = _convert_array_to_list(self.endwalls.shroud)

            self.cad_endwalls =  self._cq.Workplane("XY").spline(hub_pts).polyline([hub_pts[-1], shroud_pts[-1]]).spline(
                shroud_pts[::-1]).polyline(
                [shroud_pts[0], hub_pts[0]]).close().revolve(360.0, self.axis[0], self.axis[1])
        else:
            self.cad_endwalls = self._cq.importers.importStep(self.endwalls.step_fname)

    def export(self,entity:str,fname_out:str)->None:
        """Export an entity from this class to a CAD output format.

        If the entity does not exist in the instance of this class that no export will occur

        Args:
            entity : name of entity to export
            fname_out : output file name with suffix to denote the desired output type


        """
        to_export = getattr(self,entity)
        if to_export:
            self._cq.exporters.export(to_export, str(fname_out))


    def create_periodic(self):
        """Create a CAD object to represent the periodic fluid domain."""
        z_min = self.cad_endwalls.objects[0].BoundingBox().zmin
        z_max = self.cad_endwalls.objects[0].BoundingBox().zmax

        edges = []

        endwall_min_r,endwall_max_r = find_radial_extent_of_axisymmetric_object(self.cad_endwalls)

        ps_sections = self._extend_sections(self.blade_def.ps_sections, endwall_min_r, endwall_max_r)
        ss_sections = self._extend_sections(self.blade_def.ss_sections, endwall_min_r, endwall_max_r)

        mid_points = geom.create_midlines(ps_sections,ss_sections,z_min,z_max,self.blade_def.pitch_angle_rad,
                                       n_out=self.midline_points)
        for i in range(len(mid_points)):
            pts = _convert_array_to_list(mid_points[i])
            edges.append(self._cq.Edge.makeSpline([self._cq.Vector(p) for p in pts]))

        per = self._cq.Solid.makeLoft(
            [self._cq.Wire.assembleEdges([edge]) for edge in edges]
        )

        self.per  = self._cq.Solid.revolve(per.Faces()[0], -np.rad2deg(self.blade_def.pitch_angle_rad), self.axis[0] , self.axis[1])


    @staticmethod
    def _extend_sections(sections, endwall_min_r, endwall_max_r):
        """Extend the periodic surface radially beyond the annulus *and* beyond the blade ends.

        The sector is later intersected with the endwalls, so any overshoot is trimmed away. Overshooting by
        10 % of the span keeps the sector's radial end caps away from the blade's end caps; coincident or
        near-coincident faces make the OpenCASCADE booleans extremely slow or fail outright.
        """
        margin = 0.1 * (endwall_max_r - endwall_min_r)
        r_first = np.min(geom.radial_of(sections[0]))
        r_last = np.max(geom.radial_of(sections[-1]))
        target_in = min(r_first, endwall_min_r) - margin
        target_out = max(r_last, endwall_max_r) + margin
        inner = geom.extrude_radially(sections[0], target_in - r_first)
        outer = geom.extrude_radially(sections[-1], target_out - r_last)
        return np.concatenate((np.reshape(inner, (1, inner.shape[0])), sections,
                               np.reshape(outer, (1, outer.shape[0]))))

    def annulus_span(self) -> float:
        """Hub-to-shroud radial height at the blade mid-chord (fpd / parametric endwalls only)."""
        z_le, z_te = self.blade_def.axial_extent()
        z_mid = 0.5 * (z_le + z_te)
        return float(geom.interp_radius(self.endwalls.shroud, z_mid) - geom.interp_radius(self.endwalls.hub, z_mid))

    def tip_gap_profile(self):
        """Return (z_le, z_te, c_le, c_te) for the requested tip clearance, or None if there is none."""
        if self.blade_def.tip_clearance <= 0:
            return None
        if self.endwalls.type == 'step':
            raise ValueError('Tip clearance needs fpd or parametric endwalls (the shroud curve is offset radially)')
        span = self.annulus_span()
        c_le, c_te = self.blade_def.clearance_at(span)
        if max(c_le, c_te) >= 0.5 * span:
            raise ValueError(f'tip clearance ({max(c_le, c_te):g}) is unreasonably large for a span of {span:g}')
        z_le, z_te = self.blade_def.axial_extent()
        return z_le, z_te, c_le, c_te

    def apply_tip_clearance(self):
        """Trim the blade against the shroud offset radially inwards by the tip clearance.

        The fluid domain keeps the resulting gap between the blade tip and the casing, so the tip-leakage
        passage is part of the CFD geometry. The blade sections must reach at least the shroud radius.
        The clearance may vary linearly from the leading edge (tip_clearance) to the trailing edge
        (tip_clearance_te); it is held constant upstream and downstream of the blade.
        """
        prof = self.tip_gap_profile()
        self.tip_clearance_applied = (0.0, 0.0)
        if prof is None:
            return
        z_le, z_te, c_le, c_te = prof
        shroud = self.endwalls.shroud
        frac = np.clip((shroud['z'] - z_le) / (z_te - z_le), 0.0, 1.0)
        gap = c_le + (c_te - c_le) * frac
        shroud_off = geom.radial_offset_curve(shroud, -gap)

        r_low = 0.5 * float(np.min(geom.radial_of(self.endwalls.hub)))
        pts = _convert_array_to_list(shroud_off)
        low_end = (r_low, 0.0, float(shroud_off['z'][-1]))
        low_start = (r_low, 0.0, float(shroud_off['z'][0]))
        self.tip_envelope = self._cq.Workplane("XY").spline(pts).polyline(
            [pts[-1], low_end, low_start, pts[0]]).close().revolve(360.0, self.axis[0], self.axis[1])

        self.blade_untrimmed = self.blade
        blade_wp = self._cq.Workplane("XY").add(self.blade)
        trimmed = blade_wp.intersect(self.tip_envelope)
        self.blade = trimmed.val()
        self.tip_clearance_applied = (float(c_le), float(c_te))

    def create_domain(self):
        """

        Create the full CFD domain.

        This will most likely be quite time consuming as it creates each aspect and then performs Boolean operations
        to create a single solid. T

        """
        self.extrude_blade()
        self.create_endwalls()
        self.apply_tip_clearance()
        self.create_periodic()

        blade_wp = self._cq.Workplane("XY").add(self.blade)
        per_wp = self._cq.Workplane("XY").add(self.per)
        per_and_endwalls = per_wp.intersect(self.cad_endwalls)
        self.domain = per_and_endwalls.cut(blade_wp)


def find_radial_extent_of_axisymmetric_object(input:cadquery.Workplane)->(float,float):
    """
    Find the radial extent of an axisymmetric object by taking a slice at Y=0 (therefore X=R).

    Args:
        input: cad query Workplane object that contains the object of interest

    Returns:
        (rmin,rmax) , the minimum and maximum radial values respectively

    """
    result2 = \
        input \
            .split(
            cadquery.Workplane()
            .add(cadquery.Face
            .makePlane(
                5, 5,
                cadquery.Vector(0, 0, 0.0),
                cadquery.Vector(1, 0, 0)))).solids(">>X")

    result3 = \
        result2 \
            .split(
            cadquery.Workplane()
            .add(cadquery.Face
            .makePlane(
                5, 5,
                cadquery.Vector(0, 0, 0.0),
                cadquery.Vector(0, 1, 0)))).solids("<<Y")

    result4 = result3.faces(">Y")

    bb = result4.objects[0].BoundingBox()
    return (bb.xmin,bb.xmax)