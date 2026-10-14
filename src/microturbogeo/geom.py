"""A set of functions to undertake geometrical manipulations.

Modified for MicroTurboGeo (derived from ProtoBlade, (c) 2024 Tesseract Software Limited, Apache-2.0).
Changes: arctan2 for angle recovery, scale-invariant midline tolerance, vectorised Voronoi duplicate
rejection, tip-side radial extension helpers, radial offset / interpolation helpers for tip clearance.
"""
import math
import pathlib
import numpy as np
from numpy.typing import NDArray
from scipy.interpolate import interp1d
from typing import Tuple
from enum import Enum

cartesian_type = np.dtype([("x", np.double), ("y", np.double), ("z", np.double)])
polar_type = np.dtype([("r", np.double), ("theta", np.double), ("z", np.double)])


def load_curves_from_fpd(fname: str or pathlib.Path) -> NDArray:
    """Load a series of curves from a formatted point data (fpd) file."""
    with open(fname, 'r') as f:
        n_pts, n_curve = list(map(int, f.readline().split()))
        pts = np.array([tuple(map(float, line.split())) for line in f.readlines()], dtype=cartesian_type)

    if n_curve > 1:
        pts = np.reshape(pts, (n_curve, n_pts))

    return pts


def convert_to_polar(pts: NDArray) -> NDArray:
    """
    Create arrays for radius , Theta and Z from cartesian array.

    Args:
        pts: array of cartesian points dtype cartesian

     Returns:
        An array which can be accessed through ['z'],['r'] and['theta' for Z, Radius and Theta respectively

    """
    out = np.empty(shape=(pts['z'].shape[0],), dtype=polar_type)

    out['z'] = pts['z']
    out['r'] = (pts['x'] ** 2.0 + pts['y'] ** 2) ** 0.5
    out['theta'] = np.arctan2(pts['y'], pts['x'])

    return out


def calculate_curve_length(x: NDArray, y: NDArray) -> NDArray:
    """Create an array of curve length based on the two input arrays x,y."""
    s = np.zeros(x.shape[0])
    for i in range(1, x.shape[0]):
        dx = x[i] - x[i - 1]
        dy = y[i] - y[i - 1]
        s[i] = s[i - 1] + np.sqrt(dx * dx + dy * dy)
    return s


def reinterpolate_curve(x: NDArray, y: NDArray, s: NDArray, base: float = 1.5, N_new: int = 0) -> Tuple[
    NDArray, NDArray, NDArray]:
    """
    Reinterpolate the arrays x,y and s based on a power law to a given base.

    Args:
        x: array for axial co-ordinate
        y: array for y co-ordinate
        s: array for curve length
        base: base for power law should be >1. The higher the number the closer the clustering
        N_new : Half the number points for the newly reintpolated array. Set to 0 to keep the same as the input
    Returns
        A tuple of three arrays newly interpolated arrays : x_new,y_new,s_new

    """
    N = x.shape[0]
    fx = interp1d(s, x)
    fy = interp1d(s, y)

    del_s = max(s) - min(s)

    start = 0.0 ** (1.0 / base)
    mid = 0.5 ** (1.0 / base)
    if N_new == 0:
        N_new = N
    s_new = np.zeros(N_new)
    for i in range(1, N_new + 1):
        s_new[i - 1] = (start + (i - 1) * (mid - start) / (N_new - 1)) ** base

    s_2 = np.flipud(s_new)

    s_3 = 1 - s_2
    s_4 = np.delete(s_3, 0)

    s_new = np.concatenate((s_new, s_4)) * del_s

    x_new = np.zeros(len(s_new))
    y_new = np.zeros(len(s_new))
    for i in range(len(s_new)):
        x_new[i] = fx(s_new[i])
        y_new[i] = fy(s_new[i])

    return x_new, y_new, s_new


def calculate_curvature() -> NDArray:
    """Caluclate curvature of x,y."""
    return


def create_midlines(ps_sections, ss_sections, z_min, z_max, pitch_angle_rad,n_resample: int = 0, n_out: int = 200):
    """
    Create curves to represent the midline of a pressure and section surfaces.

    Args:
        ps_sections: array of pressure surface sections N sections with M points
        ss_sections: array of suction surface sections N sections with M points
        z_min: minimum z value for the midline
        z_max: maximum z value for the midline
        pitch_angle_rad : pitch angle between blades in radians
        n_resample: number of points in the reinterpolated midpoint, set to 0 to skip reinterpolation
        n_out: number of points on each output midline when n_resample is 0. Fewer points give a lighter
            B-spline and much faster CAD booleans; 30-60 is plenty for a smooth periodic surface.

    Returns:
        midpoints : array of midpoint curves : Nsections with M points

    """
    mid_points = []
    N_sections = ps_sections.shape[0]
    from scipy.signal import resample
    for i in range(N_sections):
        section_ = np.concatenate((ps_sections[i], ss_sections[i][::-1]))
        section_polar = convert_to_polar(section_)
        rad = np.mean(section_polar['r'])

        z = section_polar['z']
        rt = section_polar['r'] * section_polar['theta']

        if n_resample > 0:
            z = resample(z, n_resample)
            rt = resample(rt, n_resample)
        else:
            n_resample = n_out
        # tolerance scales with the domain length (reproduces the original 1e-3 for a 0.05 domain)
        tol = 2e-2 * (z_max - z_min)

        mid_points.append(create_midline(n_resample, pitch_angle_rad, rad, rt, tol, z, z_max, z_min))
    return mid_points


def create_midline(Nout, pitch_angle_rad, rad, rt, tol, z, z_max, z_min):
    """
    Find the midline between two sections based on Voronoi's algorithm.

    Args:
        Nout: Number of points in the ouput array
        pitch_angle_rad: pitch angle in radians (angle between adjacent blades)
        rad: radius of current section
        rt: array of points in the r-Theta plane
        tol: tolerance used when finding duplicate points
        z: array of points in the z plane
        z_max: maximum z value for final mid line curve
        z_min: minimum z value for final mid line curve

    Returns:
        midpoints_cart: array of midline points in cartesian co-ordinate system

    """
    from scipy.spatial import Voronoi
    from scipy.interpolate import interp1d

    blade_pitch = rad * pitch_angle_rad
    points_single = np.zeros([len(z), 2])
    points = np.zeros([len(z) * 2 - 1, 2])
    for i in range(0, len(z)):
        points[i][0] = rt[i]
        points[i][1] = z[i]
        points_single[i][0] = rt[i]
        points_single[i][1] = z[i]
    for i in range(0, len(z)):
        points[len(z) - 1 + i][0] = rt[i] + blade_pitch
        points[len(z) - 1 + i][1] = z[i]

    # look at voronoi
    vor_single = Voronoi(points_single)
    vor = Voronoi(points)

    v = vor.vertices
    vs = vor_single.vertices
    # A vertex is a duplicate if it coincides (in z and r-theta) with a vertex of the single-blade diagram,
    # either as-is or shifted by one blade pitch.
    close_z = np.abs(v[:, 1][:, None] - vs[:, 1][None, :]) < tol
    close_rt = np.abs(v[:, 0][:, None] - vs[:, 0][None, :]) < tol
    close_rt_shift = np.abs(v[:, 0][:, None] - (vs[:, 0] + blade_pitch)[None, :]) < tol
    duplicate = np.any(close_z & (close_rt | close_rt_shift), axis=1)
    keep = (~duplicate) & (v[:, 1] > z_min) & (v[:, 1] < z_max) & (v[:, 0] > np.min(points[:, 0]) + blade_pitch * 0.25)
    mid_rt = list(v[keep, 0])
    mid_z = list(v[keep, 1])
    N_mid = len(mid_rt)
    midpoints = np.zeros([N_mid, 2])
    for i in range(0, N_mid):
        midpoints[i][0] = mid_rt[i]
        midpoints[i][1] = mid_z[i]
    midpoints = midpoints[np.argsort(midpoints[:, 1])]

    # interpolate midpoints
    z_int = np.linspace(z_min, z_max, Nout)
    f2 = interp1d(midpoints[:, 1], midpoints[:, 0], fill_value='extrapolate')
    midpoints_out = np.zeros([Nout, 2])
    midpoints_cart = np.zeros([Nout, 3])
    for i in range(Nout):
        midpoints_out[i][1] = z_int[i]
        midpoints_out[i][0] = f2(z_int[i])

        midpoints_cart[i][2] = z_int[i]
        t = f2(z_int[i]) / rad
        midpoints_cart[i][0] = rad * np.cos(t)
        midpoints_cart[i][1] = rad * np.sin(t)

    return midpoints_cart


def extrude_radially(input: NDArray, delta_r: float) -> NDArray:
    """Take an array of X,Y,Z points and extrude radially inwards or outwards.

    Args
        input : input array to be extruded
        delta_r : amount to be extruded in metres. Positive for increase in radius, negative for decrease

    Returns
        output : output array at new radial location

    Raises
        ValueError : if new radial position is less than zero
    """

    new_r = np.hypot(input['x'], input['y']) + delta_r

    if np.where(new_r < 0)[0].size > 0:
        raise ValueError('New radius has gone below zero')

    t = np.arctan2(input['y'], input['x'])

    output = np.zeros(input.size, dtype=cartesian_type)

    output['x'] = new_r * np.cos(t)
    output['y'] = new_r * np.sin(t)
    output['z'] = input['z']

    return output



def radial_of(curve: NDArray) -> NDArray:
    """Return radius of every point of a cartesian-typed array."""
    return np.hypot(curve['x'], curve['y'])


def interp_radius(curve: NDArray, z) -> NDArray:
    """Radius of an axisymmetric meridional curve (hub or shroud) at axial position(s) z."""
    order = np.argsort(curve['z'])
    return np.interp(z, curve['z'][order], radial_of(curve)[order])


def radial_offset_curve(curve: NDArray, delta) -> NDArray:
    """Shift a meridional curve radially by delta (scalar or per-point array). Negative moves towards the axis."""
    delta = np.broadcast_to(np.asarray(delta, dtype=float), curve.shape)
    new_r = radial_of(curve) + delta
    if np.any(new_r <= 0):
        raise ValueError('Radial offset moves curve onto or through the axis')
    t = np.arctan2(curve['y'], curve['x'])
    out = np.zeros(curve.shape, dtype=cartesian_type)
    out['x'] = new_r * np.cos(t)
    out['y'] = new_r * np.sin(t)
    out['z'] = curve['z']
    return out


def make_parametric_meridional(r_in: float, r_out: float, z_in: float, z_out: float, n_points: int = 21,
                               blend: str = 'cosine') -> NDArray:
    """Build a hub or shroud meridional curve (y=0 plane) running between two radii.

    Args:
        r_in: radius at z_in
        r_out: radius at z_out
        z_in: axial position of the domain inlet
        z_out: axial position of the domain outlet
        n_points: number of points on the curve (>=3)
        blend: 'linear' (conical) or 'cosine' (smooth S-shaped transition with zero slope at both ends)
    """
    if n_points < 3:
        raise ValueError('n_points must be at least 3')
    if z_out <= z_in:
        raise ValueError('z_outlet must be greater than z_inlet')
    t = np.linspace(0.0, 1.0, n_points)
    if blend == 'cosine':
        f = 0.5 * (1.0 - np.cos(np.pi * t))
    elif blend == 'linear':
        f = t
    else:
        raise ValueError("blend must be 'linear' or 'cosine'")
    out = np.zeros(n_points, dtype=cartesian_type)
    out['x'] = r_in + (r_out - r_in) * f
    out['z'] = z_in + (z_out - z_in) * t
    out['z'][0], out['z'][-1] = z_in, z_out          # exact end planes for the inlet/outlet faces
    out['x'][0], out['x'][-1] = r_in, r_out
    return out
