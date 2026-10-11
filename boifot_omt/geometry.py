"""CSXCAD geometry construction for the Boifot junction.

The structure is built by filling the whole computational domain with metal
and then carving the waveguide volumes out of it with higher-priority "air"
primitives.  The septum and the compensation pins are then put back as metal
at a still higher priority.  That ordering makes the junction geometry a
handful of boxes instead of an assembly of walls.
"""
from __future__ import annotations

from typing import List
import numpy as np

from .params import BoifotParams

PRIO_FILL = 1
PRIO_AIR = 10
PRIO_METAL = 20


def _add_box(prop, start, stop, priority):
    prop.AddBox(priority=priority, start=start, stop=stop)


def build(CSX, p: BoifotParams) -> dict:
    """Add all primitives for the junction to `CSX`.

    Returns a dict of the landmark coordinates the caller needs for ports and
    mesh generation.
    """
    fill = CSX.AddMetal("shell")
    air = CSX.AddMaterial("air", epsilon=1.0, mue=1.0)
    septum = CSX.AddMetal("septum")

    xs, ys, zs = p.x_side_end, p.a_sq / 2.0, None
    z0, z1 = p.z_min, p.z_max

    # --- solid metal block covering the whole domain ----------------------
    _add_box(fill, [-xs, -ys, z0], [xs, ys, z1], PRIO_FILL)

    # --- square common guide, up to the start of the main-arm steps -------
    _add_box(
        air,
        [-p.a_sq / 2.0, -ys, z0],
        [p.a_sq / 2.0, ys, p.z_step0],
        PRIO_AIR,
    )

    # --- E-plane matching steps into the main arm -------------------------
    for i, w in enumerate(p.step_widths):
        za = p.z_step0 + i * p.step_len
        zb = za + p.step_len
        _add_box(air, [-w, -ys, za], [w, ys, zb], PRIO_AIR)

    # --- straight main arm ------------------------------------------------
    _add_box(
        air,
        [-p.b_main / 2.0, -ys, p.z_main_start],
        [p.b_main / 2.0, ys, z1],
        PRIO_AIR,
    )

    # --- side arms (and, where they meet the common guide, the apertures) --
    for sgn in (+1, -1):
        x_in = sgn * p.a_sq / 2.0
        x_out = sgn * xs
        _add_box(
            air,
            [min(x_in, x_out), -p.b_side / 2.0, p.z_ap0],
            [max(x_in, x_out), p.b_side / 2.0, p.z_ap1],
            PRIO_AIR,
        )

    # --- septum -----------------------------------------------------------
    if p.t_sep > 0:
        _add_box(
            septum,
            [-p.t_sep / 2.0, -ys, p.z_sep_start],
            [p.t_sep / 2.0, ys, p.z_sep_end],
            PRIO_METAL,
        )
    else:
        # zero-thickness PEC sheet on the x=0 mesh plane
        _add_box(septum, [0.0, -ys, p.z_sep_start], [0.0, ys, p.z_sep_end], PRIO_METAL)

    # --- compensation pins ------------------------------------------------
    pin_coords = []
    if p.pins_enabled and p.d_pin > 0:
        for sgn in (+1, -1):
            for frac in p.pin_z_frac:
                zc = p.z_ap0 + frac * p.a_side
                xc = sgn * p.a_sq / 2.0
                septum.AddCylinder(
                    priority=PRIO_METAL,
                    start=[xc, -p.b_side / 2.0, zc],
                    stop=[xc, p.b_side / 2.0, zc],
                    radius=p.d_pin / 2.0,
                )
                pin_coords.append((xc, zc))

    return {
        "x_lim": (-xs, xs),
        "y_lim": (-ys, ys),
        "z_lim": (z0, z1),
        "pins": pin_coords,
    }


def mesh_lines(p: BoifotParams) -> dict:
    """Critical mesh lines for the junction.

    Returns, per axis, a sorted line list.  Lines come in two classes:

    * hard - material boundaries, which must be present or the geometry is
      discretised wrongly;
    * soft - local refinement around the septum, the aperture edges and the
      pins, which may be dropped.

    Soft lines are only kept when they stay at least ``tol`` away from every
    line already kept.  Without that filter, independently generated
    refinement ladders can land a few microns apart, and since the FDTD
    timestep is set by the smallest cell anywhere in the mesh, a single
    stray sliver cell slows the whole simulation by orders of magnitude.
    """
    ys = p.a_sq / 2.0
    fine = p.mesh_res_fine
    tol = 0.5 * fine

    # ---------------- hard lines: material boundaries --------------------
    hx = {-p.x_side_end, p.x_side_end, -p.a_sq / 2.0, p.a_sq / 2.0,
          -p.b_main / 2.0, p.b_main / 2.0}
    for w in p.step_widths:
        hx.update({-w, w})
    if p.t_sep > 0:
        hx.update({-p.t_sep / 2.0, p.t_sep / 2.0})
    else:
        hx.add(0.0)
    if p.pins_enabled and p.d_pin > 0:
        for sgn in (+1, -1):
            hx.update({sgn * (p.a_sq / 2.0 - p.d_pin / 2.0),
                       sgn * (p.a_sq / 2.0 + p.d_pin / 2.0)})

    hy = {-ys, ys, -p.b_side / 2.0, p.b_side / 2.0}

    hz = {p.z_min, p.z_max, p.z_ap0, p.z_ap1, p.z_main_start,
          p.z_sep_start, p.z_sep_end}
    for i in range(p.n_steps + 1):
        hz.add(p.z_step0 + i * p.step_len)
    if p.pins_enabled and p.d_pin > 0:
        for frac in p.pin_z_frac:
            zc = p.z_ap0 + frac * p.a_side
            hz.update({zc - p.d_pin / 2.0, zc + p.d_pin / 2.0})

    # ---------------- soft lines: local refinement -----------------------
    sx, sy, sz = set(), set(), set()
    # across the septum
    for k in range(1, 5):
        sx.update({-p.t_sep / 2.0 - k * fine, p.t_sep / 2.0 + k * fine})
    # either side of each aperture plane
    for sgn in (+1, -1):
        for k in range(1, 4):
            sx.add(sgn * (p.a_sq / 2.0 - k * fine))
    # around the side-arm broad walls
    for k in range(1, 3):
        sy.update({-p.b_side / 2.0 - k * fine, p.b_side / 2.0 + k * fine,
                   -p.b_side / 2.0 + k * fine, p.b_side / 2.0 - k * fine})
    # aperture edges, septum leading edge, pin surroundings
    for k in range(1, 4):
        sz.update({p.z_ap0 + k * fine, p.z_ap1 - k * fine,
                   p.z_sep_start - k * fine, p.z_sep_start + k * fine})
    if p.pins_enabled and p.d_pin > 0:
        for frac in p.pin_z_frac:
            zc = p.z_ap0 + frac * p.a_side
            sz.add(zc)
            for k in range(1, 3):
                sz.update({zc - p.d_pin / 2.0 - k * fine,
                           zc + p.d_pin / 2.0 + k * fine})

    def assemble(hard, soft, lo, hi):
        keep = sorted({round(v, 9) for v in hard if lo - 1e-9 <= v <= hi + 1e-9})
        # a hard-line pair closer than tol is a geometry problem, not a mesh
        # one: report it rather than silently merging a feature away.
        for a_, b_ in zip(keep, keep[1:]):
            if b_ - a_ < 1e-6:
                continue
            if b_ - a_ < tol:
                raise ValueError(
                    f"material boundaries {a_:.4f} and {b_:.4f} mm are "
                    f"{(b_-a_)*1000:.1f} um apart, below the {tol*1000:.0f} um "
                    "mesh tolerance; adjust the dimensions or mesh_res_fine"
                )
        for v in sorted(soft):
            if not (lo - 1e-9 <= v <= hi + 1e-9):
                continue
            import bisect
            i = bisect.bisect_left(keep, v)
            near = min(
                [abs(v - keep[j]) for j in (i - 1, i) if 0 <= j < len(keep)] or [1e9]
            )
            if near >= tol:
                bisect.insort(keep, round(v, 9))
        return keep

    return {
        "x": assemble(hx, sx, -p.x_side_end, p.x_side_end),
        "y": assemble(hy, sy, -ys, ys),
        "z": assemble(hz, sz, p.z_min, p.z_max),
    }
