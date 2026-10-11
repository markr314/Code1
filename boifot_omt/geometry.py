"""CSXCAD geometry for the Bøifot junction (ALMA Memo 425 arrangement).

The domain is filled with metal and the guide volumes are carved out of it
with higher-priority "air" primitives; the septum and pins then go back in as
metal at a higher priority still.
"""
from __future__ import annotations

from typing import List
import bisect

import numpy as np

from .params import BoifotParams

PRIO_FILL = 1
PRIO_AIR = 10
PRIO_METAL = 20


def build(CSX, p: BoifotParams) -> dict:
    """Add every primitive for the junction to `CSX`."""
    fill = CSX.AddMetal("shell")
    air = CSX.AddMaterial("air", epsilon=1.0, mue=1.0)
    metal = CSX.AddMetal("septum_pins")

    ys = p.a_o / 2.0  # the common arm and the side arms are full height in y
    xs = p.x_side_end

    # --- solid metal block over the whole domain -------------------------
    fill.AddBox(priority=PRIO_FILL, start=[-xs, -ys, p.z_min], stop=[xs, ys, p.z_max])

    # --- square common arm, up to the start of the transformer ------------
    air.AddBox(
        priority=PRIO_AIR,
        start=[-p.a_o / 2, -ys, p.z_min],
        stop=[p.a_o / 2, ys, p.z_trans0],
    )

    # --- four-step Chebyshev height transformer ---------------------------
    za = p.z_trans0
    for h, L in zip(p.step_heights, p.step_lengths):
        air.AddBox(
            priority=PRIO_AIR,
            start=[-p.a_o / 2, -h / 2.0, za],
            stop=[p.a_o / 2, h / 2.0, za + L],
        )
        za += L

    # --- straight main arm ------------------------------------------------
    air.AddBox(
        priority=PRIO_AIR,
        start=[-p.a_o / 2, -p.b_o / 2, p.z_trans1],
        stop=[p.a_o / 2, p.b_o / 2, p.z_max],
    )

    # --- side arms; where they meet the common arm they form the apertures
    for sgn in (+1, -1) if p.side_arms_enabled else ():
        x_in, x_out = sgn * p.a_o / 2.0, sgn * xs
        air.AddBox(
            priority=PRIO_AIR,
            start=[min(x_in, x_out), -ys, p.z_ap0],
            stop=[max(x_in, x_out), ys, p.z_ap1],
        )

    # --- septum: thin plate in y=0 with a mitred leading edge -------------
    # Footprint in the (z, x) plane: a wedge opening from the tip apex at
    # (z=0, x=0) back to the side walls, then a full-width rectangle.  For
    # the side-arm polarisation this edge is the pair of back-to-back miters
    # that divides the power between the two side arms.
    if p.dt_sep > 0:
        def tip_poly(dz: float):
            """Septum footprint in the (z, x) plane, shifted dz downstream."""
            if p.tip_w > 0.0 or p.tip_len > 0.0:
                z_pts = [p.z_tip + dz, p.z_tip + dz, p.z_tip_end + dz,
                         p.z_sep_end, p.z_sep_end, p.z_tip_end + dz]
                x_pts = [-p.tip_w, p.tip_w, p.a_o / 2, p.a_o / 2,
                         -p.a_o / 2, -p.a_o / 2]
            else:
                z_pts = [p.z_tip + dz, p.z_tip_end + dz, p.z_sep_end,
                         p.z_sep_end, p.z_tip_end + dz]
                x_pts = [0.0, p.a_o / 2, p.a_o / 2, -p.a_o / 2, -p.a_o / 2]
            return np.array([z_pts, x_pts])

        layers = p.taper_layers
        if not layers:
            layers = [(p.dt_sep / 2.0, 0.0)]
        for half_t, dz in layers:
            metal.AddLinPoly(
                priority=PRIO_METAL,
                points=tip_poly(dz),
                norm_dir=1,
                elevation=-half_t,
                length=2.0 * half_t,
            )

    # --- compensation pins -------------------------------------------------
    pins = []
    for sgn in (+1, -1):
        for zc in p.pin_z:
            xc = sgn * p.pin_x
            if p.pin_span == "half":
                spans = [(p.dt_sep / 2.0, ys), (-ys, -p.dt_sep / 2.0)]
            else:
                spans = [(-ys, ys)]
            for y0, y1 in spans:
                metal.AddCylinder(
                    priority=PRIO_METAL,
                    start=[xc, y0, zc],
                    stop=[xc, y1, zc],
                    radius=p.d_pin / 2.0,
                )
            pins.append((xc, zc))

    return {"x_lim": (-xs, xs), "y_lim": (-ys, ys), "z_lim": (p.z_min, p.z_max),
            "pins": pins}


def mesh_lines(p: BoifotParams) -> dict:
    """Critical mesh lines.

    Hard lines are material boundaries and must survive; soft lines are local
    refinement and are dropped when they fall within half a fine cell of a
    line already kept.  The FDTD timestep is set by the smallest cell anywhere
    in the mesh, so a stray sliver cell is expensive.
    """
    ys = p.a_o / 2.0
    fine = p.mesh_res_fine
    tol = 0.75 * fine
    res_j = p.mesh_res_junction

    # ---------------- hard: material boundaries --------------------------
    hx = {-p.x_side_end, p.x_side_end, -p.a_o / 2.0, p.a_o / 2.0}
    if p.tip_w > 0.0:
        hx.update({-p.tip_w, p.tip_w})
    if p.pins_enabled and p.n_pins > 0:
        for sgn in (+1, -1):
            hx.update({sgn * (p.pin_x - p.d_pin / 2.0), sgn * (p.pin_x + p.d_pin / 2.0)})

    hy = {-ys, ys, -p.b_o / 2.0, p.b_o / 2.0}
    for h in p.step_heights:
        hy.update({-h / 2.0, h / 2.0})
    if p.dt_sep > 0:
        hy.update({-p.dt_sep / 2.0, p.dt_sep / 2.0})
        for half_t, _ in p.taper_layers:
            hy.update({-half_t, half_t})
    else:
        hy.add(0.0)

    hz = {p.z_min, p.z_max, p.z_tip, p.z_tip_end, p.z_ap0, p.z_ap1, p.z_trans0,
          p.z_trans1, p.z_sep_end}
    za = p.z_trans0
    hz.add(za)
    for L in p.step_lengths:
        za += L
        hz.add(za)
    for zc in p.pin_z:
        hz.update({zc - p.d_pin / 2.0, zc + p.d_pin / 2.0})
    for _, dz in p.taper_layers:
        hz.update({p.z_tip + dz, p.z_tip_end + dz})

    # ---------------- soft: local refinement ------------------------------
    sx, sy, sz = set(), set(), set()

    # the junction region carries the staircased mitred edge, so resolve both
    # x and z there more finely than the straight guides
    n_j = max(int(round(p.a_o / res_j)), 2)
    for i in range(n_j + 1):
        sx.add(-p.a_o / 2.0 + i * p.a_o / n_j)
    n_jz = max(int(round((p.z_ap1 - p.z_ap0) / res_j)), 2)
    for i in range(n_jz + 1):
        sz.add(p.z_ap0 + i * (p.z_ap1 - p.z_ap0) / n_jz)

    # across the septum
    for k in range(1, 5):
        sy.update({-p.dt_sep / 2.0 - k * fine, p.dt_sep / 2.0 + k * fine})
    # around the main-arm step corners
    for h in p.step_heights + [p.b_o]:
        for k in (1, 2):
            sy.update({-h / 2.0 - k * fine, h / 2.0 + k * fine})
    # around the septum tip and the aperture edges in z
    for k in range(1, 4):
        sz.update({p.z_tip + k * fine, p.z_tip - k * fine,
                   p.z_ap1 - k * fine, p.z_ap1 + k * fine})
    # around the pins
    for zc in p.pin_z:
        sz.add(zc)
        for k in (1, 2):
            sz.update({zc - p.d_pin / 2.0 - k * fine, zc + p.d_pin / 2.0 + k * fine})
    if p.pins_enabled and p.n_pins > 0:
        # A round post is staircased, and if its surface is under-resolved the
        # steps form a weakly coupled sliver resonator that rings for the whole
        # run.  Resolve the post itself with several cells across, and keep that
        # resolution a little way out from its surface.
        step = p.d_pin / 4.0
        reach = 1.5 * p.d_pin
        n = int(np.ceil(reach / step))
        for sgn in (+1, -1):
            for k in range(-n, n + 1):
                sx.add(sgn * p.pin_x + k * step)
            for zc in p.pin_z:
                for k in range(-n, n + 1):
                    sz.add(zc + k * step)

    def assemble(hard, soft, lo, hi):
        keep = sorted({round(v, 9) for v in hard if lo - 1e-9 <= v <= hi + 1e-9})
        for a_, b_ in zip(keep, keep[1:]):
            if 1e-6 < b_ - a_ < tol:
                raise ValueError(
                    f"material boundaries {a_:.4f} and {b_:.4f} mm are "
                    f"{(b_-a_)*1000:.1f} um apart, below the {tol*1000:.0f} um mesh "
                    "tolerance; adjust the dimensions or mesh_res_fine"
                )
        for v in sorted(soft):
            if not (lo - 1e-9 <= v <= hi + 1e-9):
                continue
            i = bisect.bisect_left(keep, v)
            near = min([abs(v - keep[j]) for j in (i - 1, i) if 0 <= j < len(keep)]
                       or [1e9])
            if near >= tol:
                bisect.insort(keep, round(v, 9))
        return keep

    return {
        "x": assemble(hx, sx, -p.x_side_end, p.x_side_end),
        "y": assemble(hy, sy, -ys, ys),
        "z": assemble(hz, sz, p.z_min, p.z_max),
    }
