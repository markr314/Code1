"""Static geometry checks that do not need a solver run.

The expensive failure mode in an FDTD model of an all-metal structure is an
air volume that no port can reach: it rings as a lossless cavity and the
simulation never meets its energy end-criterion.  Rather than infer that from
a stalled run, rebuild the air regions analytically on the mesh and flood-fill
them from the port planes.
"""
from __future__ import annotations

from typing import List, Tuple

import numpy as np

from .params import BoifotParams
from . import geometry


def air_boxes(p: BoifotParams) -> List[Tuple[float, ...]]:
    """Every air volume, as (x0, x1, y0, y1, z0, z1); mirrors geometry.build."""
    ys = p.a_o / 2.0
    boxes = [
        (-p.a_o / 2, p.a_o / 2, -ys, ys, p.z_min, p.z_trans0),  # common arm
        (-p.a_o / 2, p.a_o / 2, -p.b_o / 2, p.b_o / 2, p.z_trans1, p.z_max),  # main
    ]
    for i, h in enumerate(p.step_heights):
        za = p.z_trans0 + i * p.step_len
        boxes.append((-p.a_o / 2, p.a_o / 2, -h / 2, h / 2, za, za + p.step_len))
    for sgn in (+1, -1) if p.side_arms_enabled else ():
        x_in, x_out = sgn * p.a_o / 2.0, sgn * p.x_side_end
        boxes.append(
            (min(x_in, x_out), max(x_in, x_out), -ys, ys, p.z_ap0, p.z_ap1)
        )
    return boxes


def occupancy(p: BoifotParams, lines=None):
    """Boolean air mask on the mesh cell centres, plus the cell-centre axes."""
    if lines is None:
        lines = geometry.mesh_lines(p)
    cx, cy, cz = (0.5 * (np.array(lines[d][:-1]) + np.array(lines[d][1:]))
                  for d in "xyz")
    air = np.zeros((cx.size, cy.size, cz.size), bool)
    for x0, x1, y0, y1, z0, z1 in air_boxes(p):
        ix = (cx > x0) & (cx < x1)
        iy = (cy > y0) & (cy < y1)
        iz = (cz > z0) & (cz < z1)
        air |= ix[:, None, None] & iy[None, :, None] & iz[None, None, :]

    # metal back in: septum (mitred wedge) and pins
    if p.dt_sep > 0:
        sy = np.abs(cy) < p.dt_sep / 2
        wedge = (cz[None, :] >= p.z_tip + p.miter_slope * np.maximum(
            np.abs(cx) - p.tip_w, 0.0)[:, None]) & (
            cz[None, :] <= p.z_sep_end
        )
        inx = np.abs(cx) <= p.a_o / 2
        air &= ~(inx[:, None, None] & sy[None, :, None] & wedge[:, None, :])
    for sgn in (+1, -1):
        for zc in p.pin_z:
            r = p.d_pin / 2
            d2 = (cx[:, None] - sgn * p.pin_x) ** 2 + (cz[None, :] - zc) ** 2
            air &= ~(d2 <= r * r)[:, None, :]
    return air, (cx, cy, cz)


def port_seeds(p: BoifotParams, axes) -> List[Tuple[int, int, int]]:
    """One seed cell inside each of the five modal ports' guides."""
    cx, cy, cz = axes
    ix = lambda v: int(np.argmin(np.abs(cx - v)))
    iy = lambda v: int(np.argmin(np.abs(cy - v)))
    iz = lambda v: int(np.argmin(np.abs(cz - v)))
    half_ap = 0.5 * (p.z_ap0 + p.z_ap1)
    return [
        (ix(0.0), iy(p.a_o / 4), iz(p.z_min + 0.5)),  # common arm
        (ix(0.0), iy(0.0), iz(p.z_max - 0.5)),  # main arm
        (ix(p.x_side_end - 0.5), iy(p.a_o / 4), iz(half_ap)),  # side +x
        (ix(-p.x_side_end + 0.5), iy(p.a_o / 4), iz(half_ap)),  # side -x
    ]


def check(p: BoifotParams, verbose: bool = True) -> dict:
    """Flood-fill the air from the ports; report anything left over."""
    air, axes = occupancy(p)
    cx, cy, cz = axes
    seen = np.zeros_like(air)
    stack = []
    for s in port_seeds(p, axes):
        if air[s]:
            stack.append(s)
            seen[s] = True
    while stack:
        i, j, k = stack.pop()
        for di, dj, dk in ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0),
                           (0, 0, 1), (0, 0, -1)):
            a, b, c = i + di, j + dj, k + dk
            if (0 <= a < air.shape[0] and 0 <= b < air.shape[1]
                    and 0 <= c < air.shape[2] and air[a, b, c]
                    and not seen[a, b, c]):
                seen[a, b, c] = True
                stack.append((a, b, c))

    orphan = air & ~seen
    n_air, n_orph = int(air.sum()), int(orphan.sum())
    res = {"air_cells": n_air, "orphan_cells": n_orph, "orphans": []}
    if n_orph:
        idx = np.argwhere(orphan)
        for i, j, k in idx[:: max(1, len(idx) // 8)][:8]:
            res["orphans"].append((float(cx[i]), float(cy[j]), float(cz[k])))
    if verbose:
        print(f"air cells {n_air}, unreachable from any port: {n_orph}")
        for o in res["orphans"]:
            print("  orphan near (x,y,z) = %.3f, %.3f, %.3f mm" % o)
    return res
