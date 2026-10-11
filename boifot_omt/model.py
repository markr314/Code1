"""openEMS FDTD model assembly for the 4-port Boifot junction.

Modal port numbering (the common port is dual-mode, so a 4-port junction has
a 5x5 modal S-matrix):

    0  common guide, polH  (E along x)  -> main arm
    1  common guide, polV  (E along y)  -> side arms
    2  main arm
    3  side arm at +x
    4  side arm at -x
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional
import os

import numpy as np
from CSXCAD import ContinuousStructure
from openEMS import openEMS

from .params import BoifotParams
from . import geometry

PORT_NAMES = ["common_H", "common_V", "main", "side_p", "side_m"]
PORT_LABELS = {
    "common_H": "common (polH, E||x)",
    "common_V": "common (polV, E||y)",
    "main": "main arm",
    "side_p": "side arm +x",
    "side_m": "side arm -x",
}
MM = 1e-3


@dataclass
class Model:
    FDTD: openEMS
    CSX: ContinuousStructure
    ports: List
    params: BoifotParams
    # distance (mm) from each port's measurement plane to the junction
    # reference plane, for phase de-embedding
    deembed_mm: List[float]


def build_model(p: BoifotParams, excite_port: Optional[int] = None) -> Model:
    """Assemble FDTD + geometry + the five modal ports.

    `excite_port` selects which modal port carries the excitation; pass None
    for a geometry-only build (mesh inspection, XML export).
    """
    f0 = 0.5 * (p.f_lo + p.f_hi)
    fc = 0.5 * (p.f_hi - p.f_lo)

    FDTD = openEMS(NrTS=int(p.nr_timesteps), EndCriteria=p.end_criteria)
    FDTD.SetGaussExcite(f0, fc)
    # x: side arms leave the domain; y: the guide walls themselves; z: common
    # port and main arm leave the domain.
    FDTD.SetBoundaryCond(
        ["PML_8", "PML_8", "PEC", "PEC", "PML_8", "PML_8"]
    )

    CSX = ContinuousStructure()
    FDTD.SetCSX(CSX)
    mesh = CSX.GetGrid()
    mesh.SetDeltaUnit(MM)

    geometry.build(CSX, p)

    lines = geometry.mesh_lines(p)
    for d in "xyz":
        mesh.AddLine(d, lines[d])

    pml = p.pml_cells * p.mesh_res
    ys = p.a_sq / 2.0

    # ---- port planes -----------------------------------------------------
    # excitation plane sits just inside the PML; the measurement plane is
    # five cells further into the structure.
    d_exc = pml + 3 * p.mesh_res
    d_meas = pml + 8 * p.mesh_res

    z_c_exc = p.z_min + d_exc
    z_c_meas = p.z_min + d_meas
    z_m_exc = p.z_max - d_exc
    z_m_meas = p.z_max - d_meas
    x_s_exc = p.x_side_end - d_exc
    x_s_meas = p.x_side_end - d_meas

    for d in "xz":
        pass
    mesh.AddLine("z", [z_c_exc, z_c_meas, z_m_exc, z_m_meas])
    mesh.AddLine("x", [x_s_exc, x_s_meas, -x_s_exc, -x_s_meas])

    mesh.SmoothMeshLines("x", p.mesh_res, ratio=1.4)
    mesh.SmoothMeshLines("y", p.mesh_res, ratio=1.4)
    mesh.SmoothMeshLines("z", p.mesh_res, ratio=1.4)

    exc = lambda n: (1 if excite_port == n else 0)
    ports = []

    # 0/1: common square guide, both polarisations, propagating +z
    a_m = p.a_sq * MM
    ports.append(
        FDTD.AddRectWaveGuidePort(
            0,
            [-p.a_sq / 2, -ys, z_c_exc],
            [p.a_sq / 2, ys, z_c_meas],
            "z",
            a_m,
            a_m,
            "TE01",  # E along x  -> polH
            exc(0),
        )
    )
    ports.append(
        FDTD.AddRectWaveGuidePort(
            1,
            [-p.a_sq / 2, -ys, z_c_exc],
            [p.a_sq / 2, ys, z_c_meas],
            "z",
            a_m,
            a_m,
            "TE10",  # E along y  -> polV
            exc(1),
        )
    )

    # 2: main arm, propagating -z into the structure
    ports.append(
        FDTD.AddRectWaveGuidePort(
            2,
            [-p.b_main / 2, -ys, z_m_exc],
            [p.b_main / 2, ys, z_m_meas],
            "z",
            p.b_main * MM,
            p.a_sq * MM,
            "TE01",  # E along x, cutoff set by the y dimension
            exc(2),
        )
    )

    # 3: side arm at +x, propagating -x into the structure
    ports.append(
        FDTD.AddRectWaveGuidePort(
            3,
            [x_s_exc, -p.b_side / 2, p.z_ap0],
            [x_s_meas, p.b_side / 2, p.z_ap1],
            "x",
            p.b_side * MM,
            p.a_side * MM,
            "TE01",  # E along y, cutoff set by the z dimension
            exc(3),
        )
    )

    # 4: side arm at -x, propagating +x into the structure
    ports.append(
        FDTD.AddRectWaveGuidePort(
            4,
            [-x_s_exc, -p.b_side / 2, p.z_ap0],
            [-x_s_meas, p.b_side / 2, p.z_ap1],
            "x",
            p.b_side * MM,
            p.a_side * MM,
            "TE01",
            exc(4),
        )
    )

    # distance from each measurement plane to the junction reference plane
    deembed = [
        abs(z_c_meas - p.z_ap0),
        abs(z_c_meas - p.z_ap0),
        abs(z_m_meas - p.z_ap0),
        abs(x_s_meas - p.a_sq / 2.0),
        abs(x_s_meas - p.a_sq / 2.0),
    ]

    return Model(FDTD=FDTD, CSX=CSX, ports=ports, params=p, deembed_mm=deembed)


def mesh_summary(m: Model) -> str:
    mesh = m.CSX.GetGrid()
    n = [mesh.GetQtyLines(d) for d in range(3)]
    total = n[0] * n[1] * n[2]
    lines = [
        f"mesh: {n[0]} x {n[1]} x {n[2]} = {total/1e6:.2f} M cells",
    ]
    for i, d in enumerate("xyz"):
        L = np.array(mesh.GetLines(i))
        dl = np.diff(L)
        lines.append(
            f"  {d}: {L[0]:+.3f} .. {L[-1]:+.3f} mm, "
            f"d{d} min {dl.min():.4f} max {dl.max():.4f} mm"
        )
    return "\n".join(lines)
