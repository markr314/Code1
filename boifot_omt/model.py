"""openEMS FDTD model assembly for the 4-port Bøifot junction.

Modal port numbering.  The common port is dual-mode, so a 4-port junction
has a 5x5 modal S-matrix:

    0  common arm, main-arm polarisation (E || y, normal to the septum)
    1  common arm, side-arm polarisation (E || x, in the septum plane)
    2  main arm
    3  side arm at +x
    4  side arm at -x
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional

import numpy as np
from CSXCAD import ContinuousStructure
from openEMS import openEMS

from .params import BoifotParams
from . import geometry

PORT_NAMES = ["common_M", "common_S", "main", "side_p", "side_m"]
PORT_LABELS = {
    "common_M": "common arm, main-arm pol (E||y)",
    "common_S": "common arm, side-arm pol (E||x)",
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
    deembed_mm: List[float]  # measurement plane -> junction reference plane


def build_model(p: BoifotParams, excite_port: Optional[int] = None) -> Model:
    f0 = 0.5 * (p.f_lo + p.f_hi)
    fc = 0.5 * (p.f_hi - p.f_lo)

    FDTD = openEMS(NrTS=int(p.nr_timesteps), EndCriteria=p.end_criteria)
    FDTD.SetGaussExcite(f0, fc)
    # x: the side arms leave the domain; y: the common-arm walls themselves;
    # z: common arm and main arm leave the domain.
    FDTD.SetBoundaryCond(["PML_8", "PML_8", "PEC", "PEC", "PML_8", "PML_8"])

    CSX = ContinuousStructure()
    FDTD.SetCSX(CSX)
    mesh = CSX.GetGrid()
    mesh.SetDeltaUnit(MM)

    geometry.build(CSX, p)
    lines = geometry.mesh_lines(p)
    for d in "xyz":
        mesh.AddLine(d, lines[d])

    pml = p.pml_cells * p.mesh_res
    d_exc = pml + 3 * p.mesh_res
    d_meas = pml + 8 * p.mesh_res

    z_c_exc, z_c_meas = p.z_min + d_exc, p.z_min + d_meas
    z_m_exc, z_m_meas = p.z_max - d_exc, p.z_max - d_meas
    x_s_exc, x_s_meas = p.x_side_end - d_exc, p.x_side_end - d_meas

    mesh.AddLine("z", [z_c_exc, z_c_meas, z_m_exc, z_m_meas])
    mesh.AddLine("x", [x_s_exc, x_s_meas, -x_s_exc, -x_s_meas])
    for d in "xyz":
        mesh.SmoothMeshLines(d, p.mesh_res, ratio=1.4)

    exc = lambda n: (1 if excite_port == n else 0)
    ys = p.a_o / 2.0
    a_m = p.a_o * MM
    ports = []

    # 0/1: common arm, both polarisations, propagating +z
    ports.append(
        FDTD.AddRectWaveGuidePort(
            0, [-p.a_o / 2, -ys, z_c_exc], [p.a_o / 2, ys, z_c_meas],
            "z", a_m, a_m, "TE10", exc(0),   # E || y -> main-arm pol
        )
    )
    ports.append(
        FDTD.AddRectWaveGuidePort(
            1, [-p.a_o / 2, -ys, z_c_exc], [p.a_o / 2, ys, z_c_meas],
            "z", a_m, a_m, "TE01", exc(1),   # E || x -> side-arm pol
        )
    )

    # 2: main arm, propagating -z into the structure
    ports.append(
        FDTD.AddRectWaveGuidePort(
            2, [-p.a_o / 2, -p.b_o / 2, z_m_exc], [p.a_o / 2, p.b_o / 2, z_m_meas],
            "z", a_m, p.b_o * MM, "TE10", exc(2),  # E || y
        )
    )

    # 3/4: side arms, E || z, broad dimension along y
    ports.append(
        FDTD.AddRectWaveGuidePort(
            3, [x_s_exc, -ys, p.z_ap0], [x_s_meas, ys, p.z_ap1],
            "x", a_m, p.ap_len * MM, "TE10", exc(3),
        )
    )
    ports.append(
        FDTD.AddRectWaveGuidePort(
            4, [-x_s_exc, -ys, p.z_ap0], [-x_s_meas, ys, p.z_ap1],
            "x", a_m, p.ap_len * MM, "TE10", exc(4),
        )
    )

    deembed = [
        abs(z_c_meas - p.z_tip),
        abs(z_c_meas - p.z_tip),
        abs(z_m_meas - p.z_tip),
        abs(x_s_meas - p.a_o / 2.0),
        abs(x_s_meas - p.a_o / 2.0),
    ]
    return Model(FDTD=FDTD, CSX=CSX, ports=ports, params=p, deembed_mm=deembed)


def mesh_summary(m: Model) -> str:
    mesh = m.CSX.GetGrid()
    n = [mesh.GetQtyLines(d) for d in range(3)]
    out = [f"mesh: {n[0]} x {n[1]} x {n[2]} = {n[0]*n[1]*n[2]/1e6:.2f} M cells"]
    dmin = []
    for i, d in enumerate("xyz"):
        L = np.array(mesh.GetLines(i))
        dl = np.diff(L)
        dmin.append(dl.min())
        out.append(
            f"  {d}: {L[0]:+.3f} .. {L[-1]:+.3f} mm, "
            f"d{d} min {dl.min():.4f} max {dl.max():.4f} mm"
        )
    c0 = 299792458.0
    dt = 1.0 / (c0 * np.sqrt(sum((1.0 / (v * MM)) ** 2 for v in dmin)))
    out.append(f"  Courant timestep ~ {dt*1e15:.1f} fs")
    return "\n".join(out)
