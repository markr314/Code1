#!/usr/bin/env python3
"""Validate the openEMS port/mesh/PML setup on a straight WR-28 guide.

A uniform PEC-walled guide has an analytic answer: |S21| = 1, |S11| = 0, and a
wave impedance Z0/sqrt(1-(fc/f)^2).  Any appreciable deviation means the port
planes, the PML, or the mesh are wrong, and would be mistaken for junction
behaviour later.
"""
import os
import sys
import tempfile

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from CSXCAD import ContinuousStructure
from openEMS import openEMS

from boifot_omt.waveguide import WR_SIZES, f_cutoff_te, wave_impedance_te

MM = 1e-3


def main() -> int:
    wr = WR_SIZES["WR-28"]
    a, b = wr.a, wr.b
    length = 30.0
    res = 0.3
    f_lo, f_hi = 26.5e9, 40.0e9

    FDTD = openEMS(NrTS=30000, EndCriteria=1e-5)
    FDTD.SetGaussExcite(0.5 * (f_lo + f_hi), 0.5 * (f_hi - f_lo))
    FDTD.SetBoundaryCond(["PEC", "PEC", "PEC", "PEC", "PML_8", "PML_8"])

    CSX = ContinuousStructure()
    FDTD.SetCSX(CSX)
    mesh = CSX.GetGrid()
    mesh.SetDeltaUnit(MM)
    mesh.AddLine("x", [0, a])
    mesh.AddLine("y", [0, b])
    mesh.AddLine("z", [0, length])

    d_exc, d_meas = 8 * res + 3 * res, 8 * res + 8 * res
    z = [d_exc, d_meas, length - d_exc, length - d_meas]
    mesh.AddLine("z", z)
    for d in "xyz":
        mesh.SmoothMeshLines(d, res, ratio=1.4)

    p1 = FDTD.AddRectWaveGuidePort(
        0, [0, 0, z[0]], [a, b, z[1]], "z", a * MM, b * MM, "TE10", 1
    )
    p2 = FDTD.AddRectWaveGuidePort(
        1, [a, b, z[2]], [0, 0, z[3]], "z", a * MM, b * MM, "TE10", 0
    )

    sim = os.path.join(tempfile.gettempdir(), "wr28_validate")
    cwd = os.getcwd()
    FDTD.Run(sim, cleanup=True, verbose=0, numThreads=0)
    os.chdir(cwd)

    freq = np.linspace(f_lo, f_hi, 141)
    p1.CalcPort(sim, freq)
    p2.CalcPort(sim, freq)

    s11 = p1.uf_ref / p1.uf_inc
    s21 = p2.uf_ref / p1.uf_inc
    ZL = p1.uf_tot / p1.if_tot
    fc = f_cutoff_te(a, b, 1, 0)
    ZL_a = np.array([wave_impedance_te(f, fc) for f in freq])

    # analytic transmission phase across the separation of the two
    # measurement planes
    d = (length - d_meas) - d_meas
    k0 = 2 * np.pi * freq / 299792458.0
    kc = 2 * np.pi * fc / 299792458.0
    bet = np.sqrt(k0**2 - kc**2)
    phase_a = -np.degrees(bet * d * MM) % 360
    phase_n = np.degrees(np.angle(s21)) % 360
    dphase = (phase_n - phase_a + 180) % 360 - 180

    print(f"cutoff          : {fc/1e9:.3f} GHz")
    print(f"|S11| worst     : {20*np.log10(np.abs(s11)).max():+.2f} dB")
    print(f"|S21| worst     : {20*np.log10(np.abs(s21)).min():+.4f} dB")
    print(f"ZL err worst    : {np.abs(ZL/ZL_a - 1).max()*100:.2f} %")
    print(f"S21 phase err   : {np.abs(dphase).max():.2f} deg over {d:.1f} mm")

    ok = True
    checks = [
        ("|S11| < -35 dB", 20 * np.log10(np.abs(s11)).max() < -35),
        ("|S21| > -0.05 dB", 20 * np.log10(np.abs(s21)).min() > -0.05),
        ("ZL within 2%", np.abs(ZL / ZL_a - 1).max() < 0.02),
        ("S21 phase within 3 deg", np.abs(dphase).max() < 3.0),
    ]
    for name, passed in checks:
        print(f"  [{'PASS' if passed else 'FAIL'}] {name}")
        ok &= passed
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
