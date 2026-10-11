#!/usr/bin/env python3
"""Run the 4-port Boifot junction S-parameter sweep in openEMS."""
import argparse
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from boifot_omt.params import wr28_default, BoifotParams
from boifot_omt import sweep as sw
from boifot_omt.model import PORT_NAMES


def build_args() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", default="results/wr28", help="output directory")
    ap.add_argument("--nfreq", type=int, default=141)
    ap.add_argument("--threads", type=int, default=0, help="0 = all cores")
    ap.add_argument("--verbose", type=int, default=0)
    ap.add_argument(
        "--excite",
        default="all",
        help="comma-separated modal port indices to excite, or 'all'",
    )
    ap.add_argument("--mesh-res", type=float, default=None)
    ap.add_argument("--mesh-res-fine", type=float, default=None)
    ap.add_argument("--timesteps", type=int, default=None)
    ap.add_argument("--t-sep", type=float, default=None)
    ap.add_argument("--l-sep", type=float, default=None)
    ap.add_argument("--no-pins", action="store_true")
    ap.add_argument("--d-pin", type=float, default=None)
    ap.add_argument("--pin-z", default=None, help="e.g. 0.18,0.82")
    ap.add_argument("--n-steps", type=int, default=None)
    ap.add_argument("--step-len", type=float, default=None)
    return ap


def params_from_args(a) -> BoifotParams:
    p = wr28_default()
    for attr, val in [
        ("mesh_res", a.mesh_res),
        ("mesh_res_fine", a.mesh_res_fine),
        ("nr_timesteps", a.timesteps),
        ("t_sep", a.t_sep),
        ("l_sep", a.l_sep),
        ("d_pin", a.d_pin),
        ("n_steps", a.n_steps),
        ("step_len", a.step_len),
    ]:
        if val is not None:
            setattr(p, attr, val)
    if a.no_pins:
        p.pins_enabled = False
    if a.pin_z:
        p.pin_z_frac = tuple(float(v) for v in a.pin_z.split(","))
    p.__post_init__()
    return p


def main() -> int:
    a = build_args().parse_args()
    p = params_from_args(a)
    print(p.report())
    for w in p.warnings():
        print(f"WARNING: {w}")
    print()

    freq = np.linspace(p.f_lo, p.f_hi, a.nfreq)
    exc = (
        list(range(len(PORT_NAMES)))
        if a.excite == "all"
        else [int(v) for v in a.excite.split(",")]
    )

    t0 = time.time()
    S = sw.run_sweep(
        p, freq, a.out, excitations=exc, threads=a.threads, verbose=a.verbose
    )
    print(f"\ntotal FDTD time: {time.time()-t0:.0f} s")

    sw.save(a.out, freq, S, p)
    sw.write_touchstone(os.path.join(a.out, "boifot_wr28.s5p"), freq, S)

    Sd = sw.deembed(S, freq, p, [8.0, 8.0, 22.5, 10.0, 10.0])
    m = sw.metrics(S, freq, band=(p.f_lo, p.f_hi))
    md = sw.metrics(Sd, freq, band=(p.f_lo, p.f_hi))

    print("\n=== in-band performance (at the port reference planes) ===")
    for k, v in m["summary"].items():
        print(f"  {k:46s}: {v}")
    print("\n  side-arm relative phase, de-embedded to the junction:")
    print(f"    {md['summary']['side-arm relative phase, min/max (deg)']}")

    np.savez(
        os.path.join(a.out, "metrics.npz"),
        **{k: v for k, v in m.items() if k != "summary"},
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
