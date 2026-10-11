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
    ap.add_argument("--mesh-res-junction", type=float, default=None)
    ap.add_argument("--timesteps", type=int, default=None)
    ap.add_argument("--dt-sep-norm", type=float, default=None,
                    help="septum thickness / a_o (memo: 0.02-0.05)")
    ap.add_argument("--pin-d-norm", type=float, default=None,
                    help="pin diameter / a_o (memo: 0.05-0.07)")
    ap.add_argument("--n-pins", type=int, default=None, help="pins per aperture")
    ap.add_argument("--no-pins", action="store_true")
    ap.add_argument("--miter-slope", type=float, default=None,
                    help="dz/dx of the septum tip edge; 1.0 = 45 deg")
    ap.add_argument("--sep-overhang", type=float, default=None)
    ap.add_argument("--trans-gap", type=float, default=None)
    ap.add_argument("--n-steps", type=int, default=None)
    ap.add_argument("--ap-len-norm", type=float, default=None)
    return ap


def params_from_args(a) -> BoifotParams:
    p = wr28_default()
    for attr, val in [
        ("mesh_res", a.mesh_res),
        ("mesh_res_fine", a.mesh_res_fine),
        ("mesh_res_junction", a.mesh_res_junction),
        ("nr_timesteps", a.timesteps),
        ("dt_sep_norm", a.dt_sep_norm),
        ("pin_d_norm", a.pin_d_norm),
        ("n_pins", a.n_pins),
        ("miter_slope", a.miter_slope),
        ("sep_overhang", a.sep_overhang),
        ("trans_gap", a.trans_gap),
        ("n_steps", a.n_steps),
        ("ap_len_norm", a.ap_len_norm),
    ]:
        if val is not None:
            setattr(p, attr, val)
    if a.no_pins:
        p.pins_enabled = False
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

    from boifot_omt.model import build_model
    Sd = sw.deembed(S, freq, p, build_model(p).deembed_mm)
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
