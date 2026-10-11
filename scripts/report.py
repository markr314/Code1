#!/usr/bin/env python3
"""Summarise a finished sweep: power conservation, in-band figures, plots."""
import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from boifot_omt import sweep as sw
from boifot_omt.params import wr28_default
from boifot_omt.model import build_model


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--results", default="results/wr28")
    ap.add_argument("--no-plot", action="store_true")
    a = ap.parse_args()

    d = np.load(os.path.join(a.results, "sparams.npz"))
    freq, S = d["freq"], d["S"]
    p = wr28_default()
    band = (p.f_lo, p.f_hi)

    print(sw.unitarity_report(S, freq, band))
    print()

    m = sw.metrics(S, freq, band=band)
    print("in-band performance (port reference planes):")
    for k, v in m["summary"].items():
        if isinstance(v, tuple):
            print(f"  {k:46s}: {v[0]:+.2f} .. {v[1]:+.2f}")
        else:
            print(f"  {k:46s}: {v:+.2f}")

    Sd = sw.deembed(S, freq, p, build_model(p).deembed_mm)
    md = sw.metrics(Sd, freq, band=band)
    lo, hi = md["summary"]["side-arm relative phase, min/max (deg)"]
    print(f"\n  side-arm relative phase at the junction      : "
          f"{lo:+.2f} .. {hi:+.2f} deg")

    print("\n f/GHz   f/fc   main_IL  side_sum   RL_main  RL_side"
          "   iso_s>m  iso_m>s")
    step = max(1, len(freq) // 14)
    for i in range(0, len(freq), step):
        print("%6.2f %6.3f %8.2f %9.2f %9.1f %8.1f %9.1f %8.1f" % (
            freq[i] / 1e9, freq[i] / p.fc, m["main_il_dB"][i], m["side_sum_dB"][i],
            m["S11_main_dB"][i], m["S11_side_dB"][i],
            m["iso_sidepol_to_main_dB"][i], m["iso_mainpol_to_side_p_dB"][i]))

    if not a.no_plot:
        from plot_results import plot
        for mode in ("light", "dark"):
            png = os.path.join(a.results, f"sparams_{mode}.png")
            plot(freq, S, p, png, mode)
            print(f"\nwrote {png}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
