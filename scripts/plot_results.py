#!/usr/bin/env python3
"""Plot the Bøifot junction S-parameters as reflection / transmission /
isolation panels, in the style of ALMA Memo 425 Figures 6-8."""
import argparse
import os
import sys

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from boifot_omt import sweep as sw
from boifot_omt.params import BoifotParams, wr28_default

# Validated two-slot categorical palette (see references/palette.md).
THEME = {
    "light": dict(
        surface="#fcfcfb", ink="#0b0b0b", ink2="#52514e", grid="#dcdbd6",
        s1="#2a78d6", s2="#eb6834",
    ),
    "dark": dict(
        surface="#1a1a19", ink="#f0efec", ink2="#c3c2b7", grid="#3a3a38",
        s1="#3987e5", s2="#d95926",
    ),
}


def _panel(ax, x, series, t, ylabel, ylim):
    """One panel: at most two series, direct-labelled, with a legend."""
    ax.set_ylim(*ylim)
    # headroom on the right so the direct labels sit inside the axes
    span = x[-1] - x[0]
    ax.set_xlim(x[0] - 0.01 * span, x[-1] + 0.13 * span)
    for (label, y, color) in series:
        ax.plot(x, y, color=color, linewidth=2.0, label=label,
                solid_capstyle="round")
        # label the last point that is actually on-scale, so a curve running
        # off the bottom does not drag its label out of the figure
        ok = np.isfinite(y) & (y >= ylim[0]) & (y <= ylim[1])
        if ok.any():
            j = np.nonzero(ok)[0][-1]
            ax.plot([x[j]], [y[j]], marker="o", ms=4, color=color, zorder=5)
            ax.annotate(" " + label, (x[j], y[j]), color=t["ink"], fontsize=8,
                        va="center", ha="left", xytext=(5, 0),
                        textcoords="offset points", zorder=5)
    ax.set_ylabel(ylabel, color=t["ink2"], fontsize=9)
    ax.grid(True, color=t["grid"], linewidth=0.6, alpha=0.9)
    ax.set_axisbelow(True)
    ax.set_facecolor(t["surface"])
    for sp in ax.spines.values():
        sp.set_color(t["grid"])
    ax.tick_params(colors=t["ink2"], labelsize=8)
    leg = ax.legend(loc="lower left", fontsize=8, framealpha=0.0, ncol=2)
    for txt in leg.get_texts():
        txt.set_color(t["ink"])


def plot(freq, S, p: BoifotParams, out_png: str, mode: str = "light") -> None:
    t = THEME[mode]
    m = sw.metrics(S, freq, band=(p.f_lo, p.f_hi))
    x = freq / p.fc  # normalized frequency, as the memo plots it

    fig, axes = plt.subplots(
        3, 1, figsize=(7.2, 8.4), sharex=True, constrained_layout=True
    )
    fig.patch.set_facecolor(t["surface"])

    _panel(
        axes[0], x,
        [("main arm", m["S11_main_dB"], t["s1"]),
         ("side arm", m["S11_side_dB"], t["s2"])],
        t, "Reflection (dB)", (-50, 0),
    )
    _panel(
        axes[1], x,
        [("main arm", m["main_il_dB"], t["s1"]),
         ("side arms, combined", m["side_sum_dB"], t["s2"])],
        t, "Transmission (dB)", (-4.0, 0.3),
    )
    _panel(
        axes[2], x,
        [("side pol → main arm", m["iso_sidepol_to_main_dB"], t["s1"]),
         ("main pol → side arm", m["iso_mainpol_to_side_p_dB"], t["s2"])],
        t, "Isolation (dB)", (-180, 0),
    )

    axes[2].set_xlabel(
        f"Normalized frequency  f / f$_c$     (f$_c$ = {p.fc/1e9:.2f} GHz)",
        color=t["ink2"], fontsize=9,
    )
    axes[0].set_title(
        "Bøifot orthomode junction, WR-28 (openEMS FDTD)",
        color=t["ink"], fontsize=11, loc="left",
    )
    # mark the operating band edges
    for ax in axes:
        for fe in (p.f_lo / p.fc, p.f_hi / p.fc):
            ax.axvline(fe, color=t["ink2"], linewidth=0.8, linestyle=(0, (4, 3)),
                       alpha=0.55)

    fig.savefig(out_png, dpi=170, facecolor=t["surface"])
    plt.close(fig)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--results", default="results/wr28")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    d = np.load(os.path.join(args.results, "sparams.npz"))
    freq, S = d["freq"], d["S"]
    p = wr28_default()
    base = args.out or os.path.join(args.results, "sparams")
    for mode in ("light", "dark"):
        png = f"{base}_{mode}.png"
        plot(freq, S, p, png, mode)
        print(f"wrote {png}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
