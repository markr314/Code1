"""Drive the FDTD runs and assemble the modal S-matrix."""
from __future__ import annotations

import os
from typing import Iterable, Optional, Sequence

import numpy as np

from .params import BoifotParams
from .model import build_model, mesh_summary, PORT_NAMES, PORT_LABELS
from .waveguide import beta as wg_beta, f_cutoff_te

NP = len(PORT_NAMES)


def port_cutoffs(p: BoifotParams) -> np.ndarray:
    """Cutoff frequency (Hz) of each modal port's operating mode."""
    return np.array(
        [
            p.fc_common,  # common polH
            p.fc_common,  # common polV
            p.fc_main,  # main arm
            p.fc_side,  # side arm +x
            p.fc_side,  # side arm -x
        ]
    )


def run_sweep(
    p: BoifotParams,
    freq: np.ndarray,
    out_dir: str,
    excitations: Optional[Sequence[int]] = None,
    threads: int = 0,
    verbose: int = 0,
    log=print,
) -> np.ndarray:
    """Run one FDTD simulation per excited port; return S of shape (nf, 5, 5).

    ``S[:, j, i]`` is the wave leaving port j for unit incidence on port i.
    Columns belonging to excitations that were not run are left as NaN.
    """
    if excitations is None:
        excitations = range(NP)
    os.makedirs(out_dir, exist_ok=True)
    S = np.full((len(freq), NP, NP), np.nan, dtype=complex)
    cwd = os.getcwd()
    try:
        for exc in excitations:
            m = build_model(p, excite_port=exc)
            if exc == list(excitations)[0]:
                log(mesh_summary(m))
            sim = os.path.join(os.path.abspath(out_dir), f"exc_{PORT_NAMES[exc]}")
            log(f"--- exciting port {exc} ({PORT_LABELS[PORT_NAMES[exc]]}) -> {sim}")
            m.FDTD.Run(sim, cleanup=True, verbose=verbose, numThreads=threads)
            os.chdir(cwd)
            for port in m.ports:
                port.CalcPort(sim, freq)
            inc = m.ports[exc].uf_inc
            for j, port in enumerate(m.ports):
                S[:, j, exc] = port.uf_ref / inc
    finally:
        os.chdir(cwd)
    return S


def deembed(
    S: np.ndarray, freq: np.ndarray, p: BoifotParams, dist_mm: Sequence[float]
) -> np.ndarray:
    """Shift every port's reference plane to the junction reference plane.

    Multiplies out the electrical length of the straight run between each
    port's measurement plane and the junction, using the analytic TE
    propagation constant of that port's mode.
    """
    fc = port_cutoffs(p)
    phase = np.zeros((len(freq), NP))
    for i in range(NP):
        b = np.array([wg_beta(f, fc[i]) for f in freq])
        phase[:, i] = b * dist_mm[i] * 1e-3
    out = np.empty_like(S)
    for i in range(NP):
        for j in range(NP):
            out[:, j, i] = S[:, j, i] * np.exp(1j * (phase[:, i] + phase[:, j]))
    return out


def db20(x) -> np.ndarray:
    return 20.0 * np.log10(np.abs(x))


def metrics(S: np.ndarray, freq: np.ndarray, band=None) -> dict:
    """Derived performance figures from the modal S-matrix."""
    if band is None:
        sel = np.ones(len(freq), bool)
    else:
        sel = (freq >= band[0]) & (freq <= band[1])

    def wc(x, worst="max"):
        v = x[sel]
        v = v[np.isfinite(v)]
        if v.size == 0:
            return float("nan")
        return float(np.max(v) if worst == "max" else np.min(v))

    polH_il = db20(S[:, 2, 0])  # common polH -> main arm
    side_p = S[:, 3, 1]  # common polV -> side arm +x
    side_m = S[:, 4, 1]  # common polV -> side arm -x
    polV_sum = db20(np.sqrt(np.abs(side_p) ** 2 + np.abs(side_m) ** 2))

    imbal = db20(side_p) - db20(side_m)
    rel_phase = np.angle(side_p / side_m, deg=True)

    out = {
        "f_GHz": freq / 1e9,
        "polH_main_dB": polH_il,
        "polV_side_sum_dB": polV_sum,
        "polV_side_p_dB": db20(side_p),
        "polV_side_m_dB": db20(side_m),
        "side_imbalance_dB": imbal,
        "side_rel_phase_deg": rel_phase,
        "S11_polH_dB": db20(S[:, 0, 0]),
        "S11_polV_dB": db20(S[:, 1, 1]),
        "iso_polV_to_main_dB": db20(S[:, 2, 1]),
        "iso_polH_to_side_p_dB": db20(S[:, 3, 0]),
        "iso_polH_to_side_m_dB": db20(S[:, 4, 0]),
        "xpol_common_dB": db20(S[:, 1, 0]),
        "summary": {},
    }
    out["summary"] = {
        "worst polH insertion loss (dB)": wc(polH_il, "min"),
        "worst polV (both arms) loss (dB)": wc(polV_sum, "min"),
        "worst polH return loss (dB)": wc(out["S11_polH_dB"], "max"),
        "worst polV return loss (dB)": wc(out["S11_polV_dB"], "max"),
        "worst polV->main isolation (dB)": wc(out["iso_polV_to_main_dB"], "max"),
        "worst polH->side isolation (dB)": max(
            wc(out["iso_polH_to_side_p_dB"], "max"),
            wc(out["iso_polH_to_side_m_dB"], "max"),
        ),
        "worst side-arm amplitude imbalance (dB)": wc(np.abs(imbal), "max"),
        "side-arm relative phase, min/max (deg)": (
            wc(rel_phase, "min"),
            wc(rel_phase, "max"),
        ),
    }
    return out


def save(out_dir: str, freq: np.ndarray, S: np.ndarray, p: BoifotParams) -> None:
    os.makedirs(out_dir, exist_ok=True)
    np.savez(
        os.path.join(out_dir, "sparams.npz"),
        freq=freq,
        S=S,
        port_names=np.array(PORT_NAMES),
    )
    p.to_json(os.path.join(out_dir, "params.json"))
    hdr = ["f_Hz"]
    cols = [freq]
    for i in range(NP):
        for j in range(NP):
            hdr += [f"re_S{j+1}{i+1}", f"im_S{j+1}{i+1}"]
            cols += [S[:, j, i].real, S[:, j, i].imag]
    np.savetxt(
        os.path.join(out_dir, "sparams.csv"),
        np.column_stack(cols),
        delimiter=",",
        header=",".join(hdr),
        comments="",
    )


def write_touchstone(path: str, freq: np.ndarray, S: np.ndarray) -> None:
    """Write a 5-port Touchstone (.s5p) file, magnitude/angle format."""
    with open(path, "w") as fh:
        fh.write("! 4-port Boifot orthomode junction, modal S-parameters\n")
        fh.write("! ports: 1=common polH  2=common polV  3=main  4=side+x  5=side-x\n")
        fh.write("# GHz S MA R 1\n")
        for k, f in enumerate(freq):
            fh.write(f"{f/1e9:.6f}")
            for j in range(NP):
                for i in range(NP):
                    v = S[k, j, i]
                    fh.write(f" {abs(v):.6e} {np.degrees(np.angle(v)):.3f}")
                fh.write("\n" + (" " * 10 if j < NP - 1 else ""))
