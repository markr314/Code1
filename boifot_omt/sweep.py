"""Drive the FDTD runs and assemble the modal S-matrix."""
from __future__ import annotations

import os
from typing import Iterable, Optional, Sequence

import numpy as np

from .params import BoifotParams
from .model import (build_model, mesh_summary, PORT_NAMES, PORT_LABELS,
                    port_names, n_ports)
from .waveguide import beta as wg_beta, f_cutoff_te

NP = len(PORT_NAMES)


def port_cutoffs(p: BoifotParams) -> np.ndarray:
    """Cutoff frequency (Hz) of each modal port's operating mode."""
    fc = [p.fc, p.fc, p.fc, p.fc_side, p.fc_side]
    if p.hom_monitor:
        fc.append(f_cutoff_te(p.a_o, p.a_o, 1, 1))
    return np.array(fc)


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
    np_tot = n_ports(p)
    names = port_names(p)
    S = np.full((len(freq), np_tot, np_tot), np.nan, dtype=complex)
    cwd = os.getcwd()
    try:
        for exc in excitations:
            m = build_model(p, excite_port=exc)
            if exc == list(excitations)[0]:
                log(mesh_summary(m))
            sim = os.path.join(os.path.abspath(out_dir), f"exc_{names[exc]}")
            log(f"--- exciting port {exc} ({PORT_LABELS[names[exc]]}) -> {sim}")
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
    nprt = S.shape[1]
    phase = np.zeros((len(freq), nprt))
    for i in range(nprt):
        b = np.array([wg_beta(f, fc[i]) for f in freq])
        phase[:, i] = b * dist_mm[i] * 1e-3
    out = np.empty_like(S)
    for i in range(nprt):
        for j in range(nprt):
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

    main_il = db20(S[:, 2, 0])  # common main-arm pol -> main arm
    side_p = S[:, 3, 1]  # common side-arm pol -> side arm +x
    side_m = S[:, 4, 1]  # common side-arm pol -> side arm -x
    side_sum = db20(np.sqrt(np.abs(side_p) ** 2 + np.abs(side_m) ** 2))

    imbal = db20(side_p) - db20(side_m)
    rel_phase = np.angle(side_p / side_m, deg=True)

    out = {
        "f_GHz": freq / 1e9,
        "main_il_dB": main_il,
        "side_sum_dB": side_sum,
        "side_p_dB": db20(side_p),
        "side_m_dB": db20(side_m),
        "side_imbalance_dB": imbal,
        "side_rel_phase_deg": rel_phase,
        "S11_main_dB": db20(S[:, 0, 0]),
        "S11_side_dB": db20(S[:, 1, 1]),
        "iso_sidepol_to_main_dB": db20(S[:, 2, 1]),
        "iso_mainpol_to_side_p_dB": db20(S[:, 3, 0]),
        "iso_mainpol_to_side_m_dB": db20(S[:, 4, 0]),
        "xpol_common_dB": db20(S[:, 1, 0]),
        "summary": {},
    }
    out["summary"] = {
        "worst main-arm insertion loss (dB)": wc(main_il, "min"),
        "worst side-arm (both) insertion loss (dB)": wc(side_sum, "min"),
        "worst main-arm return loss (dB)": wc(out["S11_main_dB"], "max"),
        "worst side-arm return loss (dB)": wc(out["S11_side_dB"], "max"),
        "worst isolation, side pol -> main arm (dB)": wc(
            out["iso_sidepol_to_main_dB"], "max"
        ),
        "worst isolation, main pol -> side arms (dB)": max(
            wc(out["iso_mainpol_to_side_p_dB"], "max"),
            wc(out["iso_mainpol_to_side_m_dB"], "max"),
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
        fh.write("! ports: 1=common main-pol  2=common side-pol  3=main  4=side+x  5=side-x\n")
        fh.write("# GHz S MA R 1\n")
        for k, f in enumerate(freq):
            fh.write(f"{f/1e9:.6f}")
            for j in range(NP):
                for i in range(NP):
                    v = S[k, j, i]
                    fh.write(f" {abs(v):.6e} {np.degrees(np.angle(v)):.3f}")
                fh.write("\n" + (" " * 10 if j < NP - 1 else ""))


def unitarity(S: np.ndarray) -> np.ndarray:
    """Power sum of each excited column: sum_j |S_ji|^2.

    The junction is lossless (PEC walls, air fill), so every column should sum
    to 1.  This is the acceptance test for a run: it catches truncated
    excitations, energy still trapped in the structure when the run stopped,
    a port whose mode or reference plane is wrong, and power leaving through a
    boundary that should have reflected it.
    """
    return np.nansum(np.abs(S) ** 2, axis=1)


def unitarity_report(S: np.ndarray, freq: np.ndarray, band=None) -> str:
    u = unitarity(S)
    sel = (
        np.ones(len(freq), bool)
        if band is None
        else (freq >= band[0]) & (freq <= band[1])
    )
    lines = ["power conservation, sum_j |S_ji|^2 per excited port i:"]
    names = PORT_NAMES + (["common_TE11"] if u.shape[1] > 5 else [])
    for i in range(u.shape[1]):
        col = u[sel, i]
        if not np.any(np.isfinite(col)) or np.allclose(col, 0):
            lines.append(f"  port {i} ({names[i]:12s}): not excited")
            continue
        lines.append(
            f"  port {i} ({names[i]:12s}): min {col.min():.4f}  "
            f"max {col.max():.4f}  worst error {max(abs(col.min()-1), abs(col.max()-1))*100:.2f} %"
        )
    return "\n".join(lines)


def combine_side_arms(S: np.ndarray):
    """Fold the two side-arm ports into their even/odd combinations.

    The real OMT recombines the side arms in a power combiner, so the
    physically driven side-arm port is one of the two symmetry combinations,
    not either arm on its own.  Returns (S4, labels) for the 4x4 matrix

        0 common main-pol, 1 common side-pol, 2 main arm, 3 side-arm odd mode

    and, as a diagnostic, the even-mode column.  The odd mode is the one the
    junction couples to the common arm; the even mode is the asymmetric
    combination, which breaks the junction's twofold symmetry and therefore
    couples to the overmoded common arm's TE11/TM11, modes these ports do not
    monitor.
    """
    nf = S.shape[0]
    # basis change on both the drive and response side
    T = np.zeros((5, 5))
    T[0, 0] = T[1, 1] = T[2, 2] = 1.0
    r = 1.0 / np.sqrt(2.0)
    T[3, 3], T[3, 4] = r, -r  # odd  (anti-phase) combination
    T[4, 3], T[4, 4] = r, r  # even (in-phase) combination
    out = np.empty_like(S)
    for k in range(nf):
        out[k] = T @ S[k] @ T.T
    labels = ["common_M", "common_S", "main", "side_odd", "side_even"]
    return out, labels
