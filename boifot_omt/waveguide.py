"""Analytic rectangular-waveguide relations and a small WR-size table.

Used for sizing the Boifot junction arms, for choosing the FDTD mesh, and as
the reference against which the openEMS model is validated.
"""
from __future__ import annotations

from dataclasses import dataclass
import math

C0 = 299792458.0  # m/s, exact
MU0 = 4e-7 * math.pi
EPS0 = 1.0 / (MU0 * C0 * C0)
Z0 = math.sqrt(MU0 / EPS0)  # ~376.73 ohm


@dataclass(frozen=True)
class WRSize:
    """A standard rectangular waveguide size. Dimensions in mm."""

    name: str
    a: float  # broad wall (mm)
    b: float  # narrow wall (mm)
    f_lo: float  # recommended band low edge (GHz)
    f_hi: float  # recommended band high edge (GHz)


# Inner dimensions, EIA WR designations.
WR_SIZES = {
    "WR-90": WRSize("WR-90", 22.860, 10.160, 8.20, 12.40),
    "WR-42": WRSize("WR-42", 10.668, 4.318, 18.00, 26.50),
    "WR-28": WRSize("WR-28", 7.112, 3.556, 26.50, 40.00),
    "WR-22": WRSize("WR-22", 5.690, 2.845, 33.00, 50.00),
    "WR-15": WRSize("WR-15", 3.759, 1.880, 50.00, 75.00),
    "WR-10": WRSize("WR-10", 2.540, 1.270, 75.00, 110.00),
}


def f_cutoff_te(a_mm: float, b_mm: float, m: int = 1, n: int = 0) -> float:
    """TE_mn cutoff frequency in Hz for an air-filled guide.

    `a_mm` is the dimension along which index `m` varies, `b_mm` along `n`.
    """
    a = a_mm * 1e-3
    b = b_mm * 1e-3
    terms = (m / a) ** 2 if a > 0 else 0.0
    terms += (n / b) ** 2 if b > 0 else 0.0
    return 0.5 * C0 * math.sqrt(terms)


def beta(f_hz: float, fc_hz: float) -> float:
    """Propagation constant (rad/m). Negative-imaginary branch is returned as 0
    for evanescent frequencies; use `is_propagating` to test first."""
    if f_hz <= fc_hz:
        return 0.0
    k0 = 2 * math.pi * f_hz / C0
    kc = 2 * math.pi * fc_hz / C0
    return math.sqrt(k0 * k0 - kc * kc)


def guide_wavelength(f_hz: float, fc_hz: float) -> float:
    """Guide wavelength (m)."""
    b = beta(f_hz, fc_hz)
    return float("inf") if b == 0.0 else 2 * math.pi / b


def wave_impedance_te(f_hz: float, fc_hz: float) -> float:
    """TE-mode wave impedance (ohm). Infinite at cutoff."""
    if f_hz <= fc_hz:
        return float("inf")
    return Z0 / math.sqrt(1.0 - (fc_hz / f_hz) ** 2)


def attenuation_db_per_m(f_hz: float, a_mm: float, b_mm: float, sigma: float) -> float:
    """TE10 conductor loss in dB/m for a guide of broad wall `a_mm`.

    Standard Marcuvitz/Pozar result; `sigma` in S/m (Al ~3.5e7, Cu ~5.8e7).
    """
    a = a_mm * 1e-3
    b = b_mm * 1e-3
    fc = f_cutoff_te(a_mm, b_mm, 1, 0)
    if f_hz <= fc:
        return float("inf")
    omega = 2 * math.pi * f_hz
    Rs = math.sqrt(omega * MU0 / (2 * sigma))
    ratio = (fc / f_hz) ** 2
    # alpha_c for TE10 (Np/m)
    alpha = (Rs / (a ** 3 * b * beta(f_hz, fc) * omega * MU0)) * (
        2 * b * math.pi ** 2 + a ** 3 * (beta(f_hz, fc) ** 2)
    )
    del ratio  # kept for readability of the closed form above
    return alpha * 8.685889638


def is_propagating(f_hz: float, fc_hz: float) -> bool:
    return f_hz > fc_hz
