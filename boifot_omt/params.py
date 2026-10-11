"""Design parameters for a 4-port Boifot orthomode junction.

Topology and axis convention
----------------------------
``z``  common-guide / main-arm axis.  ``x`` horizontal, ``y`` vertical.
Origin is the *junction reference plane*: the upstream edge of the side-arm
coupling apertures.

* Common port: square guide ``a_sq x a_sq`` centred on (x,y)=(0,0), carrying
  both polarisations, entering from -z.
  - ``polH``: E along x.  Transmitted through the septum into the MAIN ARM.
  - ``polV``: E along y.  Shorted by the septum, extracted through the two
    SIDE ARMS.
* Septum: thin metal plate in the x=0 plane (normal along x), spanning the
  full guide height in y.  E of polH is normal to it, so an ideally thin
  septum is transparent to polH; E of polV is tangential and in-plane, so
  polV is short-circuited.  Downstream of the septum the guide is split into
  two half-width channels in which polV is below cutoff.
* Main arm: the guide narrows in x (the E-plane dimension of polH, so cutoff
  is unaffected) from ``a_sq`` to ``b_main`` through ``n_steps`` steps.  That
  narrowing puts polV below cutoff in the main arm.
* Side arms: two rectangular guides propagating along +-x, cross-section
  ``b_side`` (y) x ``a_side`` (z), carrying E along y.  They couple to the
  common guide through apertures in the x=+-a_sq/2 walls spanning
  y in [-b_side/2, b_side/2], z in [0, a_side].
* Pins: metal posts with axis along y, standing in each aperture plane, which
  restore a low-impedance return path for the main-arm wall currents across
  the aperture and tune the side-arm reactance.

All lengths are in mm, frequencies in Hz.
"""
from __future__ import annotations

from dataclasses import dataclass, field, asdict
from typing import List, Tuple
import json

from .waveguide import WR_SIZES, f_cutoff_te, guide_wavelength


@dataclass
class BoifotParams:
    # ---- band -----------------------------------------------------------
    f_lo: float = 26.5e9
    f_hi: float = 40.0e9

    # ---- common guide ---------------------------------------------------
    a_sq: float = 7.112  # square common guide side

    # ---- main arm -------------------------------------------------------
    b_main: float = 3.556  # final main-arm dimension along x (E-plane of polH)
    n_steps: int = 2  # number of E-plane matching steps into the main arm
    step_len: float = 2.2  # length of each matching step along z
    main_len: float = 15.0  # straight main-arm run after the last step

    # ---- side arms ------------------------------------------------------
    a_side: float = 7.112  # side-arm broad dimension, along z (sets cutoff)
    b_side: float = 3.556  # side-arm narrow dimension, along y
    side_len: float = 14.0  # straight side-arm run from the aperture outward

    # ---- septum ---------------------------------------------------------
    t_sep: float = 0.18  # septum thickness along x (0 => zero-thickness sheet)
    z_sep_start: float = 0.0  # septum leading edge (junction reference plane)
    l_sep: float = 11.0  # septum length along z

    # ---- compensation pins ---------------------------------------------
    d_pin: float = 0.60  # pin diameter
    pin_z_frac: Tuple[float, ...] = (0.18, 0.82)  # pin centres, fraction of a_side
    pins_enabled: bool = True

    # ---- input guide ----------------------------------------------------
    common_len: float = 12.0  # straight square-guide run upstream of z=0

    # ---- numerics -------------------------------------------------------
    mesh_res: float = 0.25  # target cell size in the guides
    mesh_res_fine: float = 0.09  # also sets the cells across the septum  # cell size near septum / pins
    pml_cells: int = 8
    nr_timesteps: int = 60000
    end_criteria: float = 1e-5

    # ---- derived / validation ------------------------------------------
    def __post_init__(self) -> None:
        if self.t_sep < 0:
            raise ValueError("t_sep must be >= 0")
        if self.b_main <= 0 or self.b_main > self.a_sq:
            raise ValueError("b_main must be in (0, a_sq]")
        if self.pins_enabled and self.d_pin >= self.a_side * min(
            min(self.pin_z_frac), 1 - max(self.pin_z_frac)
        ) * 2:
            # pins must fit inside the aperture in z
            raise ValueError("pin diameter too large for the given pin_z_frac")

    # ---- geometric landmarks -------------------------------------------
    @property
    def z_ap0(self) -> float:
        """Upstream edge of the side-arm apertures."""
        return 0.0

    @property
    def z_ap1(self) -> float:
        """Downstream edge of the side-arm apertures."""
        return self.a_side

    @property
    def z_step0(self) -> float:
        """Start of the main-arm matching steps."""
        return self.z_ap1

    @property
    def step_widths(self) -> List[float]:
        """x half-widths of each matching step, from a_sq/2 down to b_main/2.

        A geometric (Chebyshev-like) progression between the two, which is the
        usual starting point for an E-plane stepped transformer; the step
        lengths are a quarter guide wavelength at band centre by default.
        """
        w0 = self.a_sq / 2.0
        w1 = self.b_main / 2.0
        n = self.n_steps
        if n <= 0:
            return []
        # geometric interpolation of the half-widths
        return [w0 * (w1 / w0) ** ((i + 1) / (n + 1)) for i in range(n)]

    @property
    def z_main_start(self) -> float:
        """z where the full-width main arm (b_main) begins."""
        return self.z_step0 + self.n_steps * self.step_len

    @property
    def z_max(self) -> float:
        return self.z_main_start + self.main_len

    @property
    def z_min(self) -> float:
        return -self.common_len

    @property
    def x_side_end(self) -> float:
        return self.a_sq / 2.0 + self.side_len

    @property
    def z_sep_end(self) -> float:
        return self.z_sep_start + self.l_sep

    # ---- electrical checks ---------------------------------------------
    @property
    def fc_common(self) -> float:
        """TE10/TE01 cutoff of the square common guide (Hz)."""
        return f_cutoff_te(self.a_sq, self.a_sq, 1, 0)

    @property
    def fc_common_te11(self) -> float:
        """First higher-order (TE11/TM11) cutoff of the square guide (Hz)."""
        return f_cutoff_te(self.a_sq, self.a_sq, 1, 1)

    @property
    def fc_main(self) -> float:
        """Main-arm operating-mode cutoff: set by the y dimension (Hz)."""
        return f_cutoff_te(self.a_sq, self.b_main, 1, 0)

    @property
    def fc_main_reject(self) -> float:
        """Cutoff of the rejected polarisation in the main arm (Hz)."""
        return f_cutoff_te(self.b_main, self.a_sq, 1, 0)

    @property
    def fc_side(self) -> float:
        """Side-arm operating-mode cutoff: set by a_side (Hz)."""
        return f_cutoff_te(self.a_side, self.b_side, 1, 0)

    @property
    def fc_halfguide(self) -> float:
        """Cutoff for polV in a septum-split half of the common guide (Hz).

        polV must be *below* cutoff here for the septum to reject it.
        """
        half = (self.a_sq - self.t_sep) / 2.0
        return f_cutoff_te(half, self.a_sq, 1, 0)

    def report(self) -> str:
        lines = [
            "Boifot junction, electrical summary",
            f"  band                      : {self.f_lo/1e9:.2f} - {self.f_hi/1e9:.2f} GHz",
            f"  common guide              : {self.a_sq:.3f} x {self.a_sq:.3f} mm",
            f"    TE10/TE01 cutoff        : {self.fc_common/1e9:.2f} GHz  (want < f_lo)",
            f"    TE11/TM11 cutoff        : {self.fc_common_te11/1e9:.2f} GHz",
            f"  main arm                  : {self.b_main:.3f} (x) x {self.a_sq:.3f} (y) mm",
            f"    operating cutoff        : {self.fc_main/1e9:.2f} GHz  (want < f_lo)",
            f"    rejected-pol cutoff     : {self.fc_main_reject/1e9:.2f} GHz  (want > f_hi)",
            f"  side arms                 : {self.b_side:.3f} (y) x {self.a_side:.3f} (z) mm",
            f"    operating cutoff        : {self.fc_side/1e9:.2f} GHz  (want < f_lo)",
            f"  septum                    : t={self.t_sep:.3f} mm, length={self.l_sep:.2f} mm",
            f"    split-half cutoff (polV): {self.fc_halfguide/1e9:.2f} GHz  (want > f_hi)",
            f"  main-arm steps            : {self.n_steps} x {self.step_len:.2f} mm, "
            f"half-widths {['%.3f' % w for w in self.step_widths]}",
            f"  pins                      : {'on' if self.pins_enabled else 'off'}, "
            f"d={self.d_pin:.2f} mm at z/a_side={self.pin_z_frac}",
        ]
        lg_lo = guide_wavelength(self.f_lo, self.fc_main) * 1e3
        lg_hi = guide_wavelength(self.f_hi, self.fc_main) * 1e3
        lines.append(
            f"  main-arm guide wavelength : {lg_lo:.2f} mm @ f_lo, {lg_hi:.2f} mm @ f_hi"
        )
        lines.append(
            f"  quarter-wave at f_centre  : "
            f"{guide_wavelength(0.5*(self.f_lo+self.f_hi), self.fc_main)*1e3/4:.2f} mm"
        )
        return "\n".join(lines)

    def warnings(self) -> List[str]:
        w = []
        if self.fc_common >= self.f_lo:
            w.append(
                f"common guide cutoff {self.fc_common/1e9:.2f} GHz is above f_lo "
                f"{self.f_lo/1e9:.2f} GHz"
            )
        if self.fc_common_te11 < self.f_hi:
            w.append(
                f"common square guide is overmoded above "
                f"{self.fc_common_te11/1e9:.2f} GHz (TE11/TM11); the junction's "
                "two-fold symmetry is what keeps those modes unexcited"
            )
        if self.fc_main_reject <= self.f_hi:
            w.append(
                f"main arm does not cut off the rejected polarisation over the "
                f"whole band (cutoff {self.fc_main_reject/1e9:.2f} GHz <= f_hi)"
            )
        if self.fc_halfguide <= self.f_hi:
            w.append(
                f"septum-split half-guide does not cut off polV over the whole "
                f"band (cutoff {self.fc_halfguide/1e9:.2f} GHz <= f_hi)"
            )
        if self.fc_side >= self.f_lo:
            w.append(f"side-arm cutoff {self.fc_side/1e9:.2f} GHz is above f_lo")
        if self.z_sep_end > self.z_max - 2.0:
            w.append(
                "septum reaches close to the main-arm port plane; the port mode "
                "profile assumes an unobstructed rectangle"
            )
        return w

    def to_json(self, path: str) -> None:
        with open(path, "w") as fh:
            json.dump(asdict(self), fh, indent=2)


def wr28_default() -> BoifotParams:
    """Ka-band (WR-28, 26.5-40 GHz) starting design."""
    wr = WR_SIZES["WR-28"]
    return BoifotParams(
        f_lo=wr.f_lo * 1e9,
        f_hi=wr.f_hi * 1e9,
        a_sq=wr.a,
        b_main=wr.b,
        a_side=wr.a,
        b_side=wr.b,
    )
