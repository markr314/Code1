"""Design parameters for a 4-port Bøifot orthomode junction.

Geometry follows Wollack, Grammer & Kingsley, "The Bøifot Orthomode
Junction", ALMA Memo 425 (2002).

Axis convention
---------------
``z``   common-arm / main-arm axis, wave entering from -z.
``y``   normal to the split-block plane; the septum lies in y=0.
``x``   towards the two side arms.
The origin is the septum tip apex.

How the junction works
----------------------
The junction is a turnstile junction with two of its four side ports folded
parallel to the common arm (memo, Fig. 1).  Both polarisations leave as
*pairs* of channels:

* **Main arm**, E parallel to y, i.e. normal to the septum.  An ideally thin
  septum is invisible to it, so it passes as two half-height channels that
  recombine; the septum tip tapers the guide adiabatically from 1:1 to
  (1 - dt/a_o):1, and a four-step 2:1 Chebyshev transformer then takes the
  full height from ``a_o`` down to standard ``b_o = a_o/2``.
* **Side arms**, E parallel to x, i.e. lying in the septum plane.  The septum
  shorts this polarisation, and its mitred leading edge acts as two
  back-to-back miters that divide the power into the two side arms through
  full-height apertures in the x = +-a_o/2 walls.  Downstream of the septum
  tip this polarisation is below cutoff in the half-height channels, so the
  main arm is a reflective termination for it.

The split-block plane (y=0) is simultaneously the main arm's H-plane and the
side arms' E-plane, which is why the septum is a shim clamped between the two
block halves.

Pins
----
Posts along y standing in the common arm, flush with the side-arm wall and
centred in the side-arm aperture.  For the main-arm polarisation they are
short-circuit stubs that tune out the side-arm discontinuity and restore a
low-impedance return path for the main-arm currents; for the side-arm
polarisation they load the aperture.  The memo's ranges:

    0.02 < dt_septum / a_o < 0.05        septum thickness
    0.05 < d_pin / a_o     < 0.07        pin diameter

All lengths are in mm, frequencies in Hz.
"""
from __future__ import annotations

from dataclasses import dataclass, asdict
from typing import List, Optional, Tuple
import json

from .waveguide import (
    WR_SIZES,
    f_cutoff_te,
    guide_wavelength,
    chebyshev_transformer,
)


@dataclass
class BoifotParams:
    # ---- band -----------------------------------------------------------
    f_lo: float = 26.5e9
    f_hi: float = 40.0e9

    # ---- guide sizes ----------------------------------------------------
    a_o: float = 7.112  # broadwall / square common-arm side
    b_o: float = 3.556  # standard guide height, a_o/2

    # ---- septum (memo: 0.02 < dt/a_o < 0.05) ----------------------------
    dt_sep_norm: float = 0.035
    miter_slope: float = 1.0  # dz/dx of the tip edge; 1.0 is a 45 deg miter
    # Septum tip (the memo's Sigma_1 / Sigma_2, Fig. 1): a blunt leading
    # finger ahead of the mitred wedge.  Its length and width set the
    # upper/lower band-edge resonances of the double-tuned side-arm match.
    tip_len: float = 0.0  # finger length upstream of the aperture plane
    tip_w: float = 0.0  # finger half-width in x
    sep_overhang: float = 2.0  # septum length past the end of the transformer
    # Septum tip thickness taper.  The memo: "the septum tip forms an
    # adiabatic impedance transformer from 1:1 guide to (1 - dt/a_o):1".
    # Upstream of the septum the main-arm polarisation sees the full height
    # a_o; where the septum is at full thickness it sees a_o - dt.  A blunt
    # tip steps between those, a tapered tip transitions gradually.  The
    # taper runs `tip_taper_len` along z *behind the mitred leading edge*, so
    # it follows the edge contour rather than a plane, and is built as
    # `tip_taper_steps` stacked layers.  0 disables it (a blunt tip).
    tip_taper_len: float = 0.0
    tip_taper_steps: int = 3

    # ---- side-arm aperture ---------------------------------------------
    # full height in y; the z extent is the side-arm guide height
    ap_len_norm: float = 1.0  # aperture z extent, in units of b_o
    ap_from_tip: bool = False  # alternative: let the tip drag the aperture
    tip_auto_slope: bool = True  # hold the aperture fixed, see `slope`
    side_len: float = 14.0  # straight side-arm run from the aperture outward

    # ---- main-arm transformer (memo: 4-step, 2:1, w_q=0.84, f_o/f_c=1.565)
    n_steps: int = 4
    transformer_ratio: float = 2.0
    w_q: float = 0.84
    fo_over_fc: float = 1.565
    trans_gap: float = 0.0  # straight square guide between tip and transformer
    # Length of the transformer's FIRST section, when it should differ from
    # the quarter-wave value.  Ideal Chebyshev theory treats every junction as
    # a perfect impedance step, but each step carries a shunt susceptance, and
    # the first section is the odd one out: it does not face a plain guide, it
    # faces the junction itself (aperture, septum tip, pins).  It is also the
    # phase knob on the reflective termination that the side-arm polarisation
    # looks into.  None keeps it equal to the others.
    first_step_len: float = None
    main_len: float = 14.0  # straight main-arm run after the transformer

    # ---- pins (memo: 0.05 < d_pin/a_o < 0.07) ---------------------------
    pin_d_norm: float = 0.06
    n_pins: int = 1  # pins per aperture
    pin_inset: float = 0.0  # extra |x| offset inward from the flush position
    # "full": one post spanning the whole common-arm height in y.
    # "half": a pair of posts per location, each running from the septum
    # plane out to a broad wall.  These are the same element wherever the
    # septum is absent, and different wherever it is present, which is the
    # open question about how the memo's pins are built.
    pin_span: str = "full"
    # Explicit pin centres as fractions of the aperture length along z.
    # None spreads n_pins symmetrically about the aperture centre.
    pin_z_frac: Optional[Tuple[float, ...]] = None
    pins_enabled: bool = True
    side_arms_enabled: bool = True  # diagnostic: blank off the side arms

    # ---- common arm -----------------------------------------------------
    common_len: float = 12.0  # straight square guide upstream of the tip

    # ---- numerics -------------------------------------------------------
    mesh_res: float = 0.25
    mesh_res_junction: float = 0.15  # junction region, where the mitred edge is staircased
    mesh_res_fine: float = 0.09
    hom_monitor: bool = False  # add a TE11 monitor port on the common arm
    pml_cells: int = 8
    nr_timesteps: int = 60000
    end_criteria: float = 1e-5

    def __post_init__(self) -> None:
        if not 0.0 <= self.dt_sep_norm < 0.5:
            raise ValueError("dt_sep_norm out of range")
        if self.n_pins < 0:
            raise ValueError("n_pins must be >= 0")
        if not 0.0 <= self.tip_w < self.a_o / 2.0:
            raise ValueError("tip_w must be in [0, a_o/2)")
        if self.pin_span not in ("full", "half"):
            raise ValueError("pin_span must be 'full' or 'half'")
        if self.tip_taper_len < 0:
            raise ValueError("tip_taper_len must be >= 0")
        if self.tip_taper_steps < 1:
            raise ValueError("tip_taper_steps must be >= 1")

    # ---- derived dimensions --------------------------------------------
    @property
    def dt_sep(self) -> float:
        return self.dt_sep_norm * self.a_o

    @property
    def d_pin(self) -> float:
        return self.pin_d_norm * self.a_o

    @property
    def slope(self) -> float:
        """Effective dz/dx of the mitred tip edge.

        With `tip_auto_slope` the aperture stays fixed at the standard
        side-arm height and the slope is whatever makes the wedge reach the
        side walls exactly at the aperture's downstream edge.  That keeps
        tip_len and tip_w as clean one-variable knobs: without it, blunting
        the tip also shifts the aperture upstream and shortens it, so the
        side-arm guide height moves at the same time and a tuning sweep
        confounds three changes at once.  tip_len = tip_w = 0 gives 45 deg.
        """
        if self.tip_auto_slope:
            return (self.ap_len + self.tip_len) / (self.a_o / 2.0 - self.tip_w)
        return self.miter_slope

    @property
    def ap_len(self) -> float:
        """Side-arm aperture extent along z (= side-arm guide height).

        With `ap_from_tip`, the aperture runs from the septum's leading edge
        to the point where the mitred wedge reaches the side walls, so the
        septum only ever spans the full guide width *downstream* of the
        aperture.  Letting the wedge reach full width inside the aperture
        splits it in two in y, puts both halves below cutoff and reflects the
        side-arm polarisation almost completely.  With tip_len = tip_w = 0
        this reduces to a_o/2 = b_o.
        """
        if self.ap_from_tip and not self.tip_auto_slope:
            return self.ap_len_norm * self.miter_slope * (
                self.a_o / 2.0 - self.tip_w
            )
        return self.ap_len_norm * self.b_o

    @property
    def z_tip(self) -> float:
        """Leading edge of the septum tip (upstream of the aperture plane)."""
        return -self.tip_len

    @property
    def taper_layers(self) -> List[Tuple[float, float]]:
        """(half_thickness, z_offset) per layer of the tapered tip.

        Layer k carries thickness dt*k/N and begins dt_len*k/N behind the
        leading edge, so their union staircases a linear thickness ramp from
        a knife edge up to the full septum.  An empty list means a blunt tip.
        """
        if self.tip_taper_len <= 0 or self.dt_sep <= 0:
            return []
        n = self.tip_taper_steps
        if n == 1:
            return [(0.5 * self.dt_sep, 0.0)]
        # thinnest layer starts at the apex so the tip is a knife edge, not a
        # gap; the full-thickness layer starts exactly tip_taper_len behind it
        return [
            (0.5 * self.dt_sep * k / n, self.tip_taper_len * (k - 1) / (n - 1))
            for k in range(1, n + 1)
        ]

    @property
    def z_tip_end(self) -> float:
        """z at which the mitred tip edge reaches the side walls."""
        return -self.tip_len + self.slope * (self.a_o / 2.0 - self.tip_w)

    @property
    def z_ap0(self) -> float:
        return self.z_tip if (self.ap_from_tip and not self.tip_auto_slope) else 0.0

    @property
    def z_ap1(self) -> float:
        return self.z_ap0 + self.ap_len

    @property
    def fc(self) -> float:
        """Common-arm / main-arm cutoff, set by a_o (Hz)."""
        return f_cutoff_te(self.a_o, self.a_o, 1, 0)

    @property
    def f_o(self) -> float:
        """Transformer synchronous frequency (Hz)."""
        return self.fo_over_fc * self.fc

    @property
    def step_heights(self) -> List[float]:
        """Full guide height of each transformer section, a_o -> b_o."""
        Z, _ = chebyshev_transformer(self.transformer_ratio, self.n_steps, self.w_q)
        return [self.a_o / z for z in Z[1:-1]]

    @property
    def step_ripple_db(self) -> float:
        import math

        _, gm = chebyshev_transformer(self.transformer_ratio, self.n_steps, self.w_q)
        return -20.0 * math.log10(gm)

    @property
    def step_len(self) -> float:
        """Quarter guide wavelength at the synchronous frequency (mm).

        E-plane steps do not move the cutoff, so every section has the same
        guide wavelength.
        """
        return guide_wavelength(self.f_o, self.fc) * 1e3 / 4.0

    @property
    def z_trans0(self) -> float:
        return self.z_tip_end + self.trans_gap

    @property
    def step_lengths(self) -> List[float]:
        """Length of each transformer section, first one possibly different."""
        L = [self.step_len] * self.n_steps
        if self.first_step_len is not None and self.n_steps > 0:
            L[0] = self.first_step_len
        return L

    @property
    def z_trans1(self) -> float:
        return self.z_trans0 + sum(self.step_lengths)

    @property
    def z_sep_end(self) -> float:
        return self.z_trans1 + self.sep_overhang

    @property
    def z_min(self) -> float:
        return -self.common_len

    @property
    def z_max(self) -> float:
        return self.z_trans1 + self.main_len

    @property
    def x_side_end(self) -> float:
        return self.a_o / 2.0 + self.side_len

    @property
    def pin_x(self) -> float:
        """Pin centre |x|: outer diameter flush with the side-arm wall."""
        return self.a_o / 2.0 - self.d_pin / 2.0 - self.pin_inset

    @property
    def pin_z(self) -> List[float]:
        """Pin centres along z, spread symmetrically about the aperture centre."""
        if not self.pins_enabled or self.n_pins <= 0:
            return []
        if self.pin_z_frac is not None:
            return [self.z_ap0 + f * self.ap_len for f in self.pin_z_frac]
        c = 0.5 * (self.z_ap0 + self.z_ap1)
        if self.n_pins == 1:
            return [c]
        pitch = (self.ap_len - self.d_pin) / self.n_pins
        k0 = (self.n_pins - 1) / 2.0
        return [c + (i - k0) * pitch for i in range(self.n_pins)]

    # ---- electrical checks ---------------------------------------------
    @property
    def fc_side(self) -> float:
        """Side-arm cutoff: broad dimension is the full height a_o (Hz)."""
        return f_cutoff_te(self.a_o, self.ap_len, 1, 0)

    @property
    def fc_halfheight(self) -> float:
        """Side-arm pol cutoff in a septum-split half channel (Hz).

        Must sit above the band so the main arm is a reflective termination
        for the side-arm polarisation.
        """
        half = (self.a_o - self.dt_sep) / 2.0
        return f_cutoff_te(half, self.a_o, 1, 0)

    @property
    def fc_main_reject(self) -> float:
        """Side-arm pol cutoff in the finished main arm (Hz)."""
        return f_cutoff_te(self.b_o, self.a_o, 1, 0)

    def report(self) -> str:
        L = [
            "Boifot orthomode junction (ALMA Memo 425 geometry)",
            f"  band                      : {self.f_lo/1e9:.2f} - {self.f_hi/1e9:.2f} GHz"
            f"  ({self.f_lo/self.fc:.2f} - {self.f_hi/self.fc:.2f} f_c)",
            f"  common arm                : {self.a_o:.3f} x {self.a_o:.3f} mm,"
            f" f_c = {self.fc/1e9:.2f} GHz",
            f"  main arm out              : {self.a_o:.3f} (x) x {self.b_o:.3f} (y) mm",
            f"  side arms                 : {self.a_o:.3f} (y) x {self.ap_len:.3f} (z) mm,"
            f" f_c = {self.fc_side/1e9:.2f} GHz",
            f"  septum                    : dt = {self.dt_sep:.3f} mm"
            f" (dt/a_o = {self.dt_sep_norm:.3f}), tip edge reaches the wall at"
            f" z = {self.z_tip_end:.3f} mm",
            f"  septum ends at            : z = {self.z_sep_end:.2f} mm",
            f"  septum tip                : finger {self.tip_len:.3f} mm long,"
            f" half-width {self.tip_w:.3f} mm, leading edge z = {self.z_tip:.3f} mm",
            f"  septum tip taper          : "
            + (
                f"{self.tip_taper_len:.3f} mm over {self.tip_taper_steps} layers, "
                + "thicknesses "
                + ", ".join("%.4f" % (2 * h) for h, _ in self.taper_layers)
                + " mm"
                if self.taper_layers
                else "none (blunt tip)"
            ),
            f"  side-arm aperture         : z = {self.z_ap0:.2f} .. {self.z_ap1:.2f} mm,"
            f" full height in y",
            f"  transformer               : {self.n_steps}-step,"
            f" {self.transformer_ratio:.0f}:1 Chebyshev, w_q = {self.w_q:.2f},"
            f" f_o/f_c = {self.fo_over_fc:.3f}",
            f"    section lengths         : "
            + ", ".join("%.3f" % L for L in self.step_lengths) + " mm",
            f"    synchronous frequency   : {self.f_o/1e9:.2f} GHz,"
            f" section length {self.step_len:.3f} mm",
            f"    section heights         : "
            + ", ".join("%.3f" % h for h in self.step_heights)
            + f"  (-> {self.b_o:.3f})",
            f"    design ripple           : {self.step_ripple_db:.1f} dB return loss",
            f"  pins                      : {self.n_pins} per aperture,"
            f" d = {self.d_pin:.3f} mm (d/a_o = {self.pin_d_norm:.3f}),"
            f" |x| = {self.pin_x:.3f} mm",
            f"    z positions             : "
            + ", ".join("%.3f" % z for z in self.pin_z),
            f"  side-arm pol cutoff in half-channel : {self.fc_halfheight/1e9:.2f} GHz"
            f"  (want > f_hi)",
            f"  side-arm pol cutoff in main arm     : {self.fc_main_reject/1e9:.2f} GHz"
            f"  (want > f_hi)",
        ]
        return "\n".join(L)

    def warnings(self) -> List[str]:
        w = []
        if not 0.02 <= self.dt_sep_norm <= 0.05:
            w.append(
                f"septum thickness dt/a_o = {self.dt_sep_norm:.3f} is outside the "
                "0.02-0.05 range the memo reports as compensatable"
            )
        if self.pins_enabled and not 0.05 <= self.pin_d_norm <= 0.07:
            w.append(
                f"pin diameter d/a_o = {self.pin_d_norm:.3f} is outside the memo's "
                "0.05-0.07 range"
            )
        if self.f_lo / self.fc < 1.2 or self.f_hi / self.fc > 1.9:
            w.append(
                f"band is {self.f_lo/self.fc:.2f}-{self.f_hi/self.fc:.2f} f_c; the "
                "memo designs for 1.2-1.9 f_c"
            )
        if self.fc_halfheight <= self.f_hi:
            w.append(
                f"septum half-channel cutoff {self.fc_halfheight/1e9:.2f} GHz is not "
                "above the band; the main arm will not reflect the side-arm pol"
            )
        if self.fc_main_reject <= self.f_hi:
            w.append(
                f"main-arm rejection cutoff {self.fc_main_reject/1e9:.2f} GHz is not "
                "above the band"
            )
        if self.pin_span == "half" and self.pins_enabled:
            for zc in self.pin_z:
                z_sep_here = self.z_tip + self.slope * max(
                    abs(self.pin_x) - self.tip_w, 0.0
                )
                if not (z_sep_here <= zc <= self.z_sep_end):
                    w.append(
                        f"half-height pin at z={zc:.3f} mm has no septum to land "
                        f"on at |x|={self.pin_x:.3f} mm (septum starts at "
                        f"z={z_sep_here:.3f}); its inner end is unsupported"
                    )
        if self.fc_side >= self.f_lo:
            w.append(f"side-arm cutoff {self.fc_side/1e9:.2f} GHz is above f_lo")
        if self.z_tip_end > self.z_ap1 + 1e-9:
            w.append(
                "the septum reaches full width inside the side-arm aperture; "
                "that splits the aperture in y and reflects the side-arm "
                "polarisation"
            )
        return w

    def to_json(self, path: str) -> None:
        d = asdict(self)
        d["_derived"] = {
            "dt_sep_mm": self.dt_sep,
            "d_pin_mm": self.d_pin,
            "step_heights_mm": self.step_heights,
            "step_len_mm": self.step_len,
            "f_c_Hz": self.fc,
            "f_o_Hz": self.f_o,
        }
        with open(path, "w") as fh:
            json.dump(d, fh, indent=2)


def wr28_default() -> BoifotParams:
    """Ka-band (WR-28, 26.5-40 GHz) Boifot junction, memo geometry."""
    wr = WR_SIZES["WR-28"]
    return BoifotParams(f_lo=wr.f_lo * 1e9, f_hi=wr.f_hi * 1e9, a_o=wr.a, b_o=wr.b)
