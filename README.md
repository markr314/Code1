# 4-port Bøifot orthomode junction in openEMS (WR-28, 26.5–40 GHz)

A parameterized openEMS (FDTD) model of a 4-port Bøifot orthomode junction,
built to the geometry described in E.J. Wollack, W. Grammer and J. Kingsley,
*The Bøifot Orthomode Junction*, ALMA Memo 425 (2002), and the harness that
extracts its modal S-parameters.

## How the junction works

The Bøifot junction is a turnstile junction with two of its four side ports
folded parallel to the common arm (memo, Fig. 1). **Both** polarizations leave
as *pairs* of channels, which is what makes it twofold symmetric and therefore
broadband — the symmetry keeps TE11/TM11 from being excited in the overmoded
square common arm.

```
        y  (normal to the split plane; the septum lies in y = 0)
        |      side arm  <--[ aperture ]   [ aperture ]--> side arm
        |                        \   junction   /
        +---- x                   \           /
       /                   common arm ===> main arm (4-step 2:1 transformer)
      z  (common-arm / main-arm axis)
```

* **Main arm** — E ∥ y, normal to the septum. A thin septum is invisible to
  it, so it passes as two half-height channels that recombine. The septum tip
  tapers the guide adiabatically from 1:1 to (1 − dt/a₀):1, then a four-step
  2:1 Chebyshev transformer takes the full height from a₀ down to standard
  b₀ = a₀/2.
* **Side arms** — E ∥ x, lying in the septum plane. The septum shorts this
  polarization, and its **mitred leading edge acts as two back-to-back miters**
  that divide the power between the two side arms through full-height
  apertures in the x = ±a₀/2 walls. Downstream of the tip this polarization is
  below cutoff in the half-height channels, so the main arm is a reflective
  termination for it.
* **Pins** — posts along y standing in the common arm, flush with the side-arm
  wall and centred in the aperture. For the main-arm polarization they are
  short-circuit stubs that tune out the side-arm discontinuity and restore a
  low-impedance return path for the main-arm currents.

The split-block plane (y = 0) is simultaneously the main arm's H-plane and the
side arms' E-plane — which is why the septum is a shim clamped between the two
block halves, and why the memo notes the main-arm loss is sensitive to the
quality of the block mate while the side-arm loss is not.

### The constraint that makes the geometry non-obvious

A polarization couples strongly out of a wall only where its E-field is
*normal* to that wall, while the septum shorts whichever polarization lies
*in* its plane. Those two requirements pull against each other, and the
mitred septum tip is what reconciles them: the side-arm polarization has its
E-field normal to the aperture walls and maximal right at the septum plane,
exactly where the miter steers it outward.

## Ports

The common port is dual-mode, so the 4-port junction is a 5×5 modal S-matrix:

| idx | port | mode |
|-----|------|------|
| 0 | common arm, main-arm pol | E ∥ y |
| 1 | common arm, side-arm pol | E ∥ x |
| 2 | main arm | E ∥ y |
| 3 | side arm, +x | E ∥ z, broad dimension along y |
| 4 | side arm, −x | E ∥ z, broad dimension along y |

The side arms are brought out separately; combining them (the memo's E-plane
power combiner, its Table 1 case "B") is outside this model.

## Design values (WR-28)

Scaled from the memo's normalized rules. Its design band is 1.2–1.9 f_c, and
WR-28 at 26.5–40 GHz is 1.26–1.90 f_c — a near-exact fit.

| quantity | value | memo rule |
|---|---|---|
| a₀ (square common arm, broadwall) | 7.112 mm | WR-28 |
| f_c | 21.077 GHz | |
| b₀ (main arm, side arms) | 3.556 mm | a₀/2 |
| septum thickness dt | 0.249 mm | 0.02 < dt/a₀ < 0.05 → 0.035 |
| pin diameter | 0.427 mm | 0.05 < d/a₀ < 0.07 → 0.060 |
| pins per aperture | 1, flush to the side-arm wall, centred | memo's one-pair optimum |
| septum tip miter | 45°, reaching the side walls at z = 3.556 mm | |
| side-arm aperture | full height in y × 3.556 mm in z | |
| main-arm transformer | 4 steps × 2.954 mm | 2:1 Chebyshev, w_q = 0.84, f_o/f_c = 1.565 |
| section heights | 6.648, 5.617, 4.503, 3.804 mm | → 3.556 |
| transformer design ripple | 40.4 dB return loss | |

The 45° miter and the aperture length are self-consistent: the tip edge
reaches the side wall at exactly z = a₀/2 = b₀, which is the aperture's own
length along z.

## Layout

```
boifot_omt/waveguide.py       WR table, cutoff/beta/impedance, Chebyshev synthesis
boifot_omt/params.py          every dimension in one dataclass, with electrical checks
boifot_omt/geometry.py        CSXCAD primitives and the critical mesh lines
boifot_omt/model.py           FDTD setup, boundaries, the five modal ports
boifot_omt/check_geometry.py  solver-free flood-fill check for trapped cavities
boifot_omt/sweep.py           runs, S-matrix, de-embedding, metrics, power conservation
scripts/run_omt.py            CLI sweep driver
scripts/report.py             summary + plots from a finished sweep
scripts/plot_results.py       reflection / transmission / isolation panels
scripts/validate_wr28.py      analytic validation of the port/mesh/PML setup
```

Geometry is built by filling the domain with metal and carving the guide
volumes out of it with higher-priority "air" primitives; the septum (a linear
extruded polygon, so the mitred tip is one primitive) and the pins go back in
as metal at higher priority still.

## Running it

openEMS is not on PyPI and has no Ubuntu package, so it must be built from
source; the Python bindings install into a venv:

```bash
git clone --recursive https://github.com/thliebig/openEMS-Project.git
cd openEMS-Project && ./update_openEMS.sh $PREFIX --python --disable-GUI
```

Then:

```bash
$PREFIX/venv/bin/python3 scripts/validate_wr28.py     # analytic check first
$PREFIX/venv/bin/python3 scripts/run_omt.py --out results/wr28 \
    --mesh-res 0.30 --mesh-res-junction 0.17 --mesh-res-fine 0.09 --nfreq 201
$PREFIX/venv/bin/python3 scripts/report.py --results results/wr28
```

Sweep options: `--excite 0,1`, `--dt-sep-norm`, `--pin-d-norm`, `--n-pins`,
`--no-pins`, `--miter-slope`, `--sep-overhang`, `--trans-gap`, `--n-steps`,
`--ap-len-norm`, `--timesteps`, `--threads`.

Outputs: `sparams.npz`, `sparams.csv`, `boifot_wr28.s5p` (5-port Touchstone),
`metrics.npz`, `params.json`, and light/dark PNGs of the response.

## Validation and acceptance tests

**Straight-guide validation** (`scripts/validate_wr28.py`) runs a uniform
WR-28 guide, where the answer is analytic, through the identical port / mesh /
PML machinery:

```
|S11| worst  : -44.23 dB        |S21| worst : -0.0006 dB
ZL err worst :   1.24 %         phase err   :  1.25 deg over 20.4 mm
```

**Power conservation** is the acceptance test for a junction run. The
structure is lossless, so every excited column of the S-matrix must sum to
1; `scripts/report.py` prints the worst-case error. It catches a truncated
excitation, energy still trapped when the run stopped, a port with the wrong
mode or reference plane, and power leaving through the wrong boundary.

**Trapped-cavity check** (`boifot_omt/check_geometry.py`) rebuilds the air
volumes analytically on the mesh and flood-fills them from the port planes,
with no solver involved. An air volume no port can reach rings as a lossless
cavity and the run never converges.

## Two modelling traps found the hard way

**Sliver cells from colliding refinement.** The FDTD timestep is set by the
smallest cell anywhere in the mesh. Two independently generated refinement
ladders once landed 0.02 mm apart, cutting the timestep ~20×; the excitation
then never completed within the timestep budget and the S-parameters were
garbage, including a nonphysical +2.4 dB return loss. `geometry.mesh_lines()`
now separates hard lines (material boundaries, always kept) from soft ones
(refinement, dropped within 0.75 of a fine cell of a kept line) and raises if
two material boundaries fall closer than the mesh tolerance.

**An under-resolved round post rings.** With the pins at the default mesh the
energy stalled near −33 dB and never decayed. An ablation sweep isolated it to
the pins (−52 dB without them), and refining showed it weakening monotonically
(−33.5 → −39.5 → −45.6 dB), the signature of a staircasing artifact: the
post's discretized surface, tangent to the aperture plane, forms a weakly
coupled sliver resonator. The pins now carry their own refinement of about
four cells across the diameter, held out to 1.5 diameters from the surface.

Note that *physical* trapped-mode resonances are a real feature of this
junction, not only a numerical one — the memo measures them (its Fig. 10) and
ties them to septum-tip deformation and pin asymmetry. A perfectly symmetric
model should not excite them.

## Septum-tip tuning: first sweep

Coarse-mesh sweep (0.40/0.24/0.08 mm, side-arm polarization excitation only,
81 points), reporting the worst in-band side-arm return loss:

| tip_len, tip_w (mm) | worst side RL | worst side IL | power error |
|---|---|---|---|
| 0, 0 (plain 45° wedge) | −1.80 dB | −4.71 dB | 2.6 % |
| 0.6, 0.4 | −2.63 dB | −3.55 dB | 12.4 % |
| 1.2, 0.4 | −2.73 dB | −3.27 dB | 17.8 % |
| 0.6, 0.8 | −3.57 dB | −2.64 dB | 10.2 % |
| 1.2, 0.8 | −3.60 dB | −2.62 dB | 10.3 % |
| 2.0, 0.6 | −3.07 dB | −2.76 dB | 15.8 % |

The trend is monotonic in tip width and roughly saturates in tip length, but
the improvement is small and nowhere near the ~20 dB the memo measures.

**Those power errors were truncation, not physics.** Re-running a tip case at
50 000 timesteps brings power conservation to 0.984–0.997 (1.6 % error), so
the 30 000-timestep sweep above simply had not decayed far enough and its
case-to-case ranking cannot be trusted. Any tuning sweep on this structure
needs to meet the energy criterion before its numbers mean anything.

A `hom_monitor` option adds a TE11 monitor port on the common arm, which
turns power leaving through the overmoded common arm into an accounted
channel rather than a deficit. On a symmetric (common-arm) excitation the
TE11 content is ~0 across the band, confirming the memo's claim that the
junction's twofold symmetry is what keeps TE11/TM11 unexcited. (TM11 is not
covered: openEMS's rectangular waveguide port supports TE modes only.)

There is also a real geometric trade-off between the two ways of adding a
tip, and it is not free either way:

* `tip_auto_slope` holds the aperture at the standard side-arm height and
  lets the miter angle absorb the tip, which makes tip_len and tip_w clean
  one-variable knobs. But blunting the tip then forces a slope well away
  from 45° — and the miter sweep says 45° is strongly preferred. A converged
  run at tip (1.2, 0.8) in this mode is essentially a short circuit
  (worst return loss −0.03 dB).
* `ap_from_tip` keeps the 45° miter and lets the aperture follow the tip,
  which is the better-behaved branch, but the aperture then shifts upstream
  and shortens, so the side-arm guide height moves along with the tip.

### Converged tuning results

Re-run at 50 000 timesteps (mesh 0.40 / 0.24 / 0.10 mm, side-arm polarization
excitation, 45° miter with the aperture following the tip). Only runs whose
power error is small are usable:

| case | worst side RL | power error | usable |
|---|---|---|---|
| 1 pin (baseline) | −1.64 dB | 3.3 % | yes |
| **2 pins** | **−2.90 dB** | 2.3 % | yes |
| 0 pins | +0.69 dB | 17.8 % | no |
| tip (1.2, 0.8), 1 pin | −3.16 dB | 21.6 % | no |

Two pins beat one by 1.3 dB on worst-case side-arm return loss, and both runs
are converged, so that comparison stands.

The two unusable runs are still informative, because *what* they fail at is
consistent: both ring. A positive return loss in dB means |S11| > 1, which is
nonphysical and is simply the signature of a run stopped while energy was
still sloshing around. Removing the pins and blunting the tip each make the
side-arm path ring enough that 50 000 timesteps is not sufficient. Note this
is the opposite of the main-arm behaviour found earlier, where it was the
*presence* of the pins that stalled energy decay — so the pins damp the
side-arm path while loading the main-arm one. That is the same trade the memo
states ("as the total number of pins is increased ... the main-arm return
loss improved while the side-arm return loss gracefully degrades"), showing
up here as a convergence signature rather than as a return-loss number.

The practical lesson for anyone extending this: convergence is
configuration-dependent. Budget timesteps per case and check the power error
before believing any comparison.

### Pin count

Main arm, matched settings, flush pins, all converged (3.7–4.7 % power
error):

| pins | worst RL (in-band) | RL @ 40 GHz | insertion loss |
|---|---|---|---|
| 0 | −12.50 dB | −13.75 dB | −0.37 dB |
| 1 | **−23.40 dB** | −26.30 dB | −0.22 dB |
| 2 | −17.09 dB | −28.49 dB | −0.19 dB |
| 3 | −14.16 dB | −28.91 dB | −0.23 dB |

Two things to read here, and they say different things.

**The memo's trend is reproduced at the top of the band.** Return loss at
40 GHz improves monotonically with pin count, −13.75 → −26.30 → −28.49 →
−28.91 dB, which is the memo's "as the total number of pins is increased
from zero to the three pairs described, the main-arm return loss improved".

**The worst case over the band peaks at one pin**, because added pins also
introduce a degradation away from the top of the band. Collapsing the band
to a single worst-case number therefore shows the reverse trend, which is
how an earlier revision of this file came to claim the model disagreed with
the memo. It does not; the metric did.

The zero-pin row is the direct confirmation of what the pins are for: with
no pins the main arm is at −12.5 dB with 0.37 dB of insertion loss, because
the side-arm apertures mismatch it. That is the memo's "the pins act as
symmetric short-circuit waveguide-stubs which tune out the discontinuity due
to the side-arm ports".

Side-arm return loss moves the other way, −1.63 dB at 1 pin to −2.90 dB at
2 pins, so the pin count is a genuine two-sided trade.

**Pin placement is first-order, not a detail.** The memo is explicit for one
pair — "the outer diameter of the pin flush to the side-arm wall (i.e. the
pins are physically in the common-arm) and centered in the side-arm
aperture" — but ambiguous for two, where "centered in the aperture" cannot
put both pins at the centre. Reading it as centred in *x*, straddling the
aperture plane, instead of flush behind it, is worth 4.6 dB of main-arm
return loss at two pins (−21.69 vs −17.09 dB). `pin_inset` covers both: 0 is
flush, −d_pin/2 straddles. Any multi-pin tuning has to state which
convention it uses.

For a single figure of merit across the whole band, one flush pin is the
best of those tried. If the top of the band matters most, more pins are
better.

### Pin geometry is not the lever for the side-arm match

Six pin variants, side-arm polarisation, 90 000 timesteps, S-parameters
saved. Return loss in dB, sampled at the memo's own Figure 2 x-axis so the
curve *shapes* can be compared rather than a worst-case scalar:

| case | 1.30 | 1.45 | 1.60 | 1.70 | 1.80 | 1.87 |
|---|---|---|---|---|---|---|
| full post, centre, flush | −21.3 | −19.6 | −14.7 | −10.7 | −6.1 | −2.9 |
| diameter 0.07 a₀ | −21.0 | −19.6 | −14.5 | −10.4 | −5.8 | −2.7 |
| moved upstream (z 0.25) | −18.6 | −17.9 | −13.9 | −10.1 | −5.3 | −1.3 |
| moved downstream (z 0.75) | −23.2 | −21.1 | −15.7 | −11.3 | −6.0 | −2.3 |
| half-height, on the septum | −20.9 | −16.5 | −11.1 | −4.8 | −0.4 | −1.3 |
| straddling the aperture plane | −21.4 | −19.9 | −14.8 | −10.6 | −5.5 | −2.3 |
| *memo Fig 2, curve 1* | *≈−22* | *≈−28* | *≈−20* | *≈−18* | *≈−22* | *≈−26* |

Every variant traces the same monotonic collapse, within about 2 dB of the
others, against a memo curve that stays flat. Diameter across the whole
allowed range is worth 0.2 dB.

The reason is orientation, and it is not a tuning matter. A post along y is
perpendicular to the side-arm polarisation's E-field in *both* arms — E ∥ x
in the common arm, E ∥ z in the side arm — so it barely couples to it. The
same post is parallel to the main-arm E-field, which is why adding the first
pin moves the main arm 10.9 dB. The pin is electrically strong; it is strong
for the other polarisation.

### The gap is in the junction, not the compensation

The comparison that matters is against Figure 2 curve 0, the memo's junction
with **no pins at all**: it sits near −30 dB at 1.6 f_c. This model, *with* a
pin, is −14.7 dB there. Their uncompensated junction beats this compensated
one through mid-band, so the shortfall was never primarily about pins.

The admittance extraction agrees. Looking into the junction from the common
port, de-embedded:

| f/f_c | \|S11\| | G | B |
|---|---|---|---|
| 1.26 | −20.8 dB | 1.09 | +0.17 |
| 1.57 | −16.3 dB | 1.04 | +0.32 |
| 1.72 | −10.2 dB | 1.34 | +0.67 |
| 1.82 | −5.4 dB | 3.04 | +0.87 |

Mid-band the conductance is essentially unity and the mismatch is a
susceptance of ~0.3 that should not be there in a well-formed junction. Only
near the top of the band does G run away, and no susceptance — pin, tip or
otherwise — can pull G back to 1.

Next diagnostic is the unambiguous one: run this junction with no pins and
compare directly against curve 0, which has no placement ambiguity to argue
about. If it still falls short mid-band, the discrepancy is in the junction
proper (septum, miter, aperture), not in any compensating element.

### Septum thickness

All at 2 pins, 70 000 timesteps, converged (1–4 % power error):

| dt/a₀ | dt | worst side RL |
|---|---|---|
| 0.020 | 0.142 mm | −2.47 dB |
| 0.035 | 0.249 mm | −2.90 dB |
| 0.050 | 0.356 mm | −3.30 dB |

Monotonic, and the whole span is 0.8 dB. In *this* structure a thicker
septum is slightly better, which runs opposite to the memo's remark that
"the thinner the septa, the easier the structure is to compensate" — but
that remark is about ease of compensation in a tuned design, and this one is
under-tuned, so the two are not in direct conflict. Worth revisiting once a
proper compensation exists.

### Transformer section lengths

All four sections are one quarter guide wavelength, 2.954 mm. That is correct
for this transformer rather than a simplification: the steps change only the
height, which is the E-plane dimension, so the cutoff stays pinned by a₀ and
λg is identical in every section. An H-plane (width) transformer would need a
different length per section.

`first_step_len` exists to override the first section alone, because ideal
small-reflection theory treats each junction as a perfect impedance step
while each real step carries a shunt susceptance — and the first section is
the odd one out, since it faces the junction (aperture, septum tip, pins)
rather than a plain guide. Swept at 2 pins:

| first section | worst side RL | worst main RL | main RL @ 40 GHz |
|---|---|---|---|
| 1.5 mm | −2.50 dB (18 % err, unusable) | — | — |
| **2.954 mm (λg/4)** | **−2.90 dB** | **−17.09 dB** | **−28.49 dB** |
| 4.5 mm | −2.74 dB | −16.47 dB | −18.91 dB |
| 6.0 mm | −2.67 dB | — | — |

**Keep them equal.** Stretching the first section gains nothing on the side
arm and costs the main arm 9.6 dB at the top of the band, where breaking the
quarter-wave synthesis hurts most.

The side-arm result is the more interesting one, because it is a clean
negative. 1.5 → 6.0 mm walks the main-arm termination through more than half
a guide wavelength of phase, and the side-arm worst case moves by 0.2 dB.
So the side-arm match is *not* sensitive to where the transformer sits, which
corrects an earlier reading of the memo's "dispersive reflective termination"
as a movable short whose phase could be tuned. It cannot be: the side-arm
polarisation is evanescent in the half-channels downstream of the septum with
a 1/e depth of 1.4–2.9 mm, so by the first step boundary it has already
decayed one to two e-foldings and the steps are screened. What the side arm
sees is the cutoff region at the septum itself, and the dispersion is in the
evanescent decay, not in a repositionable short.

That rules out first-step length, `trans_gap` and transformer placement as
routes to the side-arm match. Whatever compensates it has to sit where the
field is still strong — inside the aperture, which leaves the septum tip and
the pins.

### Where the side-arm match actually stands

Every knob tried so far moves worst-case side-arm return loss by about a dB,
and they look roughly independent:

| knob | range tried | effect |
|---|---|---|
| pin count | 1 → 2 | 1.3 dB |
| septum thickness | 0.02 → 0.05 a₀ | 0.8 dB |
| miter slope | 0.6 / 1.0 / 1.25 | 45° strongly preferred |
| septum tip | blunt finger | no usable gain; makes the structure ring |

Best converged configuration so far is 2 pins with dt/a₀ = 0.05, at
−3.30 dB. The memo measures roughly −20 dB across the band, so about 17 dB
is unaccounted for, and crucially **the shape has not changed** under any of
these: matched at the band bottom, rolling off monotonically toward 1.9 f_c.
That is the signature of an uncompensated dispersive termination, not of a
mistuned resonator.

Closing it is unlikely to come from more one-at-a-time sweeps. Two routes
look better:

1. **Extract the equivalent reactance.** The junction is well matched at the
   band bottom and degrades smoothly, so fit the junction susceptance versus
   frequency from the simulated S11 and design the pin/tip compensation
   against that, using FDTD only to verify. Far cheaper than a blind search.
2. **Use a frequency-domain solver for the tuning loop.** The memo optimized
   this in HFSS. Each converged FDTD run here costs 10–20 minutes, and a
   genuine multi-parameter search over tip length, tip width, pin count, pin
   diameter and pin position needs hundreds of them. FEM or mode-matching on
   a junction this small is seconds per point; FDTD is the wrong tool for
   the optimization inner loop, even though it is a good one for the final
   verification.

Separately, the miter angle was swept at 45° ± : slope 0.60 and 1.25 are both
far worse than 1.00, supporting the 45° miter as the right choice.
