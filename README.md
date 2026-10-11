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
