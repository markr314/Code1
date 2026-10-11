# 4-port Bøifot orthomode junction in openEMS (WR-28, 26.5–40 GHz)

A parameterized openEMS (FDTD) model of a 4-port Bøifot orthomode junction and
the harness that extracts its modal S-parameters.

**Status: the simulation framework is built and validated; the junction
geometry is not yet correct.** See [Current state](#current-state) below —
the published geometry from Wollack, Grammer & Kingsley is needed to finish it.

## Port definition

A Bøifot junction is natively a 4-port device whose common port is dual-mode,
so the S-matrix is 5x5 in modal terms:

| idx | port | mode |
|----|------|------|
| 0 | common square guide, polH | E along x |
| 1 | common square guide, polV | E along y |
| 2 | main arm | — |
| 3 | side arm, +x | — |
| 4 | side arm, −x | — |

The two side arms are brought out separately; combining them (E-plane tee or
magic-T) is left outside this model.

## Band and sizing (WR-28)

| quantity | value |
|---|---|
| band | 26.5–40 GHz |
| common square guide | 7.112 x 7.112 mm |
| common guide TE10/TE01 cutoff | 21.08 GHz |
| common guide TE11/TM11 cutoff | 29.81 GHz (overmoded in band; held off by symmetry) |
| main arm | 3.556 (x) x 7.112 (y) mm |
| main arm operating cutoff | 21.08 GHz |
| main arm rejected-pol cutoff | 42.15 GHz |
| side arms | 3.556 (y) x 7.112 (z) mm |

A square common guide of side 7.112 mm steps directly into standard WR-28 and
is the natural choice for this band, but note it is overmoded above
29.81 GHz. The junction's two-fold symmetry is what keeps TE11/TM11
unexcited, so symmetry errors (machining asymmetry, a mis-centred septum)
show up as cross-polarisation.

## Layout

```
boifot_omt/waveguide.py   WR table, cutoff / beta / wave impedance / wall loss
boifot_omt/params.py      BoifotParams: every dimension, plus electrical checks
boifot_omt/geometry.py    CSXCAD primitives and the critical mesh lines
boifot_omt/model.py       FDTD setup, boundaries, the five modal ports
boifot_omt/sweep.py       per-excitation runs, S-matrix, de-embedding, metrics
scripts/run_omt.py        CLI sweep driver
scripts/validate_wr28.py  analytic validation of the port/mesh/PML setup
```

Geometry is built by filling the domain with metal and carving the guide
volumes out of it with higher-priority "air" primitives; the septum and pins
go back in as metal at higher priority still. The junction is then a handful
of boxes rather than an assembly of walls.

## Running it

openEMS is not on PyPI and has no Ubuntu package; it has to be built from
source (`openEMS-Project`, `./update_openEMS.sh <prefix> --python
--disable-GUI`), which installs the Python bindings into a venv under
`<prefix>/venv`. Use that interpreter:

```bash
$PREFIX/venv/bin/python3 scripts/validate_wr28.py            # analytic check
$PREFIX/venv/bin/python3 scripts/run_omt.py --out results/wr28
```

Useful options: `--excite 0,1` (subset of excitations), `--mesh-res`,
`--mesh-res-fine`, `--t-sep`, `--l-sep`, `--no-pins`, `--d-pin`, `--pin-z`,
`--n-steps`, `--step-len`, `--nfreq`, `--threads`.

Outputs: `sparams.npz`, `sparams.csv`, `boifot_wr28.s5p` (Touchstone),
`metrics.npz`, `params.json`.

At the default mesh (0.25 mm, 0.09 mm locally) the model is ~1.2 M cells and
each excitation runs in a few minutes on 4 cores.

### Validation

`scripts/validate_wr28.py` runs a straight WR-28 guide, where the answer is
analytic, through the identical port/mesh/PML machinery:

```
|S11| worst     : -44.23 dB
|S21| worst     : -0.0006 dB
ZL err worst    : 1.24 %
S21 phase err   : 1.25 deg over 20.4 mm
```

So the ports, PML and mesh are sound; anything the junction model shows is
the geometry, not the solver.

### A mesh trap worth knowing about

The FDTD timestep is set by the smallest cell anywhere in the mesh. Two
independently generated refinement ladders landed 0.02 mm apart in the first
version, which cut the timestep ~20x; the excitation then never completed
inside the timestep budget and the resulting S-parameters were garbage
(including a nonphysical +2.4 dB return loss). `geometry.mesh_lines()` now
separates hard lines (material boundaries, always kept) from soft ones (local
refinement, dropped when they fall within half a fine cell of a kept line),
and raises if two material boundaries are closer than the mesh tolerance.

## Current state

The geometry currently implemented is a physics-derived interpretation of the
Bøifot arrangement, **not** the published one:

- septum in the x=0 plane, transparent to polH (E normal to it), shorting
  polV (E tangential and in-plane);
- main arm narrowing in x, which leaves polH's cutoff untouched and puts
  polV below cutoff;
- side arms branching along ±x, coupling through apertures in the x=±a_sq/2
  walls;
- pins with axes along y standing in each aperture plane.

A coarse sweep (0.5 mm mesh, excitations 0 and 1) shows this does not work:

| quantity | result |
|---|---|
| polV return loss | ~0 dB (totally reflected) |
| polV into the two side arms | −10 to −25 dB (should be ~−3 dB each) |
| polV → main arm isolation | better than −84 dB |
| side-arm amplitude imbalance | < 1e-5 dB |
| side-arm relative phase | 0 deg |

The isolation and symmetry figures confirm the septum and the model's
symmetry work as intended. The coupling figure identifies the flaw: with a
septum at x=0, polV is simply *reflected* at the septum's leading edge rather
than guided into the side arms. Behind the septum the common guide is split
into 3.466 mm half-channels in which polV is evanescent, with a 1/e decay of
1.4–2.9 mm across the band, while the apertures run 7.112 mm along z. Most of
the aperture therefore sits where polV's field has already died away, and the
apertures are in the walls where polV's E is tangential and near zero to
begin with.

Resolving this needs the actual published arrangement — specifically which
wall the side-arm apertures are in relative to the septum plane, where the
septum's leading edge sits relative to the apertures, and the pin geometry.
Each polarisation couples strongly only through the walls where its E-field
is *normal*, and the septum shorts whichever polarisation lies in its plane;
those two constraints pull against each other, and the published geometry is
what reconciles them.

Needed from *The Bøifot Orthomode Junction* (Wollack, Grammer & Kingsley,
ALMA Memo 425): the junction figure (septum orientation and extent,
aperture position/size, pin count/diameter/location) and the normalized
dimension table, which `BoifotParams` is already shaped to accept.
