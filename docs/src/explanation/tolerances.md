# Tolerances

Spacey's symmetry decisions are not floating-point exact. They depend on two tolerances — one for the lattice geometry, one for atomic positions — and the right way to set them depends on the noise in your input. This page explains what each tolerance controls, how the defaults were chosen, and how Spacey's two-tolerance design contrasts with the single-tolerance design used by most other codes.

For the practical recipe ("I have data with noise level X, what tolerances should I use?"), see [Handle noisy real-world data](../how-to/handle-noisy-data.md). This page is the *why* behind those defaults.

## Two scales, two physical sources of noise

The lattice and the atomic positions accumulate noise from independent sources:

- **Lattice noise** comes from cell-parameter refinement (X-ray Pawley fit), DFT cell relaxation, or coordinate-frame rotation — typically smooth and small (4–5 significant figures for a relaxed structure).
- **Position noise** comes from atomic-position refinement, thermal motion in single-shot snapshots, or DFT relaxation thresholds — typically larger than lattice noise (3–4 significant figures).

A symmetry analysis that uses one tolerance for both can either accept too much (loose) or reject too much (tight) on at least one of these scales. Spacey's two-tolerance design lets each scale match the noise it actually has.

| Tolerance | Used by | Operates on | Default |
|---|---|---|---|
| `tol` (or `lattice_tol` for `spacegroup`) | `pointGroup`, `spacegroup` lattice step | volume-normalized lattice vectors | `0.01` |
| `pos_tol` | `spacegroup`, `isSpacegroupOp` | atomic positions in the user's lattice basis | `default_pos_tol(c) = 0.01 · (V/N)^(1/3)` |

## Lattice tolerance: `tol`

`tol` controls the **integer-matrix test** in the point-group pipeline (see [Algorithm overview](algorithm-overview.md) §3). For a candidate basis `B`, the matrix `U = inv(A) · B` must be close to integer, with each entry within `tol` of the nearest integer.

`tol` is **relative** because Spacey volume-normalizes the input first: every lattice is rescaled to unit volume before the comparison. This means `tol = 0.01` carries the same meaning whether your cell is in Ångström, Bohr, or femtometres — it's "1% of the natural length scale of the cell." The user doesn't have to remember which units they used.

**Why 0.01 and not something tighter or looser?** The default sits at the level of typical experimental cell-parameter ESDs (10⁻³–10⁻²) and tight DFT-relaxation residuals (10⁻⁴–10⁻³), comfortably above floating-point noise (10⁻¹⁵) and below the structural distortions that should resolve as different symmetry. For a tetragonal cell with c/a = 1 + 10⁻², the cubic point group passes the integer test (a near-miss); for c/a = 1 + 10⁻³, it doesn't. This is intentional — the structural distortion that should *not* read as cubic is the one that's larger than `tol` in magnitude.

## Position tolerance: `pos_tol`

`pos_tol` controls atom matching under candidate space-group operations: an op `(R, τ)` is accepted iff every atom maps to a same-type atom of distance < `pos_tol` modulo the lattice. Unlike `tol`, `pos_tol` is **absolute** — measured in the same units as the lattice matrix.

The default formula is

```math
\text{pos\_tol} = \alpha \cdot (V / N)^{1/3}, \quad \alpha = 0.01
```

where `V` is the cell volume and `N` is the atom count. The factor `(V/N)^{1/3}` is the *characteristic atom separation* in the cell — for a simple cubic atomic arrangement with one atom per cell, it's the lattice parameter exactly. The 1% factor (`α = 0.01`) means "atoms that move by 1% of typical bond length are still considered to be at their reference position."

The formula scales correctly under any unit choice because both `V` (length³) and `(V/N)^{1/3}` (length) rescale together. A user with a lattice in Bohr and a user with the same lattice in Ångström both get a `pos_tol` that's 1% of the characteristic bond length in their units.

## Why α = 0.01: the asymmetry argument

Choosing α was not a single trade-off but two competing failure modes. The discussion is preserved in detail in `designDiscussions.md`; the key argument:

- **Over-promotion** (loose `pos_tol`, displaced atom maps onto a parent prototype position): user reports a *higher* symmetry than the structure has. Silent — no warning, no anomalous output. Wrong space group propagates into every downstream analysis (DOS, phonons, Raman/IR selection rules).
- **Under-promotion** (tight `pos_tol`, atom that should match its image fails): user reports a *lower* symmetry than the structure has. Loud — fewer operations than expected, immediately noticeable, recoverable in one keyword override.

These two failure modes are not symmetric. The default should err toward the loud, recoverable case rather than the silent, wrong one. Concretely, at α = 0.05, `pos_tol ≈ 0.1 Å` for a typical 4 Å cell — which silently classifies BaTiO₃ at room temperature (Ti displaced 0.05 Å from cubic) as cubic Pm3̄m instead of tetragonal P4mm. At α = 0.01, `pos_tol ≈ 0.02 Å`, and the same Ti displacement registers as the symmetry-breaking it actually is.

The classes of structure that this saves:

| Class | Displacement scale | Mis-classified at α = 0.05? |
|---|---|---|
| Ferroelectrics (BaTiO₃, PbTiO₃) | 0.05–0.10 Å | yes |
| Perovskite octahedral tilts (CaTiO₃, SrTiO₃ < 105 K) | 0.02–0.08 Å | often |
| Jahn–Teller distortions (LaMnO₃, KCuF₃) | 0.05–0.20 Å | yes |
| Peierls dimerizations (1D CDW) | 0.03–0.08 Å | yes |
| Martensite distortions (NiTi, β-brass) | < 0.10 Å | yes |

The one-atom test cube with V = N = 1 confirms the formula's calibration: `default_pos_tol = 0.01` exactly.

### Where α = 0.01 is too loose

The default is not infallible — it's just less wrong than the alternatives. Specifically, α = 0.01 still over-promotes:

- Quantum paraelectrics (SrTiO₃ below ~4 K — sub-0.01 Å displacements).
- Relaxor ferroelectrics (Pb(Mg,Nb)O₃, Pb(Sc,Ta)O₃) with sub-0.01 Å off-centring.
- Tight-convergence DFT data where meaningful displacements can be 0.005 Å.
- Small-amplitude phases near a continuous transition.

For these, the user must override with `pos_tol = 1e-3` or tighter. [Detect tolerance-dependent answers](../how-to/detect-tolerance-dependence.md) describes the `verify_stable` machinery that catches the case where the chosen `pos_tol` happens to land on the wrong side of a boundary.

### Where α = 0.01 is too tight

The default *under*-promotes — loudly, fixable — for:

- Powder X-ray (typical ESDs 0.01–0.05 Å)
- Hydrogen-atom positions from X-ray (0.05–0.1 Å uncertainty)
- High-temperature snapshots with large thermal motion
- Early-iteration DFT relaxations

The user override `pos_tol = 0.05` or similar resolves these.

## How Spacey's two-tolerance design compares with other codes

Most other space-group finders use a single tolerance:

- **spglib** uses `symprec` — an absolute Euclidean distance (default `1×10⁻⁵ Å`). It's applied uniformly to atom matching during the symmetry search; lattice noise enters the same threshold without separate accounting. Spglib's iterative tightening adjusts `symprec` per call but doesn't separate the two scales.
- **AFLOW-SYM** uses an iterative tolerance on Niggli-cell parameters and atom positions together, again without a hard distinction between the two scales.
- **FINDSYM** uses one tolerance applied to position-equivalence tests.

A single tolerance that's tight enough for the lattice can be too tight for atoms (and vice versa). The two-tolerance design pays a small documentation cost (the user has to learn there are two knobs) for clearer per-scale control and more honest defaults.

Spacey's `verify_stable` flag — re-run at 1/1000 the input tolerance, warn on disagreement — is also a deliberate alternative to spglib's iterative tightening. Iterative tightening *changes* the answer until it converges; `verify_stable` returns the user's requested answer and reports whether tightening would change it. The two strategies are complementary; they both address noise but report it differently. See [Over-promotion](over-promotion.md) for the full discussion.

## See also

- Explanation: [Algorithm overview](algorithm-overview.md), [Over-promotion](over-promotion.md), [Why Minkowski reduction](why-minkowski.md)
- How-to: [Handle noisy real-world data](../how-to/handle-noisy-data.md), [Detect tolerance-dependent answers](../how-to/detect-tolerance-dependence.md)
- Reference: [`pointGroup`](../reference/point-groups.md), [`spacegroup`](../reference/space-groups.md), [`default_pos_tol`](../reference/crystals.md)
