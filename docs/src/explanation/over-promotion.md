# Over-promotion

A symmetry-finding code can fail in two directions: it can report a symmetry the structure does not have (**over-promotion**), or miss a symmetry the structure does have (**under-promotion**). The two failures look very different to a downstream user, and the asymmetry between them shapes a lot of Spacey's design.

This page explains why over-promotion is the dangerous one, how it shows up in Spacey, and what `verify_stable` does to detect it. The user-facing recipe is in [Detect tolerance-dependent answers](../how-to/detect-tolerance-dependence.md); this is the *why* behind the design.

## The two failure modes are not symmetric

Both directions are wrong, but they're wrong in different ways:

| | Over-promotion | Under-promotion |
|---|---|---|
| Symptom | More ops than the structure actually has | Fewer ops than the structure actually has |
| Visibility | **Silent** — no warning; user gets the answer they were hoping for | **Loud** — visibly fewer ops than expected; immediately suspicious |
| Recovery | Hard — user has no signal that anything is wrong; bad results propagate | Easy — user notices, loosens tolerance, gets the right answer |
| Downstream impact | Wrong space group → wrong DOS, wrong phonon dispersion, wrong selection rules | Usually re-runs before publication |

This asymmetry is what makes over-promotion the first-order bug to design against. Spacey's defaults err toward under-promotion for that reason (see [Tolerances](tolerances.md)), but for noisy input the right defaults aren't always tight enough.

## How over-promotion happens in Spacey

The mechanism has three steps, all flowing from the way the algorithm filters candidate symmetry operations (see [Algorithm overview](algorithm-overview.md) §3):

1. **Spurious operations pass the integer-matrix test.** For a near-symmetric input, the would-be-cubic operations on a tetragonal-near-cubic lattice produce small-but-nonzero residuals on the integer-matrix test. At a tolerance loose enough to absorb input noise, those residuals fall below `tol` and the operations pass.
2. **The spurious operations close into a higher-symmetry group.** The group-closure step reports the *largest* subset of accepted operations that closes under composition. If the spurious operations from step 1 are mutually consistent — which they are, by design, since they came from a near-cubic lattice — they close into the cubic point group.
3. **The user gets the higher-symmetry answer with no warning.** The standard return path of `pointGroup` and `spacegroup` is silent on this kind of disagreement.

The canonical illustration is a tetragonal lattice with `c/a = 1 + ε`. The 32 cubic-only operations (the ones that mix the c-axis with the a- or b-axes — for example, the 120° rotation about the body diagonal) appear in lattice coordinates with entries `a/c = 1 / (1 + ε)` and `c/a = 1 + ε`. Where the cubic case requires exactly `±1`, the tetragonal case has entries that differ from `±1` by `O(ε)`. At any `tol ≥ ε`, the integer test passes, and the 3-fold rotation joins the group. Eight of these spurious 3-fold operations close with the 16 tetragonal operations into the full 48-op cubic group. Result: 48 reported, 16 correct.

This isn't subtle: at `ε = 10⁻³` and the default `tol = 10⁻²`, this failure happens silently — no warning, no anomaly. The user gets "cubic" as their answer. The same mechanism applies to space groups: at `pos_tol > displacement`, a slightly-displaced atom matches its parent-prototype image and the displaced-symmetry op is gained.

## The (tol, ε) crossover

The shape of the failure region is geometric. Let `ε` be the structural distortion (how far the structure is from a higher-symmetry parent) and `tol` be the tolerance. Three regimes:

- **`tol << ε`**: tolerance is tighter than the distortion. Spurious operations fail the integer-matrix test. Correct answer (lower symmetry) reported.
- **`tol >> ε`**: tolerance loose enough that spurious operations pass. Higher-symmetry parent reported (over-promotion).
- **`tol ≈ ε`**: the boundary. Tolerance comparable to the distortion. The reported group depends sensitively on which operations happen to land just below or just above `tol`.

The boundary is where the answer is *fragile*. A 10× change in `tol` flips the reported group. The user has no signal that they're at the boundary — the same `pointGroup(A; tol=...)` call returns a different number depending on `tol`. The diagnostic at `test/nearMissBoundary.jl` plots a `(tol, ε)` heatmap of group-size-reported, with a clear diagonal failure region exactly along `tol ≈ ε`.

## What `verify_stable` does

`verify_stable=true` re-runs the algorithm at a tolerance 1000× tighter than the user's choice and emits a `@warn` if the reported group order changes. The intuition: a *real* symmetry should still be a symmetry at three orders of magnitude tighter tolerance. If it isn't, the answer was tolerance-dependent and the user is probably in the boundary regime.

The 1000× factor is calibrated to catch genuine over-promotion without false-alarming on numerical noise. `tol = 0.01` becomes `tol = 1e-5` for the verification run, which is safely above floating-point noise (1e-15) but well below typical structural distortions (1e-2 to 1e-4). On a clean cubic lattice, both runs return 48; no warning. On a tetragonal-near-cubic case where `tol = 0.01` reports 48 ops but `tol = 1e-5` reports the correct 16, the warning fires.

The cost is one extra invocation per call. For one-off symmetry queries this is microseconds. For corpus-scale work (autoGR-class uses, scanning thousands of grids) the cost matters more, and the user can disable the check after the corpus has been verified clean.

The same machinery applies to `spacegroup`'s `pos_tol`: re-run at `pos_tol/1000`, warn if the operation count changes. The canonical case is a ferroelectric like BaTiO₃ where the Ti displacement (~0.05 Å) sits below a too-loose `pos_tol` and the displaced structure is reported with the cubic parent symmetry — see the worked example in [Detect tolerance-dependent answers](../how-to/detect-tolerance-dependence.md).

## How other codes handle over-promotion

The major comparator codes take different approaches to the same problem:

- **spglib** uses *iterative tolerance reduction*: if the initial answer fails an internal consistency check (group closure, Wyckoff consistency), `symprec` is tightened and the search is re-run. The key design choice is that tolerance is "generally only decreased, not increased" — spglib will never report a higher symmetry than its tightest convergent run. This is similar in spirit to `verify_stable`, but the *reported* answer is the iterated-down one, not the user's original. The user doesn't know what the original tolerance would have returned.
- **AFLOW-SYM** uses an iterative consistency check on Niggli-cell parameters and atom positions, with similar "tighten until consistent" semantics. Boundary cases between Bravais types can occasionally produce convergence to a subgroup.
- **FINDSYM** tests against the catalog of 230 known space groups. If a candidate group fails a position consistency check at the user's tolerance, the next-largest subgroup is tried. The reported answer is the largest catalog group consistent with the input — implicit over-promotion-resistance, but no boundary detection.

`verify_stable` is structurally different: Spacey returns the answer at the user's tolerance and *separately* reports whether tightening would change it. The user retains the choice — accept the higher-symmetry answer (which may be the right one for noisy data with a real parent symmetry), tighten and re-investigate, or reject the over-promotion. This is information that spglib's iterative-tightening approach makes inaccessible (you don't know what was originally accepted before being rejected).

## What to do when you can't tighten

Sometimes the boundary is a *physical* boundary, not a numerical one. A real ferroelectric just above the Curie point genuinely has a small displacement that is the physics. A noisy experimental structure with displacements smaller than the experimental ESD genuinely cannot be resolved with the given data.

The honest behavior in those cases is not to commit to one side. Spacey's design choice — return the answer the user asked for, plus an optional warning — preserves that honesty. It does not replace careful judgment about whether the input data is even capable of resolving the question.

[Snap-to-symmetry](../how-to/snap-to-symmetry.md) is the workflow for the case where the user *intends* to accept the higher-symmetry parent and wants the lattice / atom positions cleaned up to match. It's a deliberate commitment to the parent symmetry, not a way to recover the truly-distorted structure.

## See also

- Explanation: [Tolerances](tolerances.md), [Algorithm overview](algorithm-overview.md), [Why Minkowski reduction](why-minkowski.md)
- How-to: [Detect tolerance-dependent answers](../how-to/detect-tolerance-dependence.md), [Handle noisy real-world data](../how-to/handle-noisy-data.md), [Snap a noisy lattice to symmetry](../how-to/snap-to-symmetry.md)
- Reference: [`pointGroup`](../reference/point-groups.md), [`spacegroup`](../reference/space-groups.md)
