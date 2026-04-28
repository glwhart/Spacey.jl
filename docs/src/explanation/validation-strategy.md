# Validation strategy

A symmetry-finder is hard to test against a single ground-truth source — there isn't one. The International Tables list space groups by abstract type (Pm3̄m, Fd3̄m, …) but not the specific operations any concrete representative carries. Crystallographic databases like AFLOW, ICSD, COD list space-group *labels* but not always the full operation set. Other codes (spglib, FINDSYM, AFLOW-SYM) are themselves tested against each other, with documented disagreements at boundary cases.

Spacey's response is to validate against multiple independent invariants and pin disagreements explicitly. This page explains the strategy and what each layer catches.

## Three algorithm variants for cross-validation

Three internal point-group finders coexist in the package, deliberately:

| Variant | What it does | Used as |
|---|---|---|
| `Spacey.pointGroup_simple` | Brute-force: iterates over every 3×3 matrix in `{-1, 0, 1}⁹` (= 19,683 candidates), keeps those that preserve the metric tensor with default `isapprox` | Reference truth — slow but obviously correct |
| `Spacey.pointGroup_fast` | Optimized filtering pipeline, strict `isapprox` tolerance, no `tol` knob | Production speed for clean inputs |
| `Spacey.pointGroup_robust` | Same pipeline as `_fast`, but with `tol`-controlled integer-matrix test and group-closure step | The version `pointGroup` delegates to (real-world inputs) |

The test suite cross-checks all three on every clean (non-noisy) Bravais lattice — they must agree on the operation count for each of the 14 Bravais lattice prototypes. If any of the three drifts away from the others, that's an immediate signal: either a bug in the variant, or a regression in shared code.

This is `pointGroup_simple`'s only purpose. It's slow (~100× slower than `_fast`), but its correctness is *transparent* — there's no tolerance, no integer test, no group closure — and it serves as the ground truth that the optimized paths must match.

## The 14-Bravais lattices hand-built panel

For each of the 14 Bravais lattices, the test suite has a hand-curated representative basis with the expected operation count baked in (cubic = 48, hexagonal = 24, …). Each is exercised against:

- All three `pointGroup` variants (must agree).
- Rotation invariance: applying an arbitrary 3D rotation to the basis must not change the operation count.
- Small-noise tolerance: ~10⁻⁸ random perturbation must not crash or change the result.
- High aspect ratio: `AR = 256, 500, 512, 1024` are explicitly tested.
- Snap-to-symmetry round-trip: noisy lattice → `pointGroup` → `snapToSymmetry_SVD` → snapped basis → `pointGroup` again gives the same op count, with volume preserved.

This panel is the first line of defense. A regression in any algorithm component will show up here before it hits the larger corpus tests below.

## The AFLOW prototype corpus

The AFLOW Library of Crystallographic Prototypes ([Mehl 2017, Hicks 2019, Hicks 2021](https://aflow.org/CrystalDatabase/)) catalogs 1095 distinct prototype structures from three published papers — every meaningful crystal structure observed in inorganic chemistry, with its lattice and atomic positions in the conventional setting and a label encoding the space-group number plus the Pearson symbol.

The label gives Spacey two independent invariants per structure:

- **Operation count**: the third underscore-separated segment is the space-group number (1–230); the operation count for that group is a known integer.
- **Crystal system**: the first one or two letters of the Pearson symbol encode the system (a/m/o/t/c for triclinic/monoclinic/orthorhombic/tetragonal/cubic, hP for hexagonal, hR for trigonal).

Spacey's test suite auto-generates a test crystal from each prototype (via `tools/generate_aflow_tests.jl`) and checks both invariants. **Different bugs surface in different invariants** — some structures pass op-count but fail crystal-system (lattice over-promotion that happens to round to a "right-looking" op count); others pass crystal-system but fail op-count (atomic-position mis-identification on a correctly-labeled lattice). Running both invariants catches a strictly larger deviation set than either alone. See [Crystal system vs full Bravais](crystal-system-vs-bravais.md) for the design decision behind shipping the second invariant cheaply.

The full corpus runs in CI on every push: 286 + 299 + 510 = 1095 structures × 2 invariants = 2,190 individual tests, finishing in a few seconds.

## `@test_broken` for known deviations

Of the 1095 prototypes, a small minority don't match Spacey's output for various reasons — over-promotion at the default tolerance, atomic positions at coincidental high-symmetry locations, structures with partial occupations that the test generator handles approximately. Rather than hide these failures or change the algorithm to match the labels, the test suite uses Julia's `@test_broken`: a marker that pins the *current* deviating behavior so any future regression — in either direction — is loud.

A `@test_broken` reading "expected 192, got 96" tells the maintainer:

- The current algorithm gives 96 for this structure (and does so consistently).
- The "right answer" per the AFLOW label is 192.
- A future change that makes Spacey return 192 here will *fail* the test (because `@test_broken` flips to passing, and the test author has to consciously change it from `@test_broken` to `@test`).
- A future change that makes Spacey return something other than 96 *also* fails — the test breaks the wrong way and the change has to be investigated.

This is a deliberate design choice over either dropping the failing tests (loses the deviation as a known issue) or changing the algorithm (might mask other regressions). The set of `@test_broken` cases is the project's documented deviation catalog: as of v0.7.2, ~5% of corpus prototypes are broken-pinned, with the failure mode classified inline in test comments.

## The `verify_stable` diagnostic tests

Two heatmap diagnostics live in the test directory but aren't run as part of normal CI:

- `test/nearMissBoundary.jl` plots a `(ε, tol)` heatmap of `pointGroup` operation counts on a tetragonal-near-cubic lattice, marking each cell with whether `verify_stable` would catch the over-promotion.
- `test/nearMissBoundaryCrystal.jl` plots the same for `spacegroup` on a BaTiO₃-style ferroelectric, varying the Ti displacement and `pos_tol`.

These are pinned by *unit tests* that exercise specific cells of the heatmap (a known over-promotion case fires the warning; a known stable case stays silent). The diagnostic itself produces a 2D table of marked cells when run interactively — the diagonal failure region (`tol ≈ ε`) is visually striking and, in the short term, the most direct evidence that `verify_stable` is doing what it claims.

## Pre-merge gating

The Documentation workflow doesn't run the AFLOW corpus tests on every PR by default — it would dominate CI time. The `Runtests.yml` workflow does run them on every push to main, on Ubuntu / macOS / Windows × Julia 1.11 × x64. For the documentation-only changes that make up most of the recent commits, the test suite still runs in `<30 seconds`.

The full test count as of v0.7.2: ~3,400 tests passing, with the broken-pinned set documented in test comments.

## What this catches and what it doesn't

The validation strategy is good at:

- Algorithmic regressions (a change that breaks a working code path)
- Tolerance-default regressions (a change that flips an over-promotion case)
- Cross-variant disagreements (`_simple` and `_robust` returning different counts on clean input)
- Numerical drift from external dependencies (a `MinkowskiReduction.jl` change that affects the integer-matrix test)

It is *not* good at:

- New failure modes that don't appear in the existing AFLOW corpus or 14-Bravais panel
- User-facing API regressions in error messages or warning text (only the warning *firing* is tested, not the wording)
- Performance regressions (timing tests live in `test/runTimingTests.jl` but are excluded from CI)
- Bugs in the test generator itself (`tools/generate_aflow_tests.jl` is not separately validated)

For these classes, the best defense is the experienced reader spotting weirdness in their own use case and filing it as a regression — and the project's pre-1.0 status means breaking changes are expected to be rare but possible.

## See also

- Reference: [`pointGroup`](../reference/point-groups.md), [`spacegroup`](../reference/space-groups.md), [`crystal_system`](../reference/crystals.md)
- Explanation: [Crystal system vs full Bravais](crystal-system-vs-bravais.md), [Over-promotion](over-promotion.md), [Canonicalizing τ](canonicalizing-tau.md)
- External: [AFLOW Library of Crystallographic Prototypes](https://aflow.org/CrystalDatabase/) (Mehl 2017, Hicks 2019, Hicks 2021)
