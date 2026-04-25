# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What This Package Does

Spacey.jl finds the point group and space group of crystal lattices. The central challenge is **robustness to finite-precision floating-point errors** — lattice vectors from real materials codes are noisy, high aspect-ratio, and may not satisfy exact symmetry conditions. The algorithms are specifically designed for these conditions.

## Commands

**Run the full test suite:**
```bash
julia --project test/runtests.jl
# or from the Julia REPL:
using Pkg; Pkg.test()
```

**Run a single test block** (use a `@testset` label from `test/runtests.jl`):
```julia
using Pkg; Pkg.test(test_args=["BravaisLatticeList"])
```

**Run timing benchmarks** (not part of CI):
```bash
julia --project test/runTimingTests.jl
```

**Build documentation:**
```bash
cd docs && julia --project make.jl
```

## Architecture

All package code lives in the single file `src/Spacey.jl`. There is no code split across multiple source files.

### Exported API

| Function | Purpose |
|---|---|
| `pointGroup_robust(u, v, w; tol=0.01, verify_stable=false)` | Main algorithm for noisy/real-world inputs; opt-in stability check warns on tolerance-dependent answers |
| `pointGroup(A; tol=0.1)` | Matrix wrapper for `pointGroup_robust` |
| `pointGroup_fast(a1, a2, a3)` | Production-speed variant, strict tolerances |
| `pointGroup_simple(a1, a2, a3)` | Naive brute-force, used only for validation |
| `snapToSymmetry_SVD(u, v, w, ops)` | Snap noisy lattice to exact symmetry via SVD |
| `snapToSymmetry_avg(v1, v2, v3, ops)` | Snap via averaging over group operations |
| `isagroup(members)` | Verify a set of matrices forms a group |
| `Crystal(A, r, types; coords)` | Construct a crystal (required `coords` kwarg: `:fractional` or `:cartesian`) |
| `spacegroup(c; lattice_tol=0.01, pos_tol=default_pos_tol(c), verify_stable=false)` | Find all `(R, τ)` space-group operations; opt-in stability check |
| `isSpacegroupOp(R, τ, c; tol)` | Check if `(R, τ)` is a symmetry of crystal `c` |
| `fractional(c)`, `cartesian(c)` | Atomic positions in the respective basis |
| `default_pos_tol(c)` | Default position tolerance: `0.01·(V/N)^(1/3)` |
| `crystal_system(A)` / `crystal_system(c)` | Identify Bravais system from lattice holohedry (one of `:triclinic`, `:monoclinic`, `:orthorhombic`, `:tetragonal`, `:trigonal`, `:hexagonal`, `:cubic`) |
| `toCartesian(op, A)` | Convert a `SpacegroupOp` to Cartesian `(R, τ)` tuple |

`pointGroup_robust` returns a tuple `(LG, G)` where `LG` contains operations in lattice coordinates and `G` contains Cartesian rotations. These are related by `A * LG[i] * inv(A) == G[i]`.

`spacegroup(c)` returns a `Vector{SpacegroupOp}` with `R::Matrix{Int}` and `τ::Vector{Float64}`, both in the user's original basis. Identity is guaranteed at index 1; remaining order is unspecified. The struct supports composition (`*`), inversion (`inv`), application to a fractional position (`op(r)`), and mod-1 equality (via canonicalised `τ`).

### Core Algorithm (pointGroup_robust)

1. **Minkowski-reduce** the input basis (delegates to `MinkowskiReduction.jl`) — this is provably sufficient to define the candidate search space.
2. **Generate candidates**: all integer-coefficient vectors from the 27-point {-1,0,1}³ grid applied to the reduced basis.
3. **Filter by**: norm match → volume conservation → integer-matrix test (U = inv(A)*B must be near-integer).
4. **Group closure**: from candidate operations, find the largest subset that closes under multiplication. Tries group sizes [48, 24, 16, 12, 8, 4, 2] in decreasing order.
5. **Stability check (opt-in)**: with `verify_stable=true`, the algorithm re-runs at `tol/1000` and emits a `@warn` if the group size differs — flagging that the lattice is near a symmetry boundary and the answer is tolerance-dependent. Default is `false` (silent, current behaviour). The same opt-in is available on `spacegroup` (re-run at `pos_tol/1000`).

Key implementation detail: inputs are normalized by `∛|det(A)|` before comparisons so all tolerances operate at unit scale.

### Validation Strategy

Three algorithm variants exist specifically for cross-validation:
- `pointGroup_simple`: brute-force over all integer matrices in [-1,1]³ — slow but obviously correct
- `pointGroup_fast`: optimized filter pipeline, strict tolerances — production speed
- `pointGroup_robust`: tolerance-tunable, designed for real-world noisy input

Tests verify all three agree on exact (non-noisy) inputs for all 14 Bravais lattice types.

### Test Suite Structure

`test/runtests.jl` covers, in order:
- All 14 Bravais lattice types with expected group sizes (e.g., cubic=48, triclinic=2)
- LG/G consistency: `A * LG[i] * inv(A) ≈ G[i]` for all returned operations
- Rotation invariance: applying arbitrary 3D rotations must not change group size
- Noise tolerance: small perturbations (~1e-8) must not crash or return wrong group
- High aspect ratio: AR = 256, 500, 512, 1024 are explicitly tested
- Snap-to-symmetry: verifies volume conservation and correct group size after snapping
- Parametric sweep: (tolerance, noise level) grid over all Bravais types
- **Near-boundary tetragonal**: pins point-group over-promotion at loose `tol` and the `verify_stable` warning (see `test/nearMissBoundary.jl` for the diagnostic heatmap)
- **Crystal/space-group infrastructure** (Phase 1): `Crystal` constructor validation, `isSpacegroupOp` trivial cases, `default_pos_tol` formula, helper accessors
- **`SpacegroupOp` methods** (Phase 2): composition, inverse, callable, mod-1 equality, `toCartesian`
- **`spacegroup` Phase 2 core cases**: simple cubic, tetragonal, orthorhombic, triclinic, CsCl
- **`spacegroup` Phase 3 known crystals**: NaCl (192), diamond (192), HCP (24)
- **`spacegroup` AFLOW Part 1 (seed)**: hand-curated 14-prototype subset with structural sub-checks (closure, isSpacegroupOp cross-check, non-symmorphic signatures)
- **`spacegroup` AFLOW Parts 1, 2, 3 (full corpus)**: every prototype from the three published AFLOW papers (286 + 299 + 510 = 1095) auto-generated by `tools/generate_aflow_tests.jl` and tested against the order encoded in the prototype label. Deviations are flipped to `@test_broken` to pin current behaviour
- **`crystal_system`**: 14-Bravais smoke test plus AFLOW cross-check on all three parts (independent invariant: lattice holohedry order maps to crystal system)
- **Phase 4 near-boundary crystal**: BaTiO₃-style ferroelectric, pins `verify_stable` warning behaviour for `spacegroup` (heatmap in `test/nearMissBoundaryCrystal.jl`)

Total: ~3,400 tests passing on a clean run.

### Known Stubs

- `debugging.jl` in `src/` contains exploratory/interactive code, not part of the module

## CI

GitHub Actions runs tests on Ubuntu, macOS, and Windows with Julia 1.11 (x64). Coverage is uploaded to codecov. Timing tests are deliberately excluded from CI due to VM timing variability. CompatHelper runs daily to check dependency updates.

## Key Dependencies

- `StatsBase`: used for averaging in `snapToSymmetry_avg`
- `LinearAlgebra`: core matrix operations throughout

## MinkowskiReduction.jl

This dependency is central to Spacey.jl's correctness. After Minkowski reduction, the symmetry operations of a lattice can be found by searching only the {-1,0,1}³ integer grid — this is why the candidate search space is provably finite and complete.

**API used by Spacey.jl:**

```julia
# Vector form — returns 4-tuple: (u, v, w, iterations)
u, v, w, nsteps = minkReduce(U, V, W)

# Matrix form — returns just the reduced matrix (no iteration count)
Areduced = minkReduce(A)
```

The vector form always returns 4 values. The matrix form drops the iteration count. Do not conflate these two signatures.

**Hard limit:** `minkReduce` throws an error (not a warning) if reduction has not converged after 15 iterations. This has never been hit in practice for physically realistic lattices, but adversarial inputs can trigger it.

**Testing utilities available from MinkowskiReduction** (all exported, usable in tests):

| Function | Purpose |
|---|---|
| `isMinkReduced(U,V,W)` or `isMinkReduced(M)` | Verify a basis is fully reduced (useful for assertions) |
| `orthogonalityDefect(a,b,c)` | `∏‖vᵢ‖ / |det(A)|` — ≥ 1, equality iff the basis is orthogonal. Distinct from Spacey's internal `aspectRatio`: see note below. |
| `DeviousMat(n)` | Generate an adversarial unimodular 3×3 matrix requiring many reduction steps; useful for stress-testing |
| `RandUnimodMat3(k=10)` | Random unimodular 3×3 matrix; use to generate random-but-valid lattice transformations in tests |
| `RandUnimodMat2(n)` | Same for 2×2 |
| `isPermutationMatrix(M)` | Check if M is a signed permutation of the identity |
| `GaussReduce(U,V)` | 2D Gauss reduction (used internally by `minkReduce`) |

**`orthogonalityDefect` vs `aspectRatio` — these are NOT analogous.** They measure different things and a basis can score badly on one while scoring perfectly on the other:

- `orthogonalityDefect = ∏‖vᵢ‖ / |det(A)|` measures **angular** non-orthogonality. Equal to 1 (the minimum) iff the basis vectors are mutually perpendicular, regardless of their lengths.
- `aspectRatio = max‖vᵢ‖ / min‖vᵢ‖` (after Minkowski reduction) measures **length disparity** between basis vectors, regardless of the angles between them.

Counter-example: the basis `(1,0,0), (0,1,0), (0,0,1000)` has `orthogonalityDefect = 1` (perfectly orthogonal) but `aspectRatio = 1000` (extremely elongated). Conversely, three equal-length but skewed vectors can have `aspectRatio = 1` and a large `orthogonalityDefect`.

For Spacey's purposes the two are not symmetric concerns. Spacey calls `minkReduce` before any symmetry analysis, and Minkowski reduction is exactly the operation that drives `orthogonalityDefect` down to whatever the lattice allows — so a poor pre-reduction defect is not a problem the user has to manage. Aspect ratio, by contrast, is a property of the *lattice itself* (e.g. a tetragonal cell with c ≫ a stays elongated no matter how it is reduced) and is what compresses the integer-grid filter's discriminative power. That's why `pointGroup_robust` warns on `aspectRatio > 100` but does not warn on `orthogonalityDefect`.
