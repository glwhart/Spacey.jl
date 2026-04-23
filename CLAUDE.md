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
| `pointGroup_robust(u, v, w; tol=0.01)` | Main algorithm for noisy/real-world inputs |
| `pointGroup(A; tol=0.1)` | Matrix wrapper for `pointGroup_robust` |
| `pointGroup_fast(a1, a2, a3)` | Production-speed variant, strict tolerances |
| `pointGroup_simple(a1, a2, a3)` | Naive brute-force, used only for validation |
| `snapToSymmetry_SVD(u, v, w, ops)` | Snap noisy lattice to exact symmetry via SVD |
| `snapToSymmetry_avg(v1, v2, v3, ops)` | Snap via averaging over group operations |
| `aspectRatio(a1, a2, a3)` | Compute lattice aspect ratio |
| `isagroup(members)` | Verify a set of matrices forms a group |
| `threeDrotation(...)` | Generate rotated test lattices |

`pointGroup_robust` returns a tuple `(LG, G)` where `LG` contains operations in lattice coordinates and `G` contains Cartesian rotations. These are related by `A * LG[i] * inv(A) == G[i]`.

### Core Algorithm (pointGroup_robust)

1. **Minkowski-reduce** the input basis (delegates to `MinkowskiReduction.jl`) — this is provably sufficient to define the candidate search space.
2. **Generate candidates**: all integer-coefficient vectors from the 27-point {-1,0,1}³ grid applied to the reduced basis.
3. **Filter by**: norm match → volume conservation → integer-matrix test (U = inv(A)*B must be near-integer).
4. **Group closure**: from candidate operations, find the largest subset that closes under multiplication. Tries group sizes [48, 24, 16, 12, 8, 4, 2] in decreasing order.
5. **Degenerate tolerance handling**: if best group is not uniquely determined, try multiple tolerance levels and take the largest valid group.

Key implementation detail: inputs are normalized by `∛|det(A)|` before comparisons so all tolerances operate at unit scale.

### Validation Strategy

Three algorithm variants exist specifically for cross-validation:
- `pointGroup_simple`: brute-force over all integer matrices in [-1,1]³ — slow but obviously correct
- `pointGroup_fast`: optimized filter pipeline, strict tolerances — production speed
- `pointGroup_robust`: tolerance-tunable, designed for real-world noisy input

Tests verify all three agree on exact (non-noisy) inputs for all 14 Bravais lattice types.

### Test Suite Structure

`test/runtests.jl` covers:
- All 14 Bravais lattice types with expected group sizes (e.g., cubic=48, triclinic=2)
- LG/G consistency: `A * LG[i] * inv(A) ≈ G[i]` for all returned operations
- Rotation invariance: applying arbitrary 3D rotations must not change group size
- Noise tolerance: small perturbations (~1e-8) must not crash or return wrong group
- High aspect ratio: AR = 256, 500, 512, 1024 are explicitly tested
- Snap-to-symmetry: verifies volume conservation and correct group size after snapping
- Parametric sweep: (tolerance, noise level) grid over all Bravais types

### Known Stubs

- `spacegroup(A)` returns `true` — not implemented, reserved for future work
- `Crystal` struct is defined but unused
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
| `orthogonalityDefect(a,b,c)` | `∏‖vᵢ‖ / |det(A)|` — measures deviation from orthogonality; analogous to Spacey's `aspectRatio` |
| `DeviousMat(n)` | Generate an adversarial unimodular 3×3 matrix requiring many reduction steps; useful for stress-testing |
| `RandUnimodMat3(k=10)` | Random unimodular 3×3 matrix; use to generate random-but-valid lattice transformations in tests |
| `RandUnimodMat2(n)` | Same for 2×2 |
| `isPermutationMatrix(M)` | Check if M is a signed permutation of the identity |
| `GaussReduce(U,V)` | 2D Gauss reduction (used internally by `minkReduce`) |

**Relationship between `orthogonalityDefect` and `aspectRatio`:** Both measure how "skew" a lattice is. `orthogonalityDefect` is ≥ 1 with equality only for orthogonal bases. Spacey.jl's `aspectRatio` is the ratio of longest to shortest basis vector. High values of either signal cases where the algorithm may need a larger tolerance `tol`.
