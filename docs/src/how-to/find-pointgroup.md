# Find a point group

Spacey has four point-group finders. Pick by the quality of your input.

| Routine | Use when |
|---|---|
| [`pointGroup_robust`](../reference/point-groups.md) | **Default for real-world data.** Tolerance-tunable; emits warnings when results may be sensitive. |
| [`pointGroup`](../reference/point-groups.md) | A matrix-form wrapper around `pointGroup_robust` (passes `tol=0.1` by default). |
| [`pointGroup_fast`](../reference/point-groups.md) | Synthetic / clean input; production speed; strict `isapprox` tolerance. |
| [`pointGroup_simple`](../reference/point-groups.md) | Validation only. Brute-force; obviously correct but slow. |

## Recipe (real input): `pointGroup_robust`

### 1. Have three basis vectors

```jldoctest
julia> using Spacey

julia> u = [1.0, 0, 0]; v = [0.5, √3/2, 0]; w = [0.0, 0, √(8/3)];   # ideal HCP
```

### 2. Call `pointGroup_robust`

It returns the tuple `(LG, G)`:

- `LG` — operations as integer matrices in lattice coordinates.
- `G` — the same operations as Cartesian rotations.

The two are related by `A · LG[i] · inv(A) ≈ G[i]` where `A = [u v w]`.

```jldoctest hcp
julia> using Spacey

julia> u = [1.0, 0, 0]; v = [0.5, √3/2, 0]; w = [0.0, 0, √(8/3)];

julia> LG, G = pointGroup_robust(u, v, w);

julia> length(LG)   # hexagonal point group has 24 ops
24
```

### 3. Tune `tol` if needed

The default `tol = 0.01` is appropriate for clean numerical input (synthetic, post-snap, or noise ≪ 1%). For experimental data, scale to the noise level. **Looser tolerance accepts more candidate operations and risks over-promotion** (reporting a higher symmetry than is actually present).

```jldoctest
julia> using Spacey

julia> u = [1.0, 0, 0]; v = [0, 1.0, 0]; w = [0, 0, 1.0];

julia> length(pointGroup_robust(u, v, w; tol=1e-6)[1])   # tight: cubic, 48
48
```

## Matrix wrapper

For a 3×3 matrix `A` whose columns are the basis vectors, [`pointGroup`](../reference/point-groups.md) is a convenience:

```jldoctest
julia> using Spacey, LinearAlgebra

julia> A = Matrix{Float64}(I, 3, 3);

julia> length(pointGroup(A)[1])
48
```

Equivalent to `pointGroup_robust(A[:,1], A[:,2], A[:,3]; tol=0.1)`. Note the looser default tolerance — chosen so that `pointGroup(A)` "just works" for typical user matrices.

## Fast finder for clean input

When you control the input (synthetic lattices, validation suites), [`pointGroup_fast`](../reference/point-groups.md) skips the tolerance machinery and runs faster. It returns just the integer-matrix list (no `(LG, G)` tuple):

```jldoctest
julia> using Spacey

julia> u = [1.0, 0, 0]; v = [0.5, √3/2, 0]; w = [0.0, 0, √(8/3)];

julia> length(pointGroup_fast(u, v, w))
24
```

## Validation finder

[`pointGroup_simple`](../reference/point-groups.md) brute-forces over every integer matrix in `[-1, 1]³`. It is slow but obviously correct, and is used only to cross-validate the other finders on clean input.

## See also

- Explanation: [Why Minkowski reduction](../explanation/why-minkowski.md), [Over-promotion](../explanation/over-promotion.md)
- How-to: [Detect tolerance-dependent answers](detect-tolerance-dependence.md), [Handle noisy real-world data](handle-noisy-data.md)
