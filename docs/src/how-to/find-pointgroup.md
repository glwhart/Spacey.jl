# Find a point group

[`pointGroup`](../reference/point-groups.md) is the single public entry point. It accepts either three basis vectors or a 3×3 matrix and returns the tuple `(LG, G)` (lattice-coordinate integer matrices and their Cartesian rotations).

## 1. Have three basis vectors (or a matrix)

```jldoctest
julia> using Spacey

julia> u = [1.0, 0, 0]; v = [0.5, √3/2, 0]; w = [0.0, 0, √(8/3)];   # ideal HCP
```

## 2. Call `pointGroup`

```jldoctest hcp
julia> using Spacey

julia> u = [1.0, 0, 0]; v = [0.5, √3/2, 0]; w = [0.0, 0, √(8/3)];

julia> LG, G = pointGroup(u, v, w);

julia> length(LG)   # hexagonal point group has 24 ops
24
```

The two halves of the tuple are related by `A · LG[i] · inv(A) ≈ G[i]` where `A = [u v w]` — same operation, different basis.

For a lattice given as a 3×3 matrix, pass it directly:

```jldoctest
julia> using Spacey, LinearAlgebra

julia> length(pointGroup(Matrix{Float64}(I, 3, 3))[1])
48
```

## 3. Tune `tol` for noisy input

The default `tol = 0.01` is appropriate for clean numerical input (synthetic, post-snap, or noise ≪ 1%). For experimental data, scale to the noise level. **Looser tolerance accepts more candidate operations and risks over-promotion** (reporting a higher symmetry than is actually present).

```jldoctest
julia> using Spacey

julia> u = [1.0, 0, 0]; v = [0, 1.0, 0]; w = [0, 0, 1.0];

julia> length(pointGroup(u, v, w; tol=1e-6)[1])   # tight: cubic, 48
48
```

## 4. Detect tolerance-dependent answers

Pass `verify_stable=true` to re-run the algorithm at `tol/1000` and emit a `@warn` if the operation count differs. See [Detect tolerance-dependent answers](detect-tolerance-dependence.md) for the full pattern.

## Rotation invariance

The point group is intrinsic to the lattice — applying an arbitrary 3D rotation to the basis must not change the result. Take a rhombohedral basis `(1,1,2), (1,2,1), (2,1,1)`: it has 12 ops (D₃d), with the 3-fold axis along the body diagonal `(1,1,1)`. Rotating to an arbitrary orientation should still give 12. Pass through `minkReduce` first because (a) the raw rhombohedral basis above is not Minkowski-reduced — its differences are shorter than the basis vectors — and (b) the rotated basis generally won't be reduced either:

```jldoctest
julia> using Spacey, MinkowskiReduction

julia> u = [1.0, 1, 2]; v = [1.0, 2, 1]; w = [2.0, 1, 1];   # rhombohedral

julia> u, v, w = minkReduce(u, v, w)[1:3];

julia> length(pointGroup(u, v, w)[1])
12

julia> u, v, w = [1.0, 1, 2], [1.0, 2, 1], [2.0, 1, 1];   # reset to the original

julia> u, v, w = Spacey.threeDrotation(u, v, w, π/3, π/5, π/7);   # rotate by Euler angles (π/3, π/5, π/7)

julia> u, v, w = minkReduce(u, v, w)[1:3];

julia> length(pointGroup(u, v, w)[1])
12
```

Same group order, regardless of orientation. (`Spacey.threeDrotation` is an internal test helper — see the variants table below for the convention on internal names.)

## Variants reachable via the qualified name

`pointGroup` delegates to the internal robust finder (`pointGroup_robust`). Three variants are kept inside the package and reachable as `Spacey.<name>` for specialized use cases — none of them is exported, none is part of the public API contract:

| Internal name | Use case |
|---|---|
| `Spacey.pointGroup_robust` | The exact function `pointGroup` delegates to. Reach directly only when explicitly comparing variants. |
| `Spacey.pointGroup_fast` | Synthetic / clean input, production speed; strict `isapprox` tolerance with no `tol` knob. |
| `Spacey.pointGroup_simple` | Validation only. Brute-force; obviously correct but slow. |

In normal use, prefer `pointGroup`.

## See also

- Reference: [`pointGroup`](../reference/point-groups.md)
- Explanation: [Why Minkowski reduction](../explanation/why-minkowski.md), [Over-promotion](../explanation/over-promotion.md)
- How-to: [Detect tolerance-dependent answers](detect-tolerance-dependence.md), [Handle noisy real-world data](handle-noisy-data.md)
