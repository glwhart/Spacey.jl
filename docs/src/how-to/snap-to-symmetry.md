# Snap a noisy lattice to symmetry

Real-world lattices (from experimental refinements, DFT relaxations, etc.) usually do not satisfy exact symmetry conditions — a "cubic" cell may have `a = 1.0000` and `b = 1.0001`. [`snapToSymmetry_SVD`](../reference/snap.md) takes a noisy lattice and the symmetry operations found by [`pointGroup`](../reference/point-groups.md) and produces a snapped basis whose symmetry is exact.

## Recipe

### 1. Find the symmetry operations of the noisy lattice

Use a tolerance that comfortably covers the noise. For example, a near-cubic lattice with ~1% distortion (basis vectors slightly shorter or longer than 1, none of them exactly orthogonal):

```jldoctest snap
julia> using Spacey, LinearAlgebra

julia> u = [1.01, 0.0, 0.0]; v = [0.0, 0.99, 0.0]; w = [0.0, 0.0, 0.999];

julia> LG, _ = pointGroup(u, v, w; tol=0.05);   # tol=0.05 absorbs the 1% noise

julia> length(LG)                                # 48 — pointGroup sees this as cubic
48
```

`LG` is a `Vector{Matrix{Int}}` — exactly the form `snapToSymmetry_SVD` wants. (At the default `tol=0.01`, the same input only finds 8 ops because the 1% distortion now exceeds the tolerance — see [Handle noisy real-world data](handle-noisy-data.md) for how to choose `tol`.)

### 2. Snap

```jldoctest snap
julia> a, b, c, iops, rops = snapToSymmetry_SVD(u, v, w, LG);
```

Returns a 5-tuple:
- `a, b, c::Vector{Float64}` — the snapped basis vectors. Lengths and inter-vector angles are the symmetry-averaged values; **volume is preserved**.
- `iops::Vector{Matrix{Int}}` — the integer-matrix lattice operations of the snapped lattice (recomputed via [`pointGroup`](../reference/point-groups.md) on the snapped basis).
- `rops::Vector{Matrix{Float64}}` — Cartesian rotations of the snapped lattice.

After snapping, `A · iops[i] · inv(A) ≈ rops[i]` to machine precision.

### 3. Use the snapped basis

```jldoctest snap
julia> isapprox(norm(a), norm(b)) && isapprox(norm(b), norm(c))   # snapped lengths equal
true

julia> isapprox(abs(u × v ⋅ w), abs(a × b ⋅ c))                    # volume preserved
true

julia> length(iops)                                                # snapped lattice has full 48 ops
48
```

Subsequent symmetry analyses on `[a b c]` no longer need a non-trivial tolerance — the snapped basis IS exactly cubic, with all three vectors of the same length and all angles exactly 90°.

## When to skip snapping

- **Synthetic / clean input** — already exact; snapping is a no-op.
- **High-distortion inputs** where the noise approaches the structural distortion — snapping enforces the higher symmetry, masking the real (lower-symmetry) structure. Use [`verify_stable`](detect-tolerance-dependence.md) on [`pointGroup`](../reference/point-groups.md) first to detect this case.

## The internal alternative

A faster but less robust alternative, `Spacey.snapToSymmetry_avg`, averages each basis vector independently over its op images. It is reachable as `Spacey.snapToSymmetry_avg(u, v, w, ops)` but is **not exported** — `snapToSymmetry_SVD` is the recommended path. Use the avg variant only for quick comparisons or when SVD performance becomes a bottleneck (only a problem when you are doing tens of thousands of cases over a large corpus of structures).

## See also

- Reference: [`snapToSymmetry_SVD`](../reference/snap.md)
- How-to: [Find a point group](find-pointgroup.md), [Detect tolerance-dependent answers](detect-tolerance-dependence.md)
