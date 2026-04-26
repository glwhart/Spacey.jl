# Snap a noisy lattice to symmetry

Real-world lattices (from experimental refinements, DFT relaxations, etc.) usually do not satisfy exact symmetry conditions — a "cubic" cell may have `a = 1.0000` and `b = 1.0001`. [`snapToSymmetry_SVD`](../reference/snap.md) takes a noisy lattice and the symmetry operations found by [`pointGroup`](../reference/point-groups.md) and produces a snapped basis whose symmetry is exact.

## Recipe

### 1. Find the symmetry operations of the noisy lattice

Use a tolerance that comfortably covers the noise.

```julia
using Spacey

# A near-cubic lattice with ~0.1% noise
u = [1.0001, 0.0, 0.0]
v = [0.0, 0.9999, 0.0]
w = [0.0, 0.0, 1.0]

LG, _ = pointGroup(u, v, w; tol=0.01)
length(LG)  # 48 — pointGroup treats this as cubic at this tolerance
```

`LG` is a `Vector{Matrix{Int}}` — exactly the form `snapToSymmetry_SVD` wants.

### 2. Snap

```julia
a, b, c, iops, rops = snapToSymmetry_SVD(u, v, w, LG)
```

Returns a 5-tuple:
- `a, b, c::Vector{Float64}` — the snapped basis vectors. Lengths and inter-vector angles are the symmetry-averaged values; **volume is preserved**.
- `iops::Vector{Matrix{Int}}` — the integer-matrix lattice operations of the snapped lattice (recomputed via `pointGroup_robust` on the snapped basis).
- `rops::Vector{Matrix{Float64}}` — Cartesian rotations of the snapped lattice.

After snapping, `A · iops[i] · inv(A) ≈ rops[i]` to machine precision.

### 3. Use the snapped basis

```julia
norm(a), norm(b), norm(c)   # all equal to within float precision
```

Subsequent symmetry analyses on `[a b c]` no longer need a non-trivial tolerance.

## When to skip snapping

- **Synthetic / clean input** — already exact; snapping is a no-op.
- **High-distortion inputs** where the noise approaches the structural distortion — snapping enforces the higher symmetry, masking the real (lower-symmetry) structure. Use [`verify_stable`](detect-tolerance-dependence.md) on `pointGroup_robust` first to detect this case.

## The internal alternative

A faster but less robust alternative, `Spacey.snapToSymmetry_avg`, averages each basis vector independently over its op images. It is reachable as `Spacey.snapToSymmetry_avg(u, v, w, ops)` but is **not exported** — `snapToSymmetry_SVD` is the recommended path. Use the avg variant only for quick comparisons or when SVD performance becomes a bottleneck (almost never the case for 3×3 lattices).

## See also

- Reference: [`snapToSymmetry_SVD`](../reference/snap.md)
- How-to: [Find a point group](find-pointgroup.md), [Detect tolerance-dependent answers](detect-tolerance-dependence.md)
