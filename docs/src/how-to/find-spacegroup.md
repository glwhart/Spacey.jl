# Find a space group

[`spacegroup(c)`](../reference/space-groups.md) returns a `Vector{SpacegroupOp}` of every symmetry operation of crystal `c`, in the user's original basis. Identity is at index 1; the remaining order is unspecified.

## 1. Build a `Crystal`

See [Construct a crystal](construct-a-crystal.md). Example: simple cubic with one atom at the origin.

```jldoctest sc
julia> using Spacey, LinearAlgebra

julia> A = Matrix{Float64}(I, 3, 3);

julia> r = reshape([0.0, 0.0, 0.0], 3, 1);

julia> c = Crystal(A, r, [:X]; coords=:fractional);
```

## 2. Call `spacegroup`

```jldoctest sc
julia> ops = spacegroup(c);

julia> length(ops)
48

julia> ops[1] == one(SpacegroupOp)   # identity is always at index 1
true
```

## 3. Read individual operations

Each `SpacegroupOp` has integer rotation `R` and fractional translation `τ` (canonicalised to `[0, 1)`).

```jldoctest sc
julia> using LinearAlgebra

julia> any(op -> op.R == -Matrix{Int}(I, 3, 3), ops)   # inversion is in the cubic group
true
```

For a non-symmorphic example (translations are not all zero), see [Compose and apply operations](compose-and-apply-ops.md).

## 4. Tune tolerances

Two tolerances control `spacegroup`:

- `lattice_tol` (default `0.01`) — passed to `pointGroup_robust` for the lattice symmetry step.
- `pos_tol` (default `default_pos_tol(c)` ≈ 1% of the characteristic atom separation) — applied when matching atoms under each candidate operation.

```jldoctest
julia> using Spacey, LinearAlgebra

julia> A = Matrix{Float64}(I, 3, 3);

julia> r = [0.0 0.5; 0.0 0.5; 0.0 0.5];

julia> c_cscl = Crystal(A, r, [:Cs, :Cl]; coords=:fractional);

julia> length(spacegroup(c_cscl; pos_tol=1e-6))   # CsCl: 48 ops
48
```

The default `pos_tol` is computed from the cell volume and atom count via [`default_pos_tol`](../reference/crystals.md). Override only when you have a specific reason — see [Handle noisy real-world data](handle-noisy-data.md).

## 5. Sanity check the result

The number of operations must match the order of one of the 230 space groups. Quick checks:

- The op count is always a multiple of the lattice point group's order.
- For a primitive cell, the op count equals the order of the space group's point group times any internal translations.
- Conventional (non-primitive) cells multiply by the centering factor (×2 for I/C, ×4 for F).

For example, NaCl in its 8-atom conventional FCC cell returns 192 ops (48 point ops × 4 centering translations).

## See also

- Reference: [`spacegroup`](../reference/space-groups.md), [`SpacegroupOp`](../reference/space-groups.md)
- How-to: [Compose and apply operations](compose-and-apply-ops.md), [Detect tolerance-dependent answers](detect-tolerance-dependence.md)
- Explanation: [Algorithm overview](../explanation/algorithm-overview.md)
