# Find a space group

[`spacegroup(c)`](../reference/space-groups.md) returns a `Vector{SpacegroupOp}` of every symmetry operation of crystal `c`, in the user's original basis. Identity is at index 1; the remaining order is unspecified.

## 1. Build a `Crystal`

See [Construct a crystal](construct-a-crystal.md). Example: diamond carbon in its 2-atom primitive cell.

```jldoctest diamond
julia> using Spacey

julia> A = [0.0  0.5  0.5;        # FCC primitive lattice vectors
            0.5  0.0  0.5;
            0.5  0.5  0.0];

julia> r = [0.0  0.25;             # 2 C atoms in the primitive cell
            0.0  0.25;
            0.0  0.25];

julia> c = Crystal(A, r, [:C, :C]; coords=:fractional);
```

## 2. Call `spacegroup`

```jldoctest diamond
julia> ops = spacegroup(c);

julia> length(ops)                 # Fd-3m, full O_h order
48

julia> ops[1] == one(SpacegroupOp)   # identity is always at index 1
true
```

## 3. Read individual operations

Each `SpacegroupOp` has integer rotation `R` and fractional translation `τ` (canonicalized to `[0, 1)`).

```jldoctest diamond
julia> using LinearAlgebra

julia> any(op -> op.R == -Matrix{Int}(I, 3, 3), ops)   # inversion is in the diamond group
true
```

For a non-symmorphic example (translations are not all zero), see [Compose and apply operations](compose-and-apply-ops.md).

## 4. Tune tolerances

Two tolerances control `spacegroup`:

- `lattice_tol` (default `0.01`) — passed to the lattice point-group finder ([`pointGroup`](../reference/point-groups.md)) for the lattice symmetry step.
- `pos_tol` (default `default_pos_tol(c)` ≈ 1% of the characteristic atom separation) — applied when matching atoms under each candidate operation.

```jldoctest
julia> using Spacey, LinearAlgebra

julia> A = Matrix{Float64}(I, 3, 3); # Simple cubic lattice

julia> r = [0.0 0.5; 0.0 0.5; 0.0 0.5]; # Two atoms, the origin and body center

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

For example, rocksalt NaCl in its 2-atom primitive cell returns 48 ops — the full O_h lattice symmetry, with no centering translations because the FCC centering is absorbed into the primitive lattice vectors:

```jldoctest
julia> using Spacey

julia> A = [0.0  0.5  0.5;        # FCC primitive lattice vectors
            0.5  0.0  0.5;
            0.5  0.5  0.0];

julia> r = [0.0  0.5;              # Na at origin, Cl at (½,½,½) primitive coords
            0.0  0.5;
            0.0  0.5];

julia> c = Crystal(A, r, [:Na, :Cl]; coords=:fractional);

julia> length(spacegroup(c))       # Fm-3m, 48 primitive-cell ops
48
```

In the conventional 8-atom FCC cell the same crystal returns 192 = 48 × 4 ops (the four FCC centering translations multiply through).

## See also

- Reference: [`spacegroup`](../reference/space-groups.md), [`SpacegroupOp`](../reference/space-groups.md)
- How-to: [Compose and apply operations](compose-and-apply-ops.md), [Detect tolerance-dependent answers](detect-tolerance-dependence.md)
- Explanation: [Algorithm overview](../explanation/algorithm-overview.md)
