# Compose and apply operations

A [`SpacegroupOp`](../reference/space-groups.md) is the operation `r ↦ R·r + τ` in lattice (fractional) coordinates. The struct supports composition, inversion, application to a position, equality (modulo the lattice), an identity element, and conversion to Cartesian form. This page is the cookbook for using them.

## 1. Get an op (or build one by hand)

Most ops come from [`spacegroup(c)`](../reference/space-groups.md):

```jldoctest setup
julia> using Spacey, LinearAlgebra

julia> A = Matrix{Float64}(I, 3, 3);

julia> c = Crystal(A, reshape([0.0, 0.0, 0.0], 3, 1), [:X]; coords=:fractional);

julia> ops = spacegroup(c);

julia> length(ops)
48
```

You can also construct one directly. The constructor canonicalizes `τ` to `[0, 1)` automatically — `1.0`, `2.5`, and `-0.5` all fold to a value in that interval:

```jldoctest
julia> using Spacey

julia> R = [0 -1 0; 1 0 0; 0 0 1];   # 4-fold rotation about z (lattice coords)

julia> op = SpacegroupOp(R, [2.5, -0.5, 1.25]);   # any real τ accepted

julia> op.τ                                       # stored canonicalized to [0,1)
3-element Vector{Float64}:
 0.5
 0.5
 0.25
```

The identity is available via Julia's standard `one`:

```jldoctest
julia> using Spacey

julia> e = one(SpacegroupOp);

julia> e.R
3×3 Matrix{Int64}:
 1  0  0
 0  1  0
 0  0  1

julia> e.τ
3-element Vector{Float64}:
 0.0
 0.0
 0.0
```

## 2. Apply an op to a fractional position

The struct is callable. `op(r)` returns the image position, folded back to `[0, 1)`:

```jldoctest
julia> using Spacey

julia> R = [0 -1 0; 1 0 0; 0 0 1];

julia> op = SpacegroupOp(R, [0.5, 0.0, 0.0]);

julia> op([0.0, 0.0, 0.0])         # origin maps to (0.5, 0, 0)
3-element Vector{Float64}:
 0.5
 0.0
 0.0
```

## 3. Compose with `*`

`op1 * op2` is function composition: **apply `op2` first, then `op1`**. Algebraically, `(R₁, τ₁) ∘ (R₂, τ₂) = (R₁·R₂, R₁·τ₂ + τ₁)`. The result is a new `SpacegroupOp`.

```jldoctest
julia> using Spacey

julia> R = [0 -1 0; 1 0 0; 0 0 1];

julia> op1 = SpacegroupOp(R, [0.5, 0.0, 0.0]);

julia> op2 = SpacegroupOp(R, [0.0, 0.5, 0.0]);

julia> op12 = op1 * op2;

julia> op12.R == R * R       # rotational parts compose
true
```

Composition is associative, and the identity behaves as expected:

```jldoctest
julia> using Spacey

julia> R = [0 -1 0; 1 0 0; 0 0 1];

julia> e = one(SpacegroupOp);

julia> op = SpacegroupOp(R, [0.5, 0.0, 0.0]);

julia> e * op == op
true

julia> op * e == op
true
```

## 4. Invert with `inv`

`inv(op)` returns the operation that undoes `op`. Algebraically, `(R, τ)⁻¹ = (R⁻¹, -R⁻¹·τ)`. `R⁻¹` is integer because `|det R| = 1` for any lattice rotation; the constructor enforces this and errors otherwise.

```jldoctest
julia> using Spacey

julia> R = [0 -1 0; 1 0 0; 0 0 1];

julia> op = SpacegroupOp(R, [0.25, 0.0, 0.0]);

julia> op * inv(op) == one(SpacegroupOp)
true

julia> inv(inv(op)) == op
true
```

## 5. Compare for equality

Equality is exact field-by-field. Because `τ` is canonicalized at construction, ops that differ only by a lattice translation (e.g. `τ = [1.0, 0, 0]` vs `τ = [0.0, 0, 0]`) compare equal:

```jldoctest
julia> using Spacey, LinearAlgebra

julia> I3 = Matrix{Int}(I, 3, 3);

julia> SpacegroupOp(I3, [0.0, 0, 0]) == SpacegroupOp(I3, [1.0, 0, 0])
true

julia> SpacegroupOp(I3, [0.0, 0, 0]) == SpacegroupOp(I3, [0.5, 0, 0])
false
```

The matching `hash` method makes ops safe as `Set` and `Dict` keys.

## 6. Convert to Cartesian form

`spacegroup` returns ops in the user's lattice coordinates. To get the Cartesian rotation `R_c` and translation `τ_c` (e.g. for plotting or for combining with other Cartesian transforms), use [`toCartesian`](../reference/space-groups.md):

```jldoctest
julia> using Spacey, LinearAlgebra

julia> A = Matrix{Float64}(I, 3, 3);

julia> R_cart, τ_cart = toCartesian(one(SpacegroupOp), A);

julia> R_cart
3×3 Matrix{Float64}:
 1.0  0.0  0.0
 0.0  1.0  0.0
 0.0  0.0  1.0

julia> τ_cart
3-element Vector{Float64}:
 0.0
 0.0
 0.0
```

The result is a `Tuple{Matrix{Float64},Vector{Float64}}`, **not** a `SpacegroupOp` — `SpacegroupOp.R` must be `Matrix{Int}`, but the Cartesian rotation is in general not integer-valued. For non-orthonormal lattices, the conversion `(A·R·A⁻¹, A·τ)` produces a Cartesian rotation that depends on the lattice.

## See also

- Reference: [`SpacegroupOp`](../reference/space-groups.md), [`toCartesian`](../reference/space-groups.md)
- How-to: [Find a space group](find-spacegroup.md)
- Explanation: [Canonicalizing τ](../explanation/canonicalizing-tau.md)
