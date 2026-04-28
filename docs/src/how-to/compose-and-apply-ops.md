# Compose and apply operations

A [`SpacegroupOp`](../reference/space-groups.md) is the operation `r ↦ R·r + τ` in lattice (fractional) coordinates, where `R` is a rotation/reflection of the unit cell and `τ` is a translation (fractional, not a full lattice translation, and often zero). The struct supports composing two symmetry elements, or finding a symmetry's inverse. It also supports applying the symmetry to atomic positions, or checking for equality of two symmetry operations (modulo the lattice). It marks the identity element. A helper function is provided that converts operations to Cartesian form (internally they are stored as integer matrices). This page is the cookbook for using them.

!!! note "Roadmap"
    Applying an op to an entire `Crystal` (and checking equality of two crystals or lattices under an operation) is currently a manual workflow — see §2 below. A `SpacegroupOp` method that takes a `Crystal` directly, plus a `crystal_equals_under_op` predicate, are on the roadmap. Until then, the manual recipe in §2 is the path.

## 1. Get an op (or build one by hand)

Most ops come from [`spacegroup(c)`](../reference/space-groups.md). For example, diamond carbon in its 2-atom primitive cell:

```jldoctest setup
julia> using Spacey

julia> A = [0.0  0.5  0.5;        # FCC primitive lattice vectors
            0.5  0.0  0.5;
            0.5  0.5  0.0];

julia> r = [0.0  0.25;             # 2 C atoms (origin and (¼,¼,¼))
            0.0  0.25;
            0.0  0.25];

julia> c = Crystal(A, r, [:C, :C]; coords=:fractional);

julia> ops = spacegroup(c);

julia> length(ops)                 # Fd-3m, full O_h order
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

## 2. Apply an op to a fractional position (and to a whole crystal)

The struct is callable. `op(r)` returns the image of fractional position `r`, folded back to `[0, 1)`:

```jldoctest
julia> using Spacey

julia> R = [0 -1 0; 1 0 0; 0 0 1];   # 4-fold about z, lattice coords

julia> op = SpacegroupOp(R, [0.5, 0.0, 0.0]);   # add a ½-translation along x

julia> op([0.25, 0.10, 0.30])        # not just origin: rotation is visible
3-element Vector{Float64}:
 0.4
 0.25
 0.3
```

For a whole [`Crystal`](../reference/crystals.md), apply the op to each atom column of `c.r`. Spacey doesn't ship a single-shot `op(c::Crystal)` method yet (see the roadmap note above); the manual form is short:

```jldoctest cscl
julia> using Spacey, LinearAlgebra

julia> A = Matrix{Float64}(I, 3, 3);

julia> r = [0.0 0.5; 0.0 0.5; 0.0 0.5];   # CsCl: Cs at origin, Cl at body-center

julia> c = Crystal(A, r, [:Cs, :Cl]; coords=:fractional);

julia> sg = spacegroup(c);

julia> hcat([sg[2](c.r[:, i]) for i in 1:size(c.r, 2)]...)   # apply sg[2] to every atom
3×2 Matrix{Float64}:
 0.0  0.5
 0.0  0.5
 0.0  0.5
```

Because `sg[2]` is a true symmetry of CsCl, the *set* of image positions is identical to the original set (each Cs and Cl maps to its own type's position) — that's exactly what makes it a symmetry.

## 3. Compose with `*`

`op1 * op2` is function composition: **apply `op2` first, then `op1`**. Algebraically, `(R₁, τ₁) ∘ (R₂, τ₂) = (R₁·R₂, R₁·τ₂ + τ₁)`. The result is a new `SpacegroupOp`. This matches Julia's `∘` and standard function-composition semantics where `(f ∘ g)(x) = f(g(x))` — so `op1 * op2` and `op1 ∘ op2` agree on element ordering.

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

Equality is exact field-by-field on `(R, τ)`. Three things to verify:

```jldoctest
julia> using Spacey, LinearAlgebra

julia> I3 = Matrix{Int}(I, 3, 3);

julia> R = [0 -1 0; 1 0 0; 0 0 1];          # 4-fold about z

julia> SpacegroupOp(R, [0.0, 0, 0]) == SpacegroupOp(R, [1.0, 0, 0])  # τ folded mod 1
true

julia> SpacegroupOp(R, [0.5, 0, 0]) == SpacegroupOp(R, [0.0, 0, 0])  # different τ → not equal
false

julia> SpacegroupOp(R, [0.5, 0, 0]) == SpacegroupOp(I3, [0.5, 0, 0])  # different R → not equal
false
```

The matching `hash` method makes ops safe as `Set` and `Dict` keys. This is the equality `spacegroup` itself uses internally — to enforce identity-at-index-1 in the returned vector and to deduplicate ops while verifying group closure.

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

The result is a `Tuple{Matrix{Float64},Vector{Float64}}`, **not** a `SpacegroupOp` — `SpacegroupOp.R` must be `Matrix{Int}`, but the Cartesian rotation is in general not integer-valued. For an orthonormal cubic lattice (`A = I`), the Cartesian rotation `R_cart = R` literally. For a non-orthonormal lattice (e.g. hexagonal, where the basis vectors are not mutually perpendicular), `R_cart = A·R·A⁻¹` is the same operation expressed in the orthogonal Cartesian frame — its entries depend on the specific basis even though the underlying symmetry is intrinsic to the lattice.

## See also

- Reference: [`SpacegroupOp`](../reference/space-groups.md), [`toCartesian`](../reference/space-groups.md)
- How-to: [Find a space group](find-spacegroup.md)
- Explanation: [Canonicalizing τ](../explanation/canonicalizing-tau.md)
