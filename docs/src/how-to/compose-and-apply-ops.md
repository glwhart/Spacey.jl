# Compose and apply operations

A [`SpacegroupOp`](../reference/space-groups.md) is the operation `r ↦ R·r + τ` in lattice (fractional) coordinates, where `R` is a rotation/reflection of the unit cell and `τ` is a translation (fractional, not a full lattice translation, and often zero). Functions in the package support composing two symmetry elements, or finding a symmetry's inverse. The package also supports applying a symmetry (a `SpacegroupOp`) to atomic positions or checking for equality of two symmetry operations (modulo the lattice). It marks the identity element. A helper function is provided that converts operations to Cartesian form (internally they are stored as integer matrices). This page is the cookbook for using `SpacegroupOp`s.

!!! note "Roadmap"
    Applying an op to an entire `Crystal` (and checking equality of two crystals or lattices under an operation) is currently a manual workflow — see §2 below. A `SpacegroupOp` method that takes a `Crystal` directly, plus a `crystal_equals_under_op` predicate, are on the roadmap. Until then, the manual recipe in §2 is the path.

## 1. Get an op (or build one by hand)

Most ops come from finding the symmetry of a crystal `c`: [`spacegroup(c)`](../reference/space-groups.md). For example, diamond carbon in its 2-atom primitive cell:

```jldoctest setup
julia> using Spacey

julia> A = [0.0  0.5  0.5;        # FCC primitive lattice vectors
            0.5  0.0  0.5;
            0.5  0.5  0.0];

julia> r = [0.0  0.25;             # 2 carbon atoms (origin and (¼,¼,¼))
            0.0  0.25;
            0.0  0.25];

julia> c = Crystal(A, r, [:C, :C]; coords=:fractional); # :C is the label for the atoms

julia> ops = spacegroup(c);

julia> length(ops)                 # Fd-3m, full O_h order
48

julia> ops[1]                      # the identity is always at index 1
SpacegroupOp(R = [1 0 0; 0 1 0; 0 0 1], τ = [0.0, 0.0, 0.0])
```

The remaining 47 ops appear in an unspecified order. Diamond is non-symmorphic — it has glide planes — so some ops carry a non-zero `τ`. Pick out a representative one geometrically rather than by index (the post-identity order is not part of the API contract):

```jldoctest setup
julia> using LinearAlgebra

julia> ops[findfirst(op -> op.R == -Matrix{Int}(I, 3, 3), ops)]   # the unique op with R = -I
SpacegroupOp(R = [-1 0 0; 0 -1 0; 0 0 -1], τ = [0.25, 0.25, 0.25])
```

Geometrically, `R = -I` is pure inversion; the non-zero `τ = (¼, ¼, ¼)` shifts the inversion center off the origin. The fixed point is `c` solving `-c + τ = c`, i.e. `c = τ/2 = (⅛, ⅛, ⅛)` in lattice coordinates — the midpoint of the C–C bond between the two atoms in the primitive cell. That bond midpoint is the actual inversion center of diamond; with our origin choice (a C atom), the symmetry expresses itself as inversion-plus-shift rather than inversion-at-origin.

You can also construct a symmetry operation directly. The constructor canonicalizes the fractional shift `τ` to the interval `[0, 1)` automatically — `1.0`, `2.5`, and `-0.5` all fold to a value in that interval:

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

The identity is always `(R = I, τ = 0)`, even in non-symmorphic groups. By definition the identity satisfies `e * g == g` for every `g` in the group, which forces `R_e = I` and `τ_e = 0` (mod the lattice). Non-symmorphic-ness shows up in the *other* operations carrying non-zero `τ` (screw rotations, glide reflections), not in the identity. Choosing a different origin shifts the `τ` of those non-symmorphic ops — for diamond with the origin at the C–C bond midpoint, the inversion above would land at `τ = 0` and the diamond glides would pick up new shifts — but the identity itself stays `(I, 0)`.

## 2. Apply an op to a fractional position (and to a whole crystal)

The struct is callable. `op(r)` returns the image of fractional position `r`, folded back to `[0, 1)`. Pick a generic interior point so the rotational part of `op` produces a visibly different image (the origin is invariant under any pure rotation, so it would hide what `R` does):

```jldoctest
julia> using Spacey

julia> R = [0 -1 0; 1 0 0; 0 0 1];   # 4-fold about z, lattice coords

julia> op = SpacegroupOp(R, [0.5, 0.0, 0.0]);   # add a ½-translation along x

julia> op([0.25, 0.10, 0.30])
3-element Vector{Float64}:
 0.4
 0.25
 0.3
```

The 4-fold rotation sends `(x, y, z)` to `(-y, x, z)`, taking `(0.25, 0.10, 0.30)` to `(-0.10, 0.25, 0.30)`; the translation `(0.5, 0, 0)` then bumps the x-component by ½, giving `(0.40, 0.25, 0.30)`. No mod-1 folding was needed here — every component was already in `[0, 1)`.

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

julia> R1 = [0 -1 0; 1 0 0; 0 0 1];   # 4-fold rotation about z

julia> R2 = [1 0 0; 0 -1 0; 0 0 1];   # mirror in the xz-plane

julia> op1 = SpacegroupOp(R1, [0.5, 0.0, 0.0]);   # 4-fold + ½ x-shift

julia> op2 = SpacegroupOp(R2, [0.0, 0.0, 0.0]);   # pure mirror, no shift

julia> op12 = op1 * op2;

julia> op12.R == R1 * R2     # rotational parts multiply (op2 first, then op1)
true

julia> op12.R                # the composed operation
3×3 Matrix{Int64}:
 0  1  0
 1  0  0
 0  0  1
```

Geometrically the composite is a mirror in the diagonal `y = x` plane: the xz-mirror first sends `(x, y, z) → (x, -y, z)`, and the 4-fold then sends that to `(y, x, z)` — swapping the x and y axes, exactly the action of `op12.R`. The translation part of the composite, `R₁·τ₂ + τ₁`, is just `(½, 0, 0)` here because `τ₂ = 0`.

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

`inv(op)` returns the inverse operation that undoes `op`. Algebraically, `(R, τ)⁻¹ = (R⁻¹, -R⁻¹·τ)`. `R⁻¹` is integer because `|det R| = 1` for any lattice rotation, so the inverse computation is exact (no rounding error); the constructor enforces this and errors otherwise.

```jldoctest
julia> using Spacey

julia> R = [0 -1 0; 1 0 0; 0 0 1];

julia> op = SpacegroupOp(R, [0.25, 0.0, 0.0]);

julia> inv(op)                                # the inverse: -90° rotation, opposite shift
SpacegroupOp(R = [0 1 0; -1 0 0; 0 0 1], τ = [0.0, 0.25, 0.0])

julia> op * inv(op) == one(SpacegroupOp)      # gives the identity
true

julia> inv(inv(op)) == op                     # inverse applied twice returns the original op
true
```

Reading off the output: `inv(op).R = [0 1 0; -1 0 0; 0 0 1]` is the −90° rotation about z (the inverse of a +90° rotation), and `inv(op).τ = (0, ¼, 0)` is exactly `−R⁻¹·τ = −[0 1 0; -1 0 0; 0 0 1]·[¼, 0, 0] = (0, ¼, 0)`.

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

`spacegroup` returns ops in the user's lattice coordinates. To get the Cartesian rotation `R_c` and translation `τ_c` (e.g., for plotting or for combining with other Cartesian transforms), use [`toCartesian`](../reference/space-groups.md):

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

The result is a `Tuple{Matrix{Float64},Vector{Float64}}`, **not** a `SpacegroupOp` — `SpacegroupOp.R` must be `Matrix{Int}`, but the Cartesian rotation is not integer-valued in general. For an orthonormal cubic lattice (`A = I`), the Cartesian rotation `R_cart = R` literally. For a non-orthonormal lattice (e.g., hexagonal, where the basis vectors are not mutually perpendicular), `R_cart = A·R·A⁻¹` is the same operation expressed in the orthogonal Cartesian frame — its entries depend on the specific basis even though the underlying symmetry is intrinsic to the lattice.

## See also

- Reference: [`SpacegroupOp`](../reference/space-groups.md), [`toCartesian`](../reference/space-groups.md)
- How-to: [Find a space group](find-spacegroup.md)
- Explanation: [Canonicalizing τ](../explanation/canonicalizing-tau.md)
