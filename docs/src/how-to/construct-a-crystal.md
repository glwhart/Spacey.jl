# Construct a crystal

A `Crystal` bundles: a lattice, atomic positions, and per-atom type labels. The `coords` keyword is **required** at construction — there is no default — to prevent silent wrong-answer errors from misinterpreting, for example, `[0.5, 0.5, 0.5]` as Cartesian when it was meant to be fractional, or vice versa.

## 1. Pick a lattice

Either a 3×3 matrix (columns are basis vectors)…

```jldoctest
julia> using Spacey, LinearAlgebra

julia> A = Matrix{Float64}(I, 3, 3);   # cubic, side 1
```

…or three separate vectors:

```jldoctest
julia> using Spacey

julia> a1 = [1.0, 0, 0]; a2 = [0, 1.0, 0]; a3 = [0, 0, 1.0];
```

## 2. Place atoms

Positions are a `3 × N` matrix — one column per atom.

```jldoctest cscl
julia> using Spacey, LinearAlgebra

julia> A = Matrix{Float64}(I, 3, 3);

julia> r = [0.0 0.5; 0.0 0.5; 0.0 0.5];   # CsCl: Cs at origin, Cl at body-center
```

## 3. Pick type labels

Any vector of length `N`. Common choices: `Symbol`, `Int`, `String`. The type parameter is inferred from this vector — `Crystal{Symbol}`, `Crystal{Int}`, etc. Labels are used by `spacegroup` and `isSpacegroupOp` only to check that an op maps each atom to one of the *same type*.

```jldoctest cscl
julia> types = [:Cs, :Cl];
```

## 4. Construct

State explicitly whether the positions are fractional or Cartesian. Cartesian input is converted to fractional once at construction; internally, positions are always fractional.

**Matrix form, fractional input:**

```jldoctest cscl
julia> c = Crystal(A, r, types; coords=:fractional);

julia> length(spacegroup(c))
48
```

**Three-vector form, Cartesian input** (L1₀ CuAu — fcc-derived tetragonal primitive, with Cu at the origin and Au at the body center along z):

```jldoctest
julia> using Spacey

julia> a1 = [0.5,  0.5, 0.0];      # face-center vector, +y component

julia> a2 = [0.5, -0.5, 0.0];      # face-center vector, -y component

julia> a3 = [0.0,  0.0, 1.0];      # along z

julia> r_cart = [0.0  0.0;          # Cu at (0,0,0), Au at (0, ½, ½) in Cartesian
                 0.0  0.5;
                 0.0  0.5];

julia> c = Crystal(a1, a2, a3, r_cart, [:Cu, :Au]; coords=:cartesian);

julia> fractional(c)               # constructor converts Cartesian → fractional
3×2 Matrix{Float64}:
 0.0   0.5
 0.0  -0.5
 0.0   0.5
```

(Au's stored fractional position has a negative component because the constructor converts Cartesian → fractional via `inv(A) * r_cart` without folding mod 1. Translating Au by the lattice vector `a₁ - a₂ = (0, 1, 0)` gives the equivalent position `(0.5, 0.5, 0.5)`, but Spacey keeps the conversion result verbatim.)

## Common errors

- **`coords` keyword not used** → `MethodError`. Add `coords=:fractional` or `coords=:cartesian`.
- **`r` has wrong number of columns vs. atomic `types`** → constructor error with both lengths shown.
- **Singular `A`** (`det(A) ≈ 0`) → constructor rejects, with the error message that determinant is zero.
- **Empty crystal** (zero atoms) → not supported; constructor errors.
## See also

- Reference: [`Crystal`](../reference/crystals.md), [`fractional`](../reference/crystals.md), [`cartesian`](../reference/crystals.md)
- How-to: [Find a space group](find-spacegroup.md), [Classify the Bravais system](classify-bravais.md)
