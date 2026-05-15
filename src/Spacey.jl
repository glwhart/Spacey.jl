module Spacey
using MinkowskiReduction
using LinearAlgebra
using StatsBase
export pointgroup, snap_to_symmetry_svd, is_group,
       Crystal, is_spacegroup_op, fractional, cartesian, default_pos_tol,
       crystal_system, SpacegroupOp, to_cartesian, spacegroup,
       is_equiv_lattice, is_derivative,
       is_primitive, make_primitive,
       read_poscar
# Internal / not-exported (reach via `Spacey.<name>(...)`):
# - `pointgroup_robust`, `pointgroup_fast`, `pointgroup_simple` — the
#   public point-group entry point is `pointgroup` (which delegates to
#   `_robust`); `_fast` and `_simple` are validation/clean-input variants.
# - `aspect_ratio`, `rotate_basis_3d` — diagnostic / test scaffolding.
# - `snap_to_symmetry_avg` — less-robust alternative to `snap_to_symmetry_svd`.

"""
    avg_vec_over_ops(vec, ops)

Apply each operator in `ops` to `vec` and return the average over those images
that lie within 10% of the input's norm. Used internally by
[`snap_to_symmetry_avg`](@ref).

This is an internal helper; not exported.
"""
function avg_vec_over_ops(vec, ops)
    cands = [iop * vec for iop ∈ ops]
    cands = filter(x -> norm(x - vec) < 0.1 * norm(vec), cands)
    return sum(cands) / length(cands)
end

"""
    Spacey.snap_to_symmetry_avg(v1, v2, v3, ops)

Snap three basis vectors `v1, v2, v3` to a higher-symmetry triple by averaging
each vector over the images produced by `ops` (a vector of 3×3 lattice
operations). For each input vector, only images within 10% of its norm
contribute to the average — this filters the operations whose action
should fix that vector.

Internal helper — not exported. Prefer [`snap_to_symmetry_svd`](@ref), which
uses singular value decomposition of the metric tensor and is generally
more robust at high distortion. `snap_to_symmetry_avg` is kept as a faster
but less-robust alternative; reach it as `Spacey.snap_to_symmetry_avg(...)`.
Returns a tuple `(w1, w2, w3)`.
"""
function snap_to_symmetry_avg(v1, v2, v3, ops)
    w1 = avg_vec_over_ops(v1, ops)
    w2 = avg_vec_over_ops(v2, ops)
    w3 = avg_vec_over_ops(v3, ops)
    return w1, w2, w3
end

"""
    Spacey.snap_to_symmetry_avg(M, ops)

Matrix-form wrapper around `Spacey.snap_to_symmetry_avg(v1, v2, v3, ops)`:
treats the columns of `M` as the three basis vectors and returns the
snapped vectors as a 3×3 matrix. Internal helper — not exported.
"""
function snap_to_symmetry_avg(M, ops)
    w1, w2, w3 = snap_to_symmetry_avg(eachcol(M)..., ops)
    return [w1 w2 w3]
end

"""
    is_group(members)

Return `true` if `members` (a vector of square matrices) is closed under matrix
multiplication and contains no duplicates — i.e. forms a group.

Two methods are provided:
- For integer matrices, equality is exact.
- For floating-point matrices, equality uses `isapprox` with `atol`/`rtol`
  keyword arguments (default `1e-8` each).

Identity and inverses are not separately checked — they're implied by
finite closure of distinct elements (Cayley's theorem applied to the
finite case).

# Examples
```jldoctest
julia> using LinearAlgebra

julia> is_group([Matrix{Int}(I, 2, 2), -Matrix{Int}(I, 2, 2)])
true

julia> is_group([[0 1; 1 0]])         # not closed: M·M = I, which isn't in the set
false
```
"""
function is_group(members::AbstractVector{<:AbstractMatrix{<:Integer}})
    # Integer matrices: `==` and `hash` are exact element-wise, so `unique`
    # and `in` are O(n) via hashing.

    # 1) distinctness
    if length(unique(members)) < length(members)
        return false
    end

    # 2) closure (integer matrices: `C in members` uses exact `==`)
    for A in members, B in members
        C = A * B
        if !(C in members)
            return false
        end
    end

    return true
end


"""
    is_group(members::AbstractVector{<:AbstractMatrix{<:AbstractFloat}};
             atol=1e-8, rtol=1e-8)

Floating-point variant of [`is_group`](@ref): uses `isapprox(...; atol, rtol)`
for distinctness and closure checks. See the integer-matrix method for the
overall contract.
"""
function is_group(members::AbstractVector{<:AbstractMatrix{<:AbstractFloat}};
                  atol = 1e-8, rtol = 1e-8)
    # Float matrices: `isapprox` has no consistent hash, so we fall back to
    # O(n²) pairwise comparison. This is why the integer and float methods
    # diverge structurally.
    cmp(A, B) = isapprox(A, B; atol=atol, rtol=rtol)
    in_list(M, lst) = any(cmp(M, N) for N in lst)

    # 1) distinctness (approximate)
    for (k, A) in enumerate(members), B in @view members[(k + 1):end]
        if cmp(A, B)
            return false
        end
    end

    # 2) closure (approximate)
    for A in members, B in members
        C = A * B
        if !in_list(C, members)
            return false
        end
    end

    return true
end

"""
    is_equiv_lattice(A, B; tol=1e-6) -> Bool

Test whether two 3×3 lattice bases `A` and `B` describe the same lattice
— i.e., span the same set of points in space. Two bases are equivalent
iff the change-of-basis matrix `S = inv(A) * B` is unimodular (integer
entries, `|det(S)| = 1`).

This is a purely geometric test on the basis vectors; it does not look
at atoms, so it operates on bare matrices rather than [`Crystal`](@ref).

`tol` is an absolute tolerance on `S`'s deviation from the nearest
integer matrix and on `|det(S)|`'s deviation from 1.

# Examples
```jldoctest
julia> using Spacey, LinearAlgebra

julia> A = Matrix{Float64}(I, 3, 3);

julia> is_equiv_lattice(A, A * [1 1 0; 0 1 0; 0 0 1])    # unimodular shear → equivalent
true

julia> is_equiv_lattice(A, 2 * A)                         # different volume → not equivalent
false
```
"""
function is_equiv_lattice(A::AbstractMatrix, B::AbstractMatrix; tol::Real=1e-6)
    S = inv(A) * B
    return isapprox(abs(det(S)), 1.0; atol=tol) &&
           isapprox(S, round.(S); atol=tol)
end

"""
    is_derivative(parent, child; tol=1e-6) -> Bool

Test whether the lattice spanned by `child` is a **sublattice of** the lattice
spanned by `parent` — i.e. every `child`-lattice point is also a
`parent`-lattice point. The argument order is significant: this is a directed
question, *not* symmetric. `is_derivative(parent, child)` and
`is_derivative(child, parent)` answer different questions, and the second is
true only when `parent` and `child` span the same lattice (index 1). For the
symmetric "do they span the same lattice" question, use [`is_equiv_lattice`](@ref).

A child lattice is a derivative iff `S = inv(parent) * child` has integer
entries; `|det(S)|` is then the *index* of the sublattice (8 for a 2×2×2
cubic supercell, etc.). No volume constraint is imposed — that's the
distinction from `is_equiv_lattice`, which additionally requires `|det(S)| = 1`.

The standard use case is validating a candidate supercell against a known
parent: "given the primitive cell I'm enumerating from, is this HNF a valid
superlattice?" Operates on bare basis matrices, no atoms.

# Examples
```jldoctest
julia> using Spacey, LinearAlgebra

julia> parent = Matrix{Float64}(I, 3, 3);

julia> super = parent * [2 0 0; 0 2 0; 0 0 2];   # cubic 8× supercell

julia> is_derivative(parent, super)              # super ⊂ parent
true

julia> is_derivative(super, parent)              # parent ⊄ super (would need 1/2-integer entries)
false
```
"""
function is_derivative(parent::AbstractMatrix, child::AbstractMatrix; tol::Real=1e-6)
    S = inv(parent) * child
    return isapprox(S, round.(S); atol=tol)
end

"""
    Crystal{T}

A crystal structure: lattice vectors, atomic positions (fractional), and atom
type labels. `T` is the user's choice for type labels — typically `Int`,
`Symbol`, or `String`.

Fields:
- `A::Matrix{Float64}` — 3×3 lattice, columns are `a1`, `a2`, `a3`.
- `r::Matrix{Float64}` — 3×N positions in fractional coordinates, columns = atoms.
- `types::Vector{T}` — length N, atom type label per column of `r`.

Constructors:

    Crystal(A, r, types; coords)
    Crystal(a1, a2, a3, r, types; coords)

`coords` must be `:fractional` or `:cartesian` — **no default**, to prevent
silent wrong-answer errors from misinterpreting position data. See
`spacegroup_plan.md` §4.1. The constructor converts Cartesian input to
fractional once at construction; internally positions are always fractional.

Stored positions are folded into `[0, 1)` per axis. Atoms that differ by a
lattice translation therefore land on the same stored representation
regardless of how the user supplied them (Cartesian or fractional, inside
or outside the unit cell).

Numeric input is accepted as any `AbstractMatrix{<:Real}` / `AbstractVector{<:Real}`
and converted to `Float64` once at construction.

# Examples
```jldoctest
julia> using LinearAlgebra

julia> A = Matrix{Float64}(I, 3, 3);   # cubic lattice, side 1

julia> r = [0.0 0.5; 0.0 0.5; 0.0 0.5];   # CsCl: Cs at origin, Cl at body center

julia> c = Crystal(A, r, [:Cs, :Cl]; coords=:fractional);

julia> length(spacegroup(c))
48
```
"""
struct Crystal{T}
    A::Matrix{Float64}
    r::Matrix{Float64}
    types::Vector{T}
    function Crystal{T}(A::AbstractMatrix{<:Real}, r::AbstractMatrix{<:Real},
                        types::AbstractVector{T}; coords::Symbol) where T
        size(A) == (3, 3) || throw(ArgumentError("A must be 3×3"))
        size(r, 1) == 3 ||
            throw(ArgumentError("r must have 3 rows (one per spatial dimension)"))
        size(r, 2) == length(types) || throw(ArgumentError(
            "r has $(size(r, 2)) columns but types has length $(length(types))"))
        size(r, 2) > 0 ||
            throw(ArgumentError("empty crystal (no atoms) is not supported"))
        coords ∈ (:fractional, :cartesian) || throw(ArgumentError(
            "coords must be :fractional or :cartesian (got $(repr(coords)))"))
        A64 = Float64.(A)
        # Reject (near-)singular lattices. Scale-invariant test: compare |det|
        # to eps · ‖A‖³, i.e. the precision at which det is distinguishable
        # from zero for a matrix of this scale.
        abs(det(A64)) > eps(Float64) * opnorm(A64)^3 || throw(ArgumentError(
            "A is singular or near-singular (det = $(det(A64)))"))
        # Fold every position to the canonical interval [0, 1) per axis. This
        # applies whether the user supplied :fractional or :cartesian — in both
        # cases the stored representation is the equivalent atom inside the
        # primary unit cell. Without this fold, the Cartesian → fractional
        # conversion can leave components outside [0, 1) (e.g. an atom at
        # Cartesian (0, ½, ½) converted through a non-orthogonal A may produce
        # a -0.5 component verbatim), which is mathematically equivalent but
        # confusing to read off.
        r64 = coords === :cartesian ? inv(A64) * Float64.(r) : Float64.(r)
        r64 = mod.(r64, 1.0)
        new{T}(A64, r64, collect(types))
    end
end

Crystal(A::AbstractMatrix, r::AbstractMatrix, types::AbstractVector{T}; coords) where T =
    Crystal{T}(A, r, types; coords=coords)

Crystal(a1::AbstractVector, a2::AbstractVector, a3::AbstractVector,
        r::AbstractMatrix, types::AbstractVector; coords) =
    Crystal(hcat(a1, a2, a3), r, types; coords=coords)

"""
    read_poscar(path::AbstractString) -> Crystal

Read a POSCAR file (VASP 4 or VASP 5+) from `path` and return a [`Crystal`](@ref).

POSCAR layout (line numbers refer to the standard ordering):
1. comment / title — discarded
2. scaling factor — a single positive number scales the lattice; a negative
   number is treated as a target volume and the lattice is rescaled to match
3. three lines of lattice vectors (rows in POSCAR; transposed to columns
   in the returned `Crystal`)
4. (VASP 5+ only) element symbols, e.g. `Na Cl`
5. atom counts per species, e.g. `4 4`
6. (optional) `Selective dynamics` — recognised and skipped
7. coordinate type — `Direct`/`d` for fractional, `Cartesian`/`c`/`k` for Cartesian
8. one position per atom (`x y z`, with any trailing per-atom flags ignored)

VASP 4 files (which omit the element-symbol line) are detected by the first
token on line 6 being numeric; the species are then labelled `:X1`, `:X2`,
… in order. Trailing velocity / lattice-velocity blocks after the position
list are ignored.

The returned `Crystal` has `Symbol`-typed atom labels.

# Extended help

## Examples

```jldoctest
julia> using Spacey

julia> path = tempname() * ".poscar";

julia> write(path, \"\"\"
       NaCl conventional
       1.0
       5.64 0.0  0.0
       0.0  5.64 0.0
       0.0  0.0  5.64
       Na Cl
       4 4
       Direct
       0.0 0.0 0.0
       0.5 0.5 0.0
       0.5 0.0 0.5
       0.0 0.5 0.5
       0.5 0.5 0.5
       0.0 0.0 0.5
       0.0 0.5 0.0
       0.5 0.0 0.0
       \"\"\");

julia> c = read_poscar(path);

julia> c.types
8-element Vector{Symbol}:
 :Na
 :Na
 :Na
 :Na
 :Cl
 :Cl
 :Cl
 :Cl

julia> length(spacegroup(c))   # NaCl conventional cell, Pm-3m centering
192
```
"""
function read_poscar(path::AbstractString)
    lines = readlines(path)
    length(lines) ≥ 8 || throw(ArgumentError(
        "read_poscar: file too short to be a valid POSCAR ($(length(lines)) lines)"))
    idx = 2   # skip comment line

    # Line 2: scaling factor (single scalar).
    scale = parse(Float64, first(split(strip(lines[idx]))))
    idx += 1

    # Lines 3–5: lattice vectors as rows.
    rows = Matrix{Float64}(undef, 3, 3)
    for r in 1:3
        toks = split(strip(lines[idx]))
        rows[r, :] = parse.(Float64, toks[1:3])
        idx += 1
    end
    A_raw = collect(transpose(rows))   # POSCAR rows → Crystal columns
    A = if scale > 0
        scale * A_raw
    else
        target_volume = -scale
        s = (target_volume / abs(det(A_raw)))^(1/3)
        s * A_raw
    end

    # Line 6: either species symbols (VASP 5+) or counts (VASP 4).
    line6_toks = split(strip(lines[idx]))
    idx += 1
    local species, counts
    if all(isdigit, first(line6_toks))
        # VASP 4: counts only, no species names.
        counts = parse.(Int, line6_toks)
        species = [Symbol("X", i) for i in 1:length(counts)]
    else
        species = Symbol.(line6_toks)
        counts = parse.(Int, split(strip(lines[idx])))
        idx += 1
    end
    length(species) == length(counts) || throw(ArgumentError(
        "read_poscar: $(length(species)) species but $(length(counts)) counts"))

    # Optional "Selective dynamics" line — recognised by leading 'S' / 's'.
    if uppercase(first(strip(lines[idx]))) == 'S'
        idx += 1
    end

    # Coordinate type: 'D' → fractional, 'C'/'K' → Cartesian.
    coord_char = uppercase(first(strip(lines[idx])))
    idx += 1
    coords = if coord_char == 'D'
        :fractional
    elseif coord_char == 'C' || coord_char == 'K'
        :cartesian
    else
        throw(ArgumentError(
            "read_poscar: unrecognized coordinate type on line $(idx - 1): \"$(strip(lines[idx - 1]))\""))
    end

    # Positions: N atoms, each line is `x y z [trailing flags ignored]`.
    N = sum(counts)
    r = Matrix{Float64}(undef, 3, N)
    for i in 1:N
        toks = split(strip(lines[idx]))
        r[:, i] = parse.(Float64, toks[1:3])
        idx += 1
    end

    types = Symbol[]
    for (sym, n) in zip(species, counts)
        append!(types, fill(sym, n))
    end

    return Crystal(A, r, types; coords=coords)
end


"""
    fractional(c::Crystal)

Return the 3×N matrix of atomic positions in fractional (lattice) coordinates.
This is the canonical internal representation; see also [`cartesian`](@ref) for
the Cartesian view.

# Examples

When the crystal is built from Cartesian positions, `fractional` returns the
positions the constructor converted to and stored:

```jldoctest
julia> using LinearAlgebra

julia> A = Matrix{Float64}(2I, 3, 3);   # cubic lattice, edge length 2

julia> r_cart = reshape([1.0, 1.0, 1.0], 3, 1);   # atom at Cartesian (1, 1, 1)

julia> c = Crystal(A, r_cart, [:X]; coords=:cartesian);

julia> fractional(c)   # halfway along each basis vector
3×1 Matrix{Float64}:
 0.5
 0.5
 0.5
```
"""
fractional(c::Crystal) = c.r

"""
    cartesian(c::Crystal)

Return the 3×N matrix of atomic positions in Cartesian coordinates (same basis
and units as `c.A`). Each call recomputes from the stored fractional positions:
`A * r`. See also [`fractional`](@ref).

# Examples
```jldoctest
julia> using LinearAlgebra

julia> A = Matrix{Float64}(2I, 3, 3);

julia> r = reshape([0.5, 0.5, 0.5], 3, 1);

julia> c = Crystal(A, r, [:X]; coords=:fractional);

julia> cartesian(c)
3×1 Matrix{Float64}:
 1.0
 1.0
 1.0
```
"""
cartesian(c::Crystal) = c.A * c.r

"""
    default_pos_tol(c::Crystal)

Default position-matching tolerance used by [`is_spacegroup_op`](@ref) and
[`spacegroup`](@ref). Equal to `0.01 · (V/N)^(1/3)` where `V = |det(A)|` and
`N` is the number of atoms — 1% of the characteristic atom separation,
expressed in the same units as `c.A`. The formula is unit-agnostic, so it
returns a sensible default whether your lattice is in Ångström, Bohr, or
arbitrary units.

See `designDiscussions.md` for the rationale behind the 1% factor and the
class of structures it correctly classifies vs. the boundary cases it
silently over-promotes.

# Examples
```jldoctest
julia> using LinearAlgebra

julia> A = Matrix{Float64}(I, 3, 3);

julia> c = Crystal(A, reshape([0.0, 0.0, 0.0], 3, 1), [:X]; coords=:fractional);

julia> default_pos_tol(c)
0.01
```
"""
default_pos_tol(c::Crystal) = 0.01 * (abs(det(c.A)) / size(c.r, 2))^(1/3)

"""
    crystal_system(A; lattice_tol=0.01)
    crystal_system(c::Crystal; lattice_tol=0.01)

Return the Bravais-system symbol of the lattice `A` — one of
`:triclinic`, `:monoclinic`, `:orthorhombic`, `:tetragonal`, `:trigonal`,
`:hexagonal`, `:cubic`.

Identified via the order of the lattice's point group (its holohedry),
which uniquely determines the system:

    order →  system            holohedry
    ─────    ───────────       ─────────
      2      :triclinic        C_i   (-1)
      4      :monoclinic       C_2h  (2/m)
      8      :orthorhombic     D_2h  (mmm)
     12      :trigonal         D_3d  (-3m)
     16      :tetragonal       D_4h  (4/mmm)
     24      :hexagonal        D_6h  (6/mmm)
     48      :cubic            O_h   (m-3m)

Note: this reports the actual symmetry of the *lattice* Spacey sees.
If lattice parameters coincidentally match a higher-symmetry relation
(e.g. a ≈ b in an orthorhombic cell at default `lattice_tol`), the
returned system may be higher than the nominal one — same behavior as
[`pointgroup`](@ref).

# Examples
```jldoctest
julia> using LinearAlgebra

julia> crystal_system(Matrix{Float64}(I, 3, 3))
:cubic

julia> crystal_system([1.0 0 0; 0 1 0; 0 0 1.5])
:tetragonal

julia> crystal_system([1.0 0 0; 0 1.2 0; 0 0 1.5])
:orthorhombic
```
"""
function crystal_system(A::AbstractMatrix{<:Real}; lattice_tol::Real=0.01)
    A_red = minkReduce(Float64.(A))
    u, v, w = eachcol(A_red)
    LG = pointgroup_robust(u, v, w; tol=lattice_tol, auto_reduce=false)
    order = length(LG)
    order == 2  && return :triclinic
    order == 4  && return :monoclinic
    order == 8  && return :orthorhombic
    order == 12 && return :trigonal
    order == 16 && return :tetragonal
    order == 24 && return :hexagonal
    order == 48 && return :cubic
    error("unexpected lattice point-group order $order (should be 2, 4, 8, 12, 16, 24, or 48)")
end

crystal_system(c::Crystal; kwargs...) = crystal_system(c.A; kwargs...)

"""
    is_spacegroup_op(R, τ, c::Crystal; tol=default_pos_tol(c))

Return `true` if the operation `(R, τ)` is a space-group symmetry of crystal
`c` — that is, if applying `R` then translating by `τ` (all in fractional
coordinates) maps the set of atomic positions to itself, preserving types,
modulo the lattice. Returns `false` otherwise.

Each original atom's image must coincide with an (injectively matched)
original atom of the same type, with per-component distance below `tol` after
wrapping the signed difference into `(-½, ½]` (i.e. comparing modulo the
lattice).

# Examples
```jldoctest
julia> using LinearAlgebra

julia> A = Matrix{Float64}(I, 3, 3);

julia> c = Crystal(A, reshape([0.0, 0.0, 0.0], 3, 1), [:X]; coords=:fractional);

julia> I3 = Matrix{Int}(I, 3, 3);

julia> is_spacegroup_op(I3, [0.0, 0.0, 0.0], c; tol=1e-8)
true

julia> is_spacegroup_op(I3, [0.5, 0.0, 0.0], c; tol=1e-8)
false
```
"""
function is_spacegroup_op(R::AbstractMatrix{<:Real}, τ::AbstractVector{<:Real},
                          c::Crystal; tol::Real=default_pos_tol(c))
    size(R) == (3, 3) || throw(ArgumentError("R must be 3×3"))
    length(τ) == 3 || throw(ArgumentError("τ must have length 3"))
    N = size(c.r, 2)
    r_image = mod.(R * c.r .+ τ, 1.0)
    claimed = falses(N)
    for i in 1:N
        found = false
        for j in 1:N
            claimed[j] && continue
            c.types[i] == c.types[j] || continue
            Δ = mod.(r_image[:, i] .- c.r[:, j] .+ 0.5, 1.0) .- 0.5
            if all(abs.(Δ) .< tol)
                claimed[j] = true
                found = true
                break
            end
        end
        found || return false
    end
    return true
end

"""
    SpacegroupOp

A single space-group operation `r ↦ R·r + τ`, expressed in lattice (fractional)
coordinates: integer rotation `R`, fractional translation `τ`.

`τ` is canonicalized to `[0, 1)` at construction via `mod.(τ, 1.0)`, so
`SpacegroupOp(I, [0,0,0])`, `SpacegroupOp(I, [1,0,0])`, and
`SpacegroupOp(I, [2.5, 0, 0])` all produce the same stored representation.
This makes Julia's default field-by-field `==` and `hash` consistent with
the periodic-boundary semantics a user expects.

Returned by `spacegroup(c::Crystal)`. Compose with `*`, invert with `inv`,
apply to a fractional position via `op(r)`, convert to Cartesian with
`to_cartesian(op, A)`.

# Examples
```jldoctest
julia> using LinearAlgebra

julia> e = one(SpacegroupOp);   # identity

julia> R = [0 -1 0; 1 0 0; 0 0 1];   # 4-fold rotation about z (in lattice coords)

julia> op = SpacegroupOp(R, [0.5, 0.0, 0.0]);

julia> op([0.0, 0.0, 0.0])
3-element Vector{Float64}:
 0.5
 0.0
 0.0

julia> (op * inv(op)) == e
true
```
"""
struct SpacegroupOp
    R::Matrix{Int}
    τ::Vector{Float64}
    # Canonicalize τ: fold mod 1 into [0, 1), then snap each component to
    # the nearest p/q with q ≤ 12 if within 1e-6. Every τ component in an
    # ITA space group is an exact rational with small denominator (0, ½,
    # ⅓, ¼, ⅙, ¹/₁₂, …), so snapping preserves them exactly while killing
    # float drift accumulated through basis transforms and composition.
    # Earlier implementation used round(τ, digits=10), but that had a
    # silent bug: round(1/3, digits=10) = 0.3333333333 and
    # round(2/3, digits=10) = 0.6666666667, so 1/3 + 1/3 no longer
    # matched 2/3, breaking closure on trigonal groups (P3₁21 etc.).
    SpacegroupOp(R, τ) = new(R, _canonicalize_τ(τ))
end

function _canonicalize_τ(τ::AbstractVector; tol::Real=1e-6)
    out = Vector{Float64}(undef, length(τ))
    for i in eachindex(τ)
        x = mod(Float64(τ[i]), 1.0)
        snapped = x
        for q in 1:12
            p = round(Int, q * x)
            # x ≈ 1 wraps to 0 in [0, 1) semantics
            if p == q
                if abs(q * x - q) < q * tol
                    snapped = 0.0
                    break
                end
                continue
            end
            cand = p / q
            if abs(x - cand) < tol
                snapped = cand
                break
            end
        end
        out[i] = snapped
    end
    return out
end

# Composition: op1 * op2 means "apply op2 first, then op1" (function-composition
# semantics). Derivation:
#   op2: r ↦ R2·r + τ2
#   op1 applied to that: R1·(R2·r + τ2) + τ1 = R1·R2·r + R1·τ2 + τ1
Base.:*(a::SpacegroupOp, b::SpacegroupOp) =
    SpacegroupOp(a.R * b.R, a.R * b.τ + a.τ)

# Inverse: (R, τ)⁻¹ = (R⁻¹, -R⁻¹·τ). R⁻¹ is integer because |det R| = 1 for
# any lattice rotation. The rounding + sanity check guards against misuse
# with a non-lattice R.
function Base.inv(op::SpacegroupOp)
    Rinv_f = inv(Float64.(op.R))
    Rinv = round.(Int, Rinv_f)
    maximum(abs, Rinv_f .- Rinv) < 1e-8 || throw(ArgumentError(
        "inv(SpacegroupOp): R⁻¹ is not integer (det(R) ≠ ±1?)"))
    return SpacegroupOp(Rinv, -Rinv * op.τ)
end

# Apply to a fractional position vector (callable struct)
(op::SpacegroupOp)(r::AbstractVector) = mod.(op.R * r + op.τ, 1.0)

# Identity op
Base.one(::Type{SpacegroupOp}) = SpacegroupOp(Matrix{Int}(I, 3, 3), zeros(3))

# Equality and hash. Julia's default `==` for a struct with Vector/Matrix
# fields falls back to `===` (object identity), which would say two ops
# with identical content are unequal. We override with explicit field-by-
# field `==` (element-wise for R and τ). Because τ is canonicalized to
# [0, 1) at construction, this correctly treats ops with τ=[0,0,0] and
# τ=[1,0,0] as equal (both stored as [0,0,0]). The matching `hash` method
# keeps Set{SpacegroupOp} and Dict{SpacegroupOp,_} consistent with `==`.
Base.:(==)(a::SpacegroupOp, b::SpacegroupOp) = a.R == b.R && a.τ == b.τ
Base.hash(op::SpacegroupOp, h::UInt) = hash(op.τ, hash(op.R, hash(:SpacegroupOp, h)))

# Pretty printing (Julia calls this automatically for REPL, println, etc.)
Base.show(io::IO, op::SpacegroupOp) =
    print(io, "SpacegroupOp(R = ", op.R, ", τ = ", op.τ, ")")

"""
    to_cartesian(op::SpacegroupOp, A::AbstractMatrix)

Convert a lattice-coordinate space-group operation to its Cartesian form.
Returns the tuple `(R_cart, τ_cart) = (A·R·A⁻¹, A·τ)` where `A` is the lattice
matrix whose columns are the basis vectors.

The result is a `Tuple{Matrix{Float64},Vector{Float64}}` rather than a
`SpacegroupOp`, because the Cartesian rotation is in general not integer-
valued while [`SpacegroupOp`](@ref)'s `R` field must be `Matrix{Int}`.

# Examples
```jldoctest
julia> using LinearAlgebra

julia> A = Matrix{Float64}(I, 3, 3);

julia> R_cart, τ_cart = to_cartesian(one(SpacegroupOp), A);

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
"""
to_cartesian(op::SpacegroupOp, A::AbstractMatrix) =
    (A * op.R * inv(A), A * op.τ)

"""
    to_cartesian(op::AbstractMatrix{<:Integer}, A::AbstractMatrix)
    to_cartesian(LG::AbstractVector{<:AbstractMatrix{<:Integer}}, A::AbstractMatrix)

Convert a lattice-coordinate point-group operation (a single integer
matrix) — or a whole vector of them, like the result of [`pointgroup`](@ref) —
to Cartesian rotation form. Returns `A · op · inv(A)` for the single-op
method and `[A · op · inv(A) for op in LG]` for the vector method.

These overloads exist because `pointgroup` returns just the integer-matrix
form (since v0.8); use these helpers if you need the Cartesian rotations.

# Examples
```jldoctest
julia> using Spacey, LinearAlgebra

julia> A = Matrix{Float64}(I, 3, 3);

julia> LG = pointgroup(A);

julia> G = to_cartesian(LG, A);   # Cartesian rotations of the cubic point group

julia> length(G)
48

julia> G[findfirst(==(Matrix{Int}(I, 3, 3)), LG)]   # identity in Cartesian
3×3 Matrix{Float64}:
 1.0  0.0  0.0
 0.0  1.0  0.0
 0.0  0.0  1.0
```
"""
to_cartesian(op::AbstractMatrix{<:Integer}, A::AbstractMatrix) = A * op * inv(A)

to_cartesian(LG::AbstractVector{<:AbstractMatrix{<:Integer}}, A::AbstractMatrix) =
    [A * op * inv(A) for op in LG]

"""
    Spacey.rotate_basis_3d(u, v, w, α, β, γ)

Rotate the basis vectors `u, v, w` by Euler angles `α, β, γ` (the
yaw–pitch–roll convention used in test scaffolding). Returns the rotated
triple `(u', v', w')` as a tuple of three vectors.

Internal helper — not exported. Used by tests to verify that
symmetry-finding routines are invariant under arbitrary lattice
orientation. Reach as `Spacey.rotate_basis_3d(...)`.

The rotation matrix is built from successive rotations about the z, y, and
z axes (matching the order in the formula). For zero angles the identity
is returned.

# Examples
```jldoctest
julia> u, v, w = Spacey.rotate_basis_3d([1.0,0,0], [0,1.0,0], [0,0,1.0], 0.0, 0.0, 0.0);

julia> u
3-element Vector{Float64}:
 1.0
 0.0
 0.0
```
"""
function rotate_basis_3d(u, v, w, α, β, γ)
    A = [u v w]
    R = [cos(α)cos(β)  cos(α)sin(β)sin(γ) - sin(α)cos(γ)  cos(α)sin(β)cos(γ) + sin(α)sin(γ);
         sin(α)cos(β)  sin(α)sin(β)sin(γ) + cos(α)cos(γ)  sin(α)sin(β)cos(γ) - cos(α)sin(γ);
         -sin(β)       cos(β)sin(γ)                       cos(β)cos(γ)]
    A = R * A
    return A[:, 1], A[:, 2], A[:, 3]
end

"""
    Spacey.pointgroup_simple(a1, a2, a3; debug=false)

Brute-force enumeration of the point group of a 3D lattice. Iterates over
every 3×3 candidate matrix with entries in `{-1, 0, 1}` (3⁹ = 19683
matrices), retains those whose action on the basis preserves the metric
tensor, and returns the survivors as Cartesian rotations.

Internal — not exported. The simplest correct implementation; used to
validate the more efficient `Spacey.pointgroup_fast` and the public
[`pointgroup`](@ref) (which delegates to `Spacey.pointgroup_robust`).
It performs strict (`isapprox` with default tolerance) equality checks,
so it is most reliable on noiseless / synthetic input.

If `debug=true`, returns the candidate `T = UᵀU` matrices instead of the
filtered ops, for use when diagnosing failures.

# Examples
```jldoctest
julia> using Spacey

julia> u = [1.0, 0, 0]; v = [0.5, √3/2, 0]; w = [0.0, 0, √(8/3)];

julia> length(Spacey.pointgroup_simple(u, v, w))
24
```
"""
function pointgroup_simple(a1, a2, a3; debug::Bool=false)
    u, v, w = minkReduce(a1, a2, a3)
    A = [u v w]                       # Put the lattice vectors as columns in matrix A
    B = inv(A) * transpose(inv(A))    # Use this for checking for orthogonality
    # A list of all possible lattice vectors in a rotated basis
    c = [A * [i; j; k] for i ∈ (-1, 0, 1) for j ∈ (-1, 0, 1) for k ∈ (-1, 0, 1)]
    # A list of all possible bases, (i.e., all combinations of c vectors)
    R = [[x y z] for x ∈ c for y ∈ c for z ∈ c]
    RT = [transpose(M) for M ∈ R]
    # This is the U^T*U, where U transforms original basis to candidate basis
    T = [R[i] * B * RT[i] for i in eachindex(R)]
    if debug return T end
    # If T==identity then the U was a symmetry of the lattice
    idx = findall(t ≈ I(3) for t ∈ T)
    Ai = inv(A)
    return [Ai * R[i] for i in idx]
end


"""
    Spacey.pointgroup_fast(a1, a2, a3)

Production-speed point-group finder for an exact / noiseless 3D lattice.
Faster than `Spacey.pointgroup_simple` by filtering candidate basis
combinations by length and volume before checking orthogonality, but uses
strict `isapprox` tolerance and so is best suited to clean inputs.

Internal — not exported. For real-world (noisy) input use the public
[`pointgroup`](@ref), which exposes a tolerance keyword.

Returns operations as integer matrices in lattice coordinates.

# Examples
```jldoctest
julia> using Spacey

julia> u = [1.0, 0, 0]; v = [0.5, √3/2, 0]; w = [0.0, 0, √(8/3)];

julia> length(Spacey.pointgroup_fast(u, v, w))
24
```
"""
function pointgroup_fast(a1, a2, a3)
    u, v, w = minkReduce(a1, a2, a3)   # Always do this first, algorithm assumes reduced basis
    A = [u v w]                        # Define a matrix with input vectors as columns
    Ai = inv(A)
    AiAiT = Ai * transpose(Ai)         # Use this for checking for orthogonality
    norms = norm.([u, v, w])           # Compute the norms of the three input vectors
    vol = abs(u × v ⋅ w)               # Volume of the parallelipiped formed by the basis vectors

    # A list of all possible lattice vectors in a rotated basis
    # These are lattice points from the vertices of the 8 cells with a corner at the origin)
    # There are 27 of these (==3^3)
    c = [A * [i, j, k] for i ∈ (-1, 0, 1) for j ∈ (-1, 0, 1) for k ∈ (-1, 0, 1)]
    # Now keep only those vectors that have a norm matching one of the input vectors
    # efficiency: Gather three groups, according to length. This limits the candidates even more
    c1 = c[findall(norm(x) ≈ norms[1] for x ∈ c)]   # All vectors with first norm
    c2 = c[findall(norm(x) ≈ norms[2] for x ∈ c)]   # All vectors with second norm
    c3 = c[findall(norm(x) ≈ norms[3] for x ∈ c)]   # All vectors with third norm
    # Construct all possible bases, (i.e., all combinations of c vectors), skip duplicate vectors
    R = [[x y z] for x ∈ c1 for y ∈ c2 if !(x ≈ y) for z ∈ c3 if !(x ≈ z) && !(y ≈ z)]
    R = R[findall(abs(det(r)) ≈ vol for r in R)]    # Delete candidate bases with the wrong volume
    # The cross product is slightly (<1%) faster
    #R = R[findall([abs(r[1]×r[2]⋅r[3])≈vol for r in R])] # Delete candidate bases with the wrong volume
    RT = [transpose(M) for M ∈ R]
    # This is the Uᵀ ̇U, where U transforms original basis to candidate basis
    # If Tᵢ==identity then the U was a symmetry of the lattice
    T = [R[i] * AiAiT * RT[i] for i in eachindex(R)]
    # Indices of candidate T's that match the identity
    idx = findall(t ≈ I(3) for t ∈ T)
    # Convert the transformations to integer matrices (formally they should be)
    return [round.(Int, Ai * R[i]) for i in idx]
end

"""
    Spacey.pointgroup_robust(u, v, w; tol=0.01, verify_stable=false, auto_reduce=true)

Tolerance-tunable point-group finder for noisy real-world input. Returns
a `Vector{Matrix{Int}}` of the lattice-coordinate symmetry operations —
the same form the public [`pointgroup`](@ref) returns. For Cartesian
rotations, pass the result through [`to_cartesian`](@ref).

Internal — not exported. The public entry point [`pointgroup`](@ref) is a
thin wrapper around this function (with the same defaults). Reach this
form directly only when explicitly disambiguating between point-group
variants (e.g. comparing against `Spacey.pointgroup_fast`).

# Keyword arguments
- `tol::Real=0.01` — relative tolerance applied to the (volume-normalized)
  lattice. Tighter values reject more spurious candidates; looser values
  tolerate more input noise but risk over-promotion to higher symmetry.
- `verify_stable::Bool=false` — opt-in stability check. When `true`, the
  algorithm re-runs at `tol/1000` and emits a `@warn` if the operation
  count differs between the two runs (i.e. the lattice is near a
  symmetry boundary). The returned group is unchanged.
- `auto_reduce::Bool=true` — Minkowski-reduce the input automatically
  before searching for symmetries. The returned operations are still
  expressed in the user's *input* basis (a unimodular change-of-basis is
  applied internally to map the reduced-basis ops back). Set to `false`
  to assert the input is already reduced; the function then errors if
  it isn't (the pre-v0.8 strict behavior, useful as a self-check).

Algorithm: Minkowski-reduce input (or verify it's already reduced if
`auto_reduce=false`), enumerate candidate basis permutations from the
{-1,0,1}³ neighbor set, filter by norm match → volume conservation →
orthogonality, keep the largest subset that closes under multiplication,
then transform the result back to the user's basis if needed.

A `@warn` fires automatically when the input aspect ratio exceeds 100 —
results may be unreliable for ratios above ~500 due to floating-point
precision in the candidate-detection step.

# Examples
```jldoctest
julia> using Spacey

julia> u = [1.0, 0, 0]; v = [0, 1.0, 0]; w = [0, 0, 1.0];

julia> length(Spacey.pointgroup_robust(u, v, w))
48
```
"""
function pointgroup_robust(u, v, w; tol=0.01, verify_stable::Bool=false,
                                    auto_reduce::Bool=true)
    # Handle the input basis: either reduce it ourselves (default) or verify the
    # caller's reduced-basis assertion. When auto-reducing we keep the change-of-
    # basis matrix so we can map the operations back to the user's basis at the end.
    local U_basis::Matrix{Int}, Uinv_basis::Matrix{Int}
    needs_basis_change = false
    if auto_reduce
        A_orig = hcat(u, v, w)
        u, v, w = minkReduce(u, v, w)[1:3]
        A_red = hcat(u, v, w)
        # A_orig = A_red · U_basis where U_basis is unimodular (both bases span the
        # same lattice). Round to absorb floating-point noise.
        U_basis = round.(Int, inv(A_red) * A_orig)
        Uinv_basis = round.(Int, inv(Float64.(U_basis)))
        abs(det(U_basis)) == 1 ||
            error("Mink reduction yielded a non-unimodular change-of-basis (det = $(det(U_basis))).")
        needs_basis_change = U_basis != Matrix{Int}(I, 3, 3)
    else
        # Mink reduction can change the basis even when the basis is already reduced (degenerate cases). So don't do it here. But do check that no reduction is needed.
        if !(orthogonalityDefect(u, v, w) ≈ orthogonalityDefect(minkReduce(u, v, w)[1:3]...))
            throw(ArgumentError(
                "Input basis for 'pointgroup' is not Minkowski-reduced. Either pass `auto_reduce=true` (the default) or run `minkReduce` first."))
        end
    end
    input_vol = ∛(abs(u × v ⋅ w))    # Rescale the basis to have a volume of 1, avoid floating point issues
    u, v, w = u ./ input_vol, v ./ input_vol, w ./ input_vol

    norms = norm.([u, v, w])         # Compute the norms of the three input vectors
    ar = maximum(norms) / minimum(norms)
    if ar > 100
        @warn "Aspect ratio is $(round(ar, digits=1)). Results may be unreliable for ratios above ~500."
    end

    A = [u v w]                      # Define a matrix with input vectors as columns
    Ai = inv(A)
    vol = abs(u × v ⋅ w)             # Volume of the parallelipiped formed by the basis vectors

    # A list of all possible lattice vectors in a rotated basis. These are lattice points from the vertices of the 8 cells that have a corner at the origin. There are 27 of these (==3^3)
    c = [A * [i, j, k] for i ∈ (-1, 0, 1) for j ∈ (-1, 0, 1) for k ∈ (-1, 0, 1)]
    # Now keep only those vectors that have a norm close the norm one of the input vectors
    # efficiency: Gather three groups, according to length. This limits the candidates even more
    c1 = c[findall(isapprox(norms[1], norm(x), rtol=tol) for x ∈ c)]   # All vectors with first norm
    c2 = c[findall(isapprox(norms[2], norm(x), rtol=tol) for x ∈ c)]   # All vectors with second norm
    c3 = c[findall(isapprox(norms[3], norm(x), rtol=tol) for x ∈ c)]   # All vectors with third norm
    # Construct all candidate bases, Rc (i.e., all combinations of c vectors), skip duplicate vectors.
    A′ = [[x y z] for x ∈ c1 for y ∈ c2 if !(x ≈ y) for z ∈ c3 if !(x ≈ z) && !(y ≈ z)]   # All candidate bases
    A′ = A′[findall(isapprox(abs(det(M)), vol, rtol=tol * min(norms...)) for M in A′)]    # Delete candidate bases with the wrong volume
    Rc = [M * Ai for M ∈ A′]         # Compute the candidate rotations from the candidate bases

    # This is the Uᵀ ̇U, where U transforms original basis to candidate basis
    # If Tᵢ==identity then the Rc is orthogonal and is a symmetry of the lattice
    T = [transpose(rc) * rc for rc ∈ Rc]
    # Indices of candidate T's that match the identity
    idx = findall(isapprox(t, I(3), rtol=tol) for t ∈ T)
    Rc = Rc[idx]
    T = T[idx]
    # Convert the transformations to lattice coordinates representation (round to integer matrices; formally they should be)
    ops = [round.(Int, Ai * M * A) for M in Rc]   # Need the 'Int' so integers are returned
    # Get norms of deviation from orthogonal case
    tn = [norm(t - I(3)) for t ∈ T]
    tp = sortperm(tn)                # Sort by deviation
    # Find the largest number of (sorted) ops that form a group.
    best_order = 48
    for n ∈ [48, 24, 16, 12, 8, 4, 2]   # These are the only possible group sizes for a 3D lattice
        if n > length(idx) continue end
        if is_group(ops[tp[1:n]])    # Keep the largest set that is a group
            best_order = n
            break
        end
    end
    result_ops = ops[tp][1:best_order]
    if needs_basis_change
        # Map ops from the reduced-basis representation back to the user's basis:
        # if R_cart = A_red · M_red · inv(A_red) = A_orig · M_orig · inv(A_orig)
        # and A_orig = A_red · U_basis, then M_orig = U_basis⁻¹ · M_red · U_basis.
        result_ops = [Uinv_basis * op * U_basis for op in result_ops]
    end
    if verify_stable
        tight_tol = tol / 1000
        # Skip auto_reduce in the recursion: at this point the local u,v,w are
        # already reduced (and rescaled), and verify_stable only inspects the
        # length of the returned group — which is invariant under basis choice.
        tight_ops = pointgroup_robust(u, v, w; tol=tight_tol, verify_stable=false,
                                                auto_reduce=false)
        if length(tight_ops) != length(result_ops)
            @warn "pointgroup_robust: group size depends on tolerance — lattice is near a symmetry boundary." tol group_at_tol=length(result_ops) tight_tol group_at_tight_tol=length(tight_ops)
        end
    end
    return result_ops
end

"""
    Spacey._probe_atoms(c::Crystal) -> (probe_indices, i0)

Pick the atom type with the fewest representatives in `c` as the probe
type. Returns the indices of all atoms of that type plus the first such
index. Used by [`spacegroup`](@ref) and [`_find_self_translations`](@ref)
to keep the per-rotation τ-enumeration small.

Internal — not exported.
"""
function _probe_atoms(c::Crystal)
    types_unique = unique(c.types)
    counts = [count(==(t), c.types) for t in types_unique]
    probe_type = types_unique[argmin(counts)]
    probe_indices = findall(==(probe_type), c.types)
    return probe_indices, probe_indices[1]
end

"""
    Spacey._find_translations_for_rotation(R, c, probe_indices, i0; pos_tol) -> Vector{Vector{Float64}}

For a candidate rotation `R` (integer matrix in `c.A`'s basis), enumerate
the fractional translations `τ` such that the operation `(R, τ)` is a
symmetry of `c`. Each surviving τ has been verified by [`is_spacegroup_op`](@ref).

This is the inner loop of [`spacegroup`](@ref) and the core of
[`_find_self_translations`](@ref) (which fixes `R = I`).

Internal — not exported.
"""
function _find_translations_for_rotation(R, c::Crystal, probe_indices, i0;
                                         pos_tol::Real)
    image_i0 = R * c.r[:, i0]
    τs = Vector{Float64}[]
    for j in probe_indices
        τ = mod.(c.r[:, j] .- image_i0, 1.0)
        if is_spacegroup_op(R, τ, c; tol=pos_tol)
            push!(τs, τ)
        end
    end
    return τs
end

"""
    Spacey._find_self_translations(c::Crystal; pos_tol=default_pos_tol(c))
        -> Vector{Vector{Float64}}

Return the fractional translations `τ` such that `(R = I, τ)` is a symmetry
of crystal `c`. The list always contains `τ = (0, 0, 0)` (the identity);
extra entries are the centering translations of a non-primitive cell.

`length(_find_self_translations(c)) == 1` iff `c` is primitive — see
[`is_primitive`](@ref) and [`make_primitive`](@ref).

Internal — not exported.
"""
function _find_self_translations(c::Crystal; pos_tol::Real=default_pos_tol(c))
    probe_indices, i0 = _probe_atoms(c)
    R = Matrix{Int}(I, 3, 3)
    return _find_translations_for_rotation(R, c, probe_indices, i0; pos_tol)
end

"""
    is_primitive(c::Crystal; pos_tol=default_pos_tol(c)) -> Bool

Return `true` iff the crystal `c` is described in a primitive cell — i.e.
the only fractional translation that maps `c` onto itself is `(0, 0, 0)`.
A `false` return means `c` has a centering translation; pass it to
[`make_primitive`](@ref) to reduce.

The criterion is purely geometric (no symmetry operations applied) and
runs in `O(N²)` for `N` atoms.

# Examples
```jldoctest
julia> using Spacey, LinearAlgebra

julia> A = Matrix{Float64}(I, 3, 3);

julia> c1 = Crystal(A, reshape([0.0, 0, 0], 3, 1), [:X]; coords=:fractional);   # primitive cubic

julia> is_primitive(c1)
true

julia> c2 = Crystal(A, [0.0 0.5; 0.0 0.5; 0.0 0.5], [:X, :X]; coords=:fractional);  # BCC, conventional cell

julia> is_primitive(c2)
false
```
"""
is_primitive(c::Crystal; pos_tol::Real=default_pos_tol(c)) =
    length(_find_self_translations(c; pos_tol)) == 1

"""
    make_primitive(c::Crystal; pos_tol=default_pos_tol(c)) -> (Crystal, Vector{Int})

Return a primitive description of crystal `c`, plus the indices of atoms
that were dropped (one representative kept per centering equivalence class).
If `c` is already primitive, returns `(c, Int[])`.

Algorithm: find the fractional self-translations of `c` (the centering
group). Build a candidate Cartesian basis pool from the original lattice
columns plus the non-zero centering vectors, then search triples for one
whose volume is `|det(c.A)| / k` (where `k` is the centering multiplicity)
and which expresses every candidate vector as an integer combination —
that triple is a primitive basis. Atomic positions are folded into the
primitive cell and deduplicated.

# Examples
```jldoctest
julia> using Spacey, LinearAlgebra

julia> A = Matrix{Float64}(I, 3, 3);

julia> c_bcc = Crystal(A, [0.0 0.5; 0.0 0.5; 0.0 0.5], [:X, :X]; coords=:fractional);

julia> c_prim, removed = make_primitive(c_bcc);

julia> size(c_prim.r, 2)            # one atom per primitive cell
1

julia> isapprox(abs(det(c_prim.A)), abs(det(A)) / 2)   # half the conventional volume
true

julia> length(removed)              # one of the two original atoms was dropped
1
```
"""
function make_primitive(c::Crystal; pos_tol::Real=default_pos_tol(c))
    T = _find_self_translations(c; pos_tol)
    k = length(T)
    if k == 1
        return c, Int[]
    end

    # Candidate primitive lattice vectors — all guaranteed to lie in the
    # primitive lattice. Original cell columns expressed in fractional coords
    # are the standard basis vectors; non-zero centering translations are
    # already fractional. Convert all to Cartesian for the volume search.
    cands_frac = Vector{Float64}[[1.0, 0, 0], [0, 1.0, 0], [0, 0, 1.0]]
    for τ in T
        if norm(τ) > pos_tol     # skip the identity translation
            push!(cands_frac, collect(τ))
        end
    end
    cands_cart = [c.A * f for f in cands_frac]
    perm = sortperm(norm.(cands_cart))   # prefer short bases
    cands_cart = cands_cart[perm]
    n = length(cands_cart)

    V_orig = abs(det(c.A))
    V_target = V_orig / k

    A_new = nothing
    for i in 1:n-2, j in i+1:n-1, l in j+1:n
        B = hcat(cands_cart[i], cands_cart[j], cands_cart[l])
        d = abs(det(B))
        if d < 1e-12 * V_orig
            continue   # degenerate triple
        end
        if !isapprox(d, V_target; rtol=1e-6)
            continue
        end
        # Verify every candidate is an integer combination of B's columns.
        Binv = inv(B)
        all_integer = true
        for v in cands_cart
            f = Binv * v
            if !isapprox(f, round.(f); atol=1e-6)
                all_integer = false
                break
            end
        end
        if all_integer
            A_new = B
            break    # first valid triple — shortest by sort order
        end
    end

    A_new === nothing &&
        error("make_primitive: failed to find a primitive basis from the candidate set (this is a bug — please open an issue with the input crystal).")

    # Fold atomic positions into the new (smaller) primitive cell, then
    # deduplicate. The centering translations guarantee N_new = N / k.
    r_new = mod.(inv(A_new) * c.A * c.r, 1.0)
    keep = Int[]
    removed = Int[]
    for i in 1:size(r_new, 2)
        is_dup = false
        for j in keep
            c.types[i] == c.types[j] || continue
            # Signed mod-1 difference, then back to Cartesian for tolerance check.
            diff_frac = mod.(r_new[:, i] .- r_new[:, j] .+ 0.5, 1.0) .- 0.5
            if norm(A_new * diff_frac) < pos_tol
                is_dup = true
                break
            end
        end
        if is_dup
            push!(removed, i)
        else
            push!(keep, i)
        end
    end

    return Crystal(A_new, r_new[:, keep], c.types[keep]; coords=:fractional), removed
end

"""
    spacegroup(c::Crystal; lattice_tol=0.01, pos_tol=default_pos_tol(c),
                           verify_stable=false)

Find all space-group operations `(R, τ)` of crystal `c`. Returns a
`Vector{SpacegroupOp}` in the user's original basis. The identity op is
guaranteed to be at index 1; the remaining order is unspecified.

`verify_stable=true` opts into an additional consistency check: the
computation is re-run at `pos_tol / 1000` and a warning is issued if the
operation count changes between the two tolerances (i.e. the crystal is
near a position-symmetry boundary and the returned group depends on how
permissive `pos_tol` is set).

Algorithm: Minkowski-reduce the lattice, find the lattice point group
(`pointgroup_robust`) in the reduced basis, enumerate candidate τ per R
via probe-atom differences, verify with `is_spacegroup_op`, then transform
surviving ops back to the user's basis via the integer change-of-basis
matrix.

# Examples
```jldoctest
julia> using LinearAlgebra

julia> A = Matrix{Float64}(I, 3, 3);

julia> c = Crystal(A, reshape([0.0, 0.0, 0.0], 3, 1), [:X]; coords=:fractional);

julia> length(spacegroup(c))
48

julia> spacegroup(c)[1] == one(SpacegroupOp)
true
```
"""
function spacegroup(c::Crystal; lattice_tol::Real=0.01,
                                pos_tol::Real=default_pos_tol(c),
                                verify_stable::Bool=false)
    # 1. Minkowski-reduce the lattice
    A_red = minkReduce(c.A)
    u_red, v_red, w_red = eachcol(A_red)

    # 2. Change-of-basis integer matrices.
    #   c.A · U_ro = A_red   (U_ro: reduced-coord → original-coord)
    #   A_red · U_or = c.A   (U_or: original-coord → reduced-coord; inv(U_ro))
    U_ro = round.(Int, inv(c.A) * A_red)
    U_or = round.(Int, inv(A_red) * c.A)
    abs(det(U_ro)) == 1 ||
        error("change-of-basis not unimodular: det(U_ro) = $(det(U_ro))")
    U_ro * U_or == Matrix{Int}(I, 3, 3) ||
        error("change-of-basis round-trip failed")
    norm(c.A * U_ro - A_red) < 1e-8 * opnorm(c.A) ||
        error("reduced basis does not agree with integer transform of original")

    # 3. Transform atomic positions to the reduced basis and build a Crystal
    #    on that basis. `Crystal` constructor folds positions mod 1.
    r_red = U_or * c.r
    c_red = Crystal(A_red, r_red, c.types; coords=:fractional)

    # 4. Point group of the reduced lattice
    LG_red = pointgroup_robust(u_red, v_red, w_red; tol=lattice_tol,
                                                     auto_reduce=false)

    # 5. Choose the probe atom type — the one with the fewest atoms, so the
    #    per-R candidate-τ set is as small as possible. Ties broken by
    #    first-appearance.
    probe_indices, i0 = _probe_atoms(c_red)

    # 6. For each R_red, enumerate candidate τ_red via probe-atom differences,
    #    test each with `is_spacegroup_op`, collect surviving (R_red, τ_red).
    ops_red = Tuple{Matrix{Int}, Vector{Float64}}[]
    for R in LG_red
        for τ in _find_translations_for_rotation(R, c_red, probe_indices, i0; pos_tol)
            push!(ops_red, (Matrix{Int}(R), τ))
        end
    end

    # 7. Transform ops back to the user's basis. SpacegroupOp constructor
    #    canonicalises τ to [0, 1).
    ops_out = [SpacegroupOp(U_ro * R_red * U_or, U_ro * τ_red)
               for (R_red, τ_red) in ops_red]

    # 8. Sort identity to the front (per §6.3 decision).
    e = one(SpacegroupOp)
    id_idx = findfirst(==(e), ops_out)
    id_idx === nothing &&
        error("identity operation missing from computed space group — bug")
    if id_idx != 1
        ops_out[1], ops_out[id_idx] = ops_out[id_idx], ops_out[1]
    end

    # 9. Opt-in stability check: re-run at tighter pos_tol and warn if the
    # operation count changes. Mirrors `pointgroup_robust`'s verify_stable.
    # Catches "near-miss" crystal cases — e.g. a ferroelectric where a
    # small atom displacement below pos_tol causes silent over-promotion
    # to the parent high-symmetry structure.
    if verify_stable
        tight_pos_tol = pos_tol / 1000
        tight_ops = spacegroup(c; lattice_tol, pos_tol=tight_pos_tol,
                                   verify_stable=false)
        if length(tight_ops) != length(ops_out)
            @warn "spacegroup: operation count depends on pos_tol — crystal is near a position-symmetry boundary." pos_tol ops_at_pos_tol=length(ops_out) tight_pos_tol ops_at_tight_pos_tol=length(tight_ops)
        end
    end

    return ops_out
end

"""
    snap_to_symmetry_svd(u, v, w, ops)

Snap a noisy lattice to its exact-symmetry form via singular value
decomposition of the symmetry-averaged metric tensor. Given basis vectors
`u, v, w` and lattice operations `ops` returned by [`pointgroup`](@ref)
(in lattice / integer-matrix form, the first element of its `(LG, G)`
tuple), produces:

    (a, b, c, iops, rops)

where:
- `a, b, c::Vector{Float64}` — the snapped basis vectors. Lengths and
  inter-vector angles are the symmetry-averaged values; volume is preserved.
- `iops::Vector{Matrix{Int}}` — the integer-matrix lattice operations of
  the snapped lattice (recomputed via [`pointgroup`](@ref) on the snapped
  basis).
- `rops::Vector{Matrix{Float64}}` — Cartesian rotations of the snapped lattice.

After snapping, the integer ops should satisfy `A · iops[i] · inv(A) == rops[i]`
to machine precision. Compare to the lighter [`snap_to_symmetry_avg`](@ref),
which averages each basis vector independently and is faster but less
robust at high distortion.

For accuracy-critical work — extracting symmetry operations from
experimental refinements, post-processing DFT-relaxed structures, etc. —
`pointgroup(...; tol)` followed by `snap_to_symmetry_svd(..., LG)` (using
the integer-matrix half of the `(LG, G)` tuple) gives lattice vectors and
rotations that are as exact as possible while remaining consistent with
the input.

For trusted/clean input (purely synthetic or already-snapped), this routine
is unnecessary.
"""
function snap_to_symmetry_svd(u, v, w, ops)
    A = [u v w]                       # Take the lattice basis as a matrix
    Ap = [A * k for k ∈ ops]          # Apply the integer transforms to get new basis vectors
    lengths = mean([[norm(col) for col ∈ eachcol(b)] for b ∈ Ap])
    angles = mean([[acos(p ⋅ q / norm(p) / norm(q))
                    for p ∈ eachcol(b) for q ∈ eachcol(b) if q < p] for b ∈ Ap])

    B = diagm(lengths .^ 2)
    n = length(lengths)
    # Fill in the off-diagonal components in the B matrix
    # get the "Cartesian indices" of the lower off-diagonal elements
    off_diag = [(i, j) for i ∈ 1:n for j ∈ 1:n if j < i]
    # for each index, assign the proper cos(angle)|a||b|==a⋅b
    for (i, idx) ∈ enumerate(off_diag)
        B[idx[1], idx[2]] = cos(angles[i]) * lengths[idx[1]] * lengths[idx[2]]
        B[idx[2], idx[1]] = B[idx[1], idx[2]]   # Symmetric matrix, copy elements across diagonal
    end
    s = svd(B)                        # Averaged metric matrix
    Anew = diagm(sqrt.(s.S)) * s.V'   # Getting back to a basis matrix
    T = A * inv(Anew)                 # Finding the transformation to get from old basis to new
    # This transformation contains a rotational component and a distortion component
    t = svd(T)
    rescale = cbrt(abs(det(A) / det(Anew)))
    Afinal = t.U * t.V' * Anew * rescale   # use the ortho transform of the svd to get rid of the distortion component
    u, v, w = (Afinal[:, i] for i in 1:length(u))
    if det([u v w]) < 0
        u, v, w = v, u, w
    end
    iops = pointgroup_robust(u, v, w)
    rops = to_cartesian(iops, hcat(u, v, w))
    return u, v, w, iops, rops
end

"""
    pointgroup(u, v, w; tol=0.01, verify_stable=false, auto_reduce=true)
    pointgroup(A; tol=0.01, verify_stable=false, auto_reduce=true)

Find the point group of the 3D lattice spanned by basis vectors `u, v, w`,
or equivalently by the columns of a 3×3 matrix `A`.

Returns a `Vector{Matrix{Int}}` of the symmetry operations expressed in
lattice coordinates. Each entry is a 3×3 integer matrix; if the basis
matrix is `A` then the Cartesian rotation corresponding to op `M` is
`A · M · inv(A)`. Use [`to_cartesian`](@ref) to convert when needed.

By default the input is Minkowski-reduced internally, so any non-reduced
basis is accepted. The returned operations are still expressed in the
*input* basis (a unimodular change-of-basis maps them back). Pass
`auto_reduce=false` to assert the input is already reduced — useful as a
self-check; the function will error if the assertion fails.

# Keyword arguments
- `tol::Real=0.01` — relative tolerance applied to the (volume-normalized)
  lattice. Tighter values reject more spurious candidates; looser values
  tolerate more input noise but risk over-promotion to higher symmetry.
- `verify_stable::Bool=false` — opt-in stability check. When `true`, the
  algorithm re-runs at `tol/1000` and emits a `@warn` if the operation
  count differs between the two runs (i.e. the lattice is near a
  symmetry boundary). The returned group is unchanged.

# Examples
```jldoctest
julia> using LinearAlgebra

julia> u = [1.0, 0, 0]; v = [0, 1.0, 0]; w = [0, 0, 1.0];

julia> length(pointgroup(u, v, w))
48

julia> length(pointgroup(Matrix{Float64}(I, 3, 3)))
48
```

# Extended help

This is the public entry point. It delegates to the internal
[`Spacey.pointgroup_robust`](@ref) which is tolerance-tunable and designed
for real-world noisy input. For other variants reachable via the qualified
name, see [`Spacey.pointgroup_fast`](@ref) (clean input, production speed)
and [`Spacey.pointgroup_simple`](@ref) (validation only, brute force).

The pre-v0.8 API returned a `(LG, G)` tuple of lattice and Cartesian
forms together. The Cartesian half was almost always discarded by callers
and the LG/G ambiguity was a documented footgun — see the v0.8 release
notes. Convert with `to_cartesian(LG, A)` if you actually need `G`.
"""
pointgroup(u::AbstractVector, v::AbstractVector, w::AbstractVector;
           tol::Real=0.01, verify_stable::Bool=false, auto_reduce::Bool=true) =
    pointgroup_robust(u, v, w; tol=tol, verify_stable=verify_stable,
                                auto_reduce=auto_reduce)

pointgroup(A::AbstractMatrix; tol::Real=0.01, verify_stable::Bool=false,
                              auto_reduce::Bool=true) =
    pointgroup_robust(eachcol(A)...; tol=tol, verify_stable=verify_stable,
                                      auto_reduce=auto_reduce)

"""
    Spacey.aspect_ratio(a1, a2, a3)

Return the lattice aspect ratio: longest / shortest basis vector after
Minkowski reduction. A useful diagnostic — high aspect ratios degrade the
numerical reliability of [`pointgroup`](@ref), and the underlying
`Spacey.pointgroup_robust` emits a `@warn` when the ratio exceeds 100.

Internal helper — not exported. Reach as `Spacey.aspect_ratio(...)`.

# Examples
```jldoctest
julia> Spacey.aspect_ratio([1.0, 0, 0], [0, 1.0, 0], [0, 0, 2.0])
2.0
```
"""
function aspect_ratio(a1, a2, a3)
    a = minkReduce(a1, a2, a3)[1:3]
    return maximum(norm, a) / minimum(norm, a)
end

"""
    Spacey.aspect_ratio(A)

Matrix-form wrapper around `Spacey.aspect_ratio(a1, a2, a3)`: treats the
columns of `A` as the three basis vectors. Internal helper — not exported.
"""
aspect_ratio(A) = aspect_ratio(eachcol(A)...)

end
