module Spacey
using MinkowskiReduction
using LinearAlgebra
using StatsBase
export pointGroup, snapToSymmetry_SVD, isagroup,
       Crystal, isSpacegroupOp, fractional, cartesian, default_pos_tol,
       crystal_system, SpacegroupOp, toCartesian, spacegroup
# Internal / not-exported (reach via `Spacey.<name>(...)`):
# - `pointGroup_robust`, `pointGroup_fast`, `pointGroup_simple` — the
#   public point-group entry point is `pointGroup` (which delegates to
#   `_robust`); `_fast` and `_simple` are validation/clean-input variants.
# - `aspectRatio`, `threeDrotation` — diagnostic / test scaffolding.
# - `snapToSymmetry_avg` — less-robust alternative to `snapToSymmetry_SVD`.

"""
    avgVecOverOps(vec, ops)

Apply each operator in `ops` to `vec` and return the average over those images
that lie within 10% of the input's norm. Used internally by
[`snapToSymmetry_avg`](@ref).

This is an internal helper; not exported.
"""
function avgVecOverOps(vec,ops)
     cands = [iop*vec for iop ∈ ops]
     cands = filter(i->norm(i-vec)<.1*norm(vec),cands)
     avgVec = sum(cands)/length(cands)
     return avgVec
end

"""
    Spacey.snapToSymmetry_avg(v1, v2, v3, ops)

Snap three basis vectors `v1, v2, v3` to a higher-symmetry triple by averaging
each vector over the images produced by `ops` (a vector of 3×3 lattice
operations). For each input vector, only images within 10% of its norm
contribute to the average — this filters the operations whose action
should fix that vector.

Internal helper — not exported. Prefer [`snapToSymmetry_SVD`](@ref), which
uses singular value decomposition of the metric tensor and is generally
more robust at high distortion. `snapToSymmetry_avg` is kept as a faster
but less-robust alternative; reach it as `Spacey.snapToSymmetry_avg(...)`.
Returns a tuple `(w1, w2, w3)`.
"""
function snapToSymmetry_avg(v1,v2,v3,ops)
     w1 = avgVecOverOps(v1,ops)
     w2 = avgVecOverOps(v2,ops)
     w3 = avgVecOverOps(v3,ops)
     return w1,w2,w3
end

"""
    Spacey.snapToSymmetry_avg(M, ops)

Matrix-form wrapper around `Spacey.snapToSymmetry_avg(v1, v2, v3, ops)`:
treats the columns of `M` as the three basis vectors and returns the
snapped vectors as a 3×3 matrix. Internal helper — not exported.
"""
function snapToSymmetry_avg(M,ops)
     res = snapToSymmetry_avg(M[:,1],M[:,2],M[:,3],ops)
     return [res[1] res[2] res[3]]
end

"""
    isagroup(members)

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

julia> isagroup([Matrix{Int}(I, 2, 2), -Matrix{Int}(I, 2, 2)])
true

julia> isagroup([[0 1; 1 0]])         # not closed: M·M = I, which isn't in the set
false
```
"""
function isagroup(members::AbstractVector{<:AbstractMatrix{<:Integer}})
    # 1) distinctness
    if length(unique(members)) < length(members)
        return false
    end

    # 2) closure
    for A in members, B in members
        C = A * B
        if !(C in members)                  # relies on exact == underneath
            return false
        end
    end

    return true
end


"""
    isagroup(members::AbstractVector{<:AbstractMatrix{<:AbstractFloat}};
             atol=1e-8, rtol=1e-8)

Floating-point variant of [`isagroup`](@ref): uses `isapprox(...; atol, rtol)`
for distinctness and closure checks. See the integer-matrix method for the
overall contract.
"""
function isagroup(members::AbstractVector{<:AbstractMatrix{<:AbstractFloat}}; atol = 1e-8, rtol = 1e-8)
    cmp(A, B) = isapprox(A, B; atol=atol, rtol=rtol)
    in_list(M, lst) = any(cmp(M, N) for N in lst)

    # 1) distinctness (approximate)
    for (k, A) in enumerate(members), B in @view members[(k+1):end]
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
        size(A) == (3, 3) || error("A must be 3×3")
        size(r, 1) == 3 || error("r must have 3 rows (one per spatial dimension)")
        size(r, 2) == length(types) ||
            error("r has $(size(r,2)) columns but types has length $(length(types))")
        size(r, 2) > 0 || error("empty crystal (no atoms) is not supported")
        coords ∈ (:fractional, :cartesian) ||
            error("coords must be :fractional or :cartesian (got $(repr(coords)))")
        A64 = Float64.(A)
        # Reject (near-)singular lattices. Scale-invariant test: compare |det|
        # to eps · ‖A‖³, i.e. the precision at which det is distinguishable
        # from zero for a matrix of this scale.
        abs(det(A64)) > eps(Float64) * opnorm(A64)^3 ||
            error("A is singular or near-singular (det = $(det(A64)))")
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

Default position-matching tolerance used by [`isSpacegroupOp`](@ref) and
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
[`pointGroup`](@ref).

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
    LG = pointGroup_robust(u, v, w; tol=lattice_tol, auto_reduce=false)
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
    isSpacegroupOp(R, τ, c::Crystal; tol=default_pos_tol(c))

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

julia> isSpacegroupOp(I3, [0.0, 0.0, 0.0], c; tol=1e-8)
true

julia> isSpacegroupOp(I3, [0.5, 0.0, 0.0], c; tol=1e-8)
false
```
"""
function isSpacegroupOp(R::AbstractMatrix{<:Real}, τ::AbstractVector{<:Real},
                        c::Crystal; tol::Real=default_pos_tol(c))
    size(R) == (3, 3) || error("R must be 3×3")
    length(τ) == 3 || error("τ must have length 3")
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
`toCartesian(op, A)`.

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

function _canonicalize_τ(τ::AbstractVector, tol::Real=1e-6)
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
    maximum(abs, Rinv_f .- Rinv) < 1e-8 ||
        error("inv(SpacegroupOp): R⁻¹ is not integer (det(R) ≠ ±1?)")
    SpacegroupOp(Rinv, -Rinv * op.τ)
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
    toCartesian(op::SpacegroupOp, A::AbstractMatrix)

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
"""
toCartesian(op::SpacegroupOp, A::AbstractMatrix) =
    (A * op.R * inv(A), A * op.τ)

"""
    toCartesian(op::AbstractMatrix{<:Integer}, A::AbstractMatrix)
    toCartesian(LG::AbstractVector{<:AbstractMatrix{<:Integer}}, A::AbstractMatrix)

Convert a lattice-coordinate point-group operation (a single integer
matrix) — or a whole vector of them, like the result of [`pointGroup`](@ref) —
to Cartesian rotation form. Returns `A · op · inv(A)` for the single-op
method and `[A · op · inv(A) for op in LG]` for the vector method.

These overloads exist because `pointGroup` returns just the integer-matrix
form (since v0.8); use these helpers if you need the Cartesian rotations.

# Examples
```jldoctest
julia> using Spacey, LinearAlgebra

julia> A = Matrix{Float64}(I, 3, 3);

julia> LG = pointGroup(A);

julia> G = toCartesian(LG, A);   # Cartesian rotations of the cubic point group

julia> length(G)
48

julia> G[findfirst(==(Matrix{Int}(I, 3, 3)), LG)]   # identity in Cartesian
3×3 Matrix{Float64}:
 1.0  0.0  0.0
 0.0  1.0  0.0
 0.0  0.0  1.0
```
"""
toCartesian(op::AbstractMatrix{<:Integer}, A::AbstractMatrix) = A * op * inv(A)

toCartesian(LG::AbstractVector{<:AbstractMatrix{<:Integer}}, A::AbstractMatrix) =
    [A * op * inv(A) for op in LG]

"""
    Spacey.threeDrotation(u, v, w, α, β, γ)

Rotate the basis vectors `u, v, w` by Euler angles `α, β, γ` (the
yaw–pitch–roll convention used in test scaffolding). Returns the rotated
triple `(u', v', w')` as a tuple of three vectors.

Internal helper — not exported. Used by tests to verify that
symmetry-finding routines are invariant under arbitrary lattice
orientation. Reach as `Spacey.threeDrotation(...)`.

The rotation matrix is built from successive rotations about the z, y, and
z axes (matching the order in the formula). For zero angles the identity
is returned.

# Examples
```jldoctest
julia> u, v, w = Spacey.threeDrotation([1.0,0,0], [0,1.0,0], [0,0,1.0], 0.0, 0.0, 0.0);

julia> u
3-element Vector{Float64}:
 1.0
 0.0
 0.0
```
"""
function threeDrotation(u,v,w,α,β,γ)
A = [u v w]
R = [[cos(α)cos(β) cos(α)sin(β)sin(γ)-sin(α)cos(γ) cos(α)sin(β)cos(γ)+sin(α)sin(γ)];
     [sin(α)cos(β) sin(α)sin(β)sin(γ)+cos(α)cos(γ) sin(α)sin(β)cos(γ)-cos(α)sin(γ)];
     [-sin(β)      cos(β)sin(γ)                    cos(β)cos(γ)                   ]]
A = R*A
return A[:,1],A[:,2],A[:,3] 
end

"""
    Spacey.pointGroup_simple(a1, a2, a3, debug=false)

Brute-force enumeration of the point group of a 3D lattice. Iterates over
every 3×3 candidate matrix with entries in `{-1, 0, 1}` (3⁹ = 19683
matrices), retains those whose action on the basis preserves the metric
tensor, and returns the survivors as Cartesian rotations.

Internal — not exported. The simplest correct implementation; used to
validate the more efficient `Spacey.pointGroup_fast` and the public
[`pointGroup`](@ref) (which delegates to `Spacey.pointGroup_robust`).
It performs strict (`isapprox` with default tolerance) equality checks,
so it is most reliable on noiseless / synthetic input.

If `debug=true`, returns the candidate `T = UᵀU` matrices instead of the
filtered ops, for use when diagnosing failures.

# Examples
```jldoctest
julia> using Spacey

julia> u = [1.0, 0, 0]; v = [0.5, √3/2, 0]; w = [0.0, 0, √(8/3)];

julia> length(Spacey.pointGroup_simple(u, v, w))
24
```
"""
function pointGroup_simple(a1,a2,a3,debug=false)
u,v,w = minkReduce(a1,a2,a3)
A = [u v w] # Put the lattice vectors as columns in matrix A
B = inv(A)*transpose(inv(A)) # Use this for checking for orthogonality
# A list of all possible lattice vectors in a rotated basis
c = [A*[i;j;k] for i ∈ (-1,0,1) for j ∈ (-1,0,1) for k ∈ (-1,0,1)]
# A list of all possible bases, (i.e., all combinations of c vectors)
R = [[i j k] for i ∈ c for j ∈ c for k ∈ c]
RT = [transpose(i) for i ∈ R]
# This is the U^T*U, where U transforms original basis to candidate basis
T = [R[i]*B*RT[i] for i ∈ 1:length(R)]
if debug return T end
# If T==identity then the U was a symmetry of the lattice
idx = findall([t≈I(3) for t ∈ T].==true)
Ai = inv(A)
ops = [Ai*R[i] for i in idx]
return ops
end


"""
    Spacey.pointGroup_fast(a1, a2, a3)

Production-speed point-group finder for an exact / noiseless 3D lattice.
Faster than `Spacey.pointGroup_simple` by filtering candidate basis
combinations by length and volume before checking orthogonality, but uses
strict `isapprox` tolerance and so is best suited to clean inputs.

Internal — not exported. For real-world (noisy) input use the public
[`pointGroup`](@ref), which exposes a tolerance keyword.

Returns operations as integer matrices in lattice coordinates.

# Examples
```jldoctest
julia> using Spacey

julia> u = [1.0, 0, 0]; v = [0.5, √3/2, 0]; w = [0.0, 0, √(8/3)];

julia> length(Spacey.pointGroup_fast(u, v, w))
24
```
"""
function pointGroup_fast(a1,a2,a3)
u,v,w = minkReduce(a1,a2,a3) # Always do this first, algorithm assumes reduced basis
A = [u v w] # Define a matrix with input vectors as columns
Ai = inv(A) 
AiAiT = Ai*transpose(Ai) # Use this for checking for orthogonality
norms=norm.([u,v,w]) # Compute the norms of the three input vectors
vol = abs(u×v⋅w) # Volume of the parallelipiped formed by the basis vectors


# A list of all possible lattice vectors in a rotated basis 
# These are lattice points from the vertices of the 8 cells with a corner at the origin)
# There are 27 of these (==3^3)
c = [A*[i,j,k] for i ∈ (-1,0,1) for j ∈ (-1,0,1) for k ∈ (-1,0,1)]
# Now keep only those vectors that have a norm matching one of the input vectors
# efficiency: Gather three groups, according to length. This limits the candidates even more
c1 = c[findall([norm(i)≈norms[1] for i ∈ c])] # All vectors with first norm
c2 = c[findall([norm(i)≈norms[2] for i ∈ c])] # All vectors with second norm
c3 = c[findall([norm(i)≈norms[3] for i ∈ c])] # All vectors with third norm
# Construct all possible bases, (i.e., all combinations of c vectors), skip duplicate vectors
R = [[i j k] for i ∈ c1 for j ∈ c2 if !(i≈j) for k ∈ c3 if !(i≈k) && !(j≈k)]
R = R[findall([abs(det(r))≈vol for r in R])] # Delete candidate bases with the wrong volume
# The cross product is slightly (<1%) faster
#R = R[findall([abs(r[1]×r[2]⋅r[3])≈vol for r in R])] # Delete candidate bases with the wrong volume
RT = [transpose(i) for i ∈ R]
# This is the Uᵀ ̇U, where U transforms original basis to candidate basis
# If Tᵢ==identity then the U was a symmetry of the lattice
T = [R[i]*AiAiT*RT[i] for i ∈ 1:length(R)]
# Indices of candidate T's that match the identity
idx = findall([t≈I(3) for t ∈ T])
# Convert the transformations to integer matrices (formally they should be)
ops = [round.(Int,Ai*R[i]) for i in idx]
return ops
end

"""
    Spacey.pointGroup_robust(u, v, w; tol=0.01, verify_stable=false, auto_reduce=true)

Tolerance-tunable point-group finder for noisy real-world input. Returns
a `Vector{Matrix{Int}}` of the lattice-coordinate symmetry operations —
the same form the public [`pointGroup`](@ref) returns. For Cartesian
rotations, pass the result through [`toCartesian`](@ref).

Internal — not exported. The public entry point [`pointGroup`](@ref) is a
thin wrapper around this function (with the same defaults). Reach this
form directly only when explicitly disambiguating between point-group
variants (e.g. comparing against `Spacey.pointGroup_fast`).

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

julia> length(Spacey.pointGroup_robust(u, v, w))
48
```
"""
function pointGroup_robust(u, v, w; tol=0.01, verify_stable::Bool=false,
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
    if !(orthogonalityDefect(u,v,w)≈orthogonalityDefect(minkReduce(u,v,w)[1:3]...))
        error("Input basis for 'pointGroup' is not Minkowski-reduced. Either pass `auto_reduce=true` (the default) or run `minkReduce` first.")
    end
end
inputVol = ∛(abs(u×v⋅w)) # Rescale the basis to have a volume of 1, avoid floating point issues
u, v, w = u ./ inputVol, v ./ inputVol, w ./ inputVol

norms=norm.([u,v,w]) # Compute the norms of the three input vectors
ar = maximum(norms) / minimum(norms)
if ar > 100
    @warn "Aspect ratio is $(round(ar,digits=1)). Results may be unreliable for ratios above ~500."
end

A = [u v w] # Define a matrix with input vectors as columns
Ai = inv(A)
vol = abs(u×v⋅w) # Volume of the parallelipiped formed by the basis vectors

# A list of all possible lattice vectors in a rotated basis. These are lattice points from the vertices of the 8 cells that have a corner at the origin. There are 27 of these (==3^3)
c = [A*[i,j,k] for i ∈ (-1,0,1) for j ∈ (-1,0,1) for k ∈ (-1,0,1)]
# Now keep only those vectors that have a norm close the norm one of the input vectors
# efficiency: Gather three groups, according to length. This limits the candidates even more
c1 = c[findall([isapprox(norms[1],norm(i),rtol=tol) for i ∈ c])] # All vectors with first norm
c2 = c[findall([isapprox(norms[2],norm(i),rtol=tol) for i ∈ c])] # All vectors with second norm
c3 = c[findall([isapprox(norms[3],norm(i),rtol=tol) for i ∈ c])] # All vectors with third norm
# Construct all candidate bases, Rc (i.e., all combinations of c vectors), skip duplicate vectors.
A′ = [[i j k] for i ∈ c1 for j ∈ c2 if !(i≈j) for k ∈ c3 if !(i≈k) && !(j≈k)] # All candidate bases
A′ = A′[findall([isapprox(abs(det(i)),vol,rtol=tol*min(norms...)) for i in A′])] # Delete candidate bases with the wrong volume
Rc = [i*Ai for i ∈ A′] # Compute the candidate rotations from the candidate bases

# This is the Uᵀ ̇U, where U transforms original basis to candidate basis
# If Tᵢ==identity then the Rc is orthogonal and is a symmetry of the lattice
T = [transpose(rc)*rc for rc ∈ Rc]
# Indices of candidate T's that match the identity
idx = findall([isapprox(t,I(3),rtol=tol) for t ∈ T])
Rc = Rc[idx]
T = T[idx]
# Convert the transformations to lattice coordinates representation (round to integer matrices; formally they should be)
ops = [round.(Int,Ai*i*A) for i in Rc] # Need the 'Int' so integers are returned
# Get norms of deviation from orthogonal case
tn = [norm(t-I(3)) for t ∈ T]
tp = sortperm(tn) # Sort by deviation
# Find the largest number of (sorted) ops that form a group.
maxl = 48
for il ∈ [48,24,16,12,8,4,2] # These are the only possible group sizes for a 3D lattice
     if il > length(idx) continue end
     if isagroup(ops[tp[1:il]]) # Keep the largest set that is a group
          maxl = il
          break
     end
end
result_ops = ops[tp][1:maxl]
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
    tight_ops = pointGroup_robust(u, v, w; tol=tight_tol, verify_stable=false,
                                            auto_reduce=false)
    if length(tight_ops) != length(result_ops)
        @warn "pointGroup_robust: group size depends on tolerance — lattice is near a symmetry boundary." tol group_at_tol=length(result_ops) tight_tol group_at_tight_tol=length(tight_ops)
    end
end
return result_ops
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
(`pointGroup_robust`) in the reduced basis, enumerate candidate τ per R
via probe-atom differences, verify with `isSpacegroupOp`, then transform
surviving ops back to the user's basis via the integer change-of-basis
matrix.

See `spacegroup_plan.md` for design notes and `phase2_plan.md` for
derivations.

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
    LG_red = pointGroup_robust(u_red, v_red, w_red; tol=lattice_tol,
                                                     auto_reduce=false)

    # 5. Choose the probe atom type — the one with the fewest atoms, so the
    #    per-R candidate-τ set is as small as possible. Ties broken by
    #    first-appearance.
    types_unique = unique(c.types)
    counts = [count(==(t), c.types) for t in types_unique]
    probe_type = types_unique[argmin(counts)]
    probe_indices = findall(==(probe_type), c.types)
    i0 = probe_indices[1]

    # 6. For each R_red, enumerate candidate τ_red via probe-atom differences,
    #    test each with `isSpacegroupOp`, collect surviving (R_red, τ_red).
    ops_red = Tuple{Matrix{Int}, Vector{Float64}}[]
    for R in LG_red
        image_i0 = R * c_red.r[:, i0]
        for j in probe_indices
            τ = mod.(c_red.r[:, j] .- image_i0, 1.0)
            if isSpacegroupOp(R, τ, c_red; tol=pos_tol)
                push!(ops_red, (Matrix{Int}(R), τ))
            end
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
    # operation count changes. Mirrors `pointGroup_robust`'s verify_stable.
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
    snapToSymmetry_SVD(u, v, w, ops)

Snap a noisy lattice to its exact-symmetry form via singular value
decomposition of the symmetry-averaged metric tensor. Given basis vectors
`u, v, w` and lattice operations `ops` returned by [`pointGroup`](@ref)
(in lattice / integer-matrix form, the first element of its `(LG, G)`
tuple), produces:

    (a, b, c, iops, rops)

where:
- `a, b, c::Vector{Float64}` — the snapped basis vectors. Lengths and
  inter-vector angles are the symmetry-averaged values; volume is preserved.
- `iops::Vector{Matrix{Int}}` — the integer-matrix lattice operations of
  the snapped lattice (recomputed via [`pointGroup`](@ref) on the snapped
  basis).
- `rops::Vector{Matrix{Float64}}` — Cartesian rotations of the snapped lattice.

After snapping, the integer ops should satisfy `A · iops[i] · inv(A) == rops[i]`
to machine precision. Compare to the lighter [`snapToSymmetry_avg`](@ref),
which averages each basis vector independently and is faster but less
robust at high distortion.

For accuracy-critical work — extracting symmetry operations from
experimental refinements, post-processing DFT-relaxed structures, etc. —
`pointGroup(...; tol)` followed by `snapToSymmetry_SVD(..., LG)` (using
the integer-matrix half of the `(LG, G)` tuple) gives lattice vectors and
rotations that are as exact as possible while remaining consistent with
the input.

For trusted/clean input (purely synthetic or already-snapped), this routine
is unnecessary.
"""
function snapToSymmetry_SVD(u,v,w,ops)
A = [u v w] # Take the lattice basis as a matrix 
Ap = [A*k for k ∈ ops] # Apply the integer tranforms to get new basis vectors
lengths = mean([[norm(i) for i ∈ eachcol(b)] for b ∈ Ap]) 
angles= mean([[acos(i⋅j/norm(i)/norm(j)) for i ∈ eachcol(b) for j ∈ eachcol(b) if j<i] for b ∈ Ap])

B = diagm(lengths.^2)
n = length(lengths)
# Fill in the off-diagonal components in the B matrix
# get the "Cartesian indices" of the lower off-diagonal elements
offDiag = [(i,j) for i ∈ 1:n for j ∈ 1:n if j < i]
# for each index, assign the proper cos(angle)|a||b|==a⋅b 
for (i,idx) ∈ enumerate(offDiag)
     B[idx[1],idx[2]] = cos(angles[i])*lengths[idx[1]]*lengths[idx[2]]
     B[idx[2],idx[1]] = B[idx[1],idx[2]] # Symmetric matrix, copy elements across diagonal
end
s = svd(B) # Averaged metric matrix
Anew = diagm(sqrt.(s.S))*s.V' # Getting back to a basis matrix
T = A*inv(Anew) # Finding the transformation to get from old basis to new
# This transformation contains a rotational component and a distortion component
t = svd(T)
rescale = cbrt(abs(det(A)/det(Anew)))
Afinal = t.U*t.V'*Anew*rescale # use the the ortho transform of the svd to get rid of the distortion component
u,v,w=[Afinal[:,i] for i ∈ 1:length(u)]
if det([u v w]) < 0 
     u,v,w = v,u,w
end
iops = pointGroup_robust(u,v,w)
rops = toCartesian(iops, hcat(u,v,w))
return u,v,w,iops,rops
end

"""
    pointGroup(u, v, w; tol=0.01, verify_stable=false, auto_reduce=true)
    pointGroup(A; tol=0.01, verify_stable=false, auto_reduce=true)

Find the point group of the 3D lattice spanned by basis vectors `u, v, w`,
or equivalently by the columns of a 3×3 matrix `A`.

Returns a `Vector{Matrix{Int}}` of the symmetry operations expressed in
lattice coordinates. Each entry is a 3×3 integer matrix; if the basis
matrix is `A` then the Cartesian rotation corresponding to op `M` is
`A · M · inv(A)`. Use [`toCartesian`](@ref) to convert when needed.

By default the input is Minkowski-reduced internally, so any non-reduced
basis is accepted. The returned operations are still expressed in the
*input* basis (a unimodular change-of-basis maps them back). Pass
`auto_reduce=false` to assert the input is already reduced — useful as a
self-check; the function will error if the assertion fails.

(The pre-v0.8 API returned a `(LG, G)` tuple of lattice and Cartesian
forms together. The Cartesian half was almost always discarded by callers
and the LG/G ambiguity was a documented footgun — see the v0.8 release
notes. Convert with `toCartesian(LG, A)` if you actually need `G`.)

# Keyword arguments
- `tol::Real=0.01` — relative tolerance applied to the (volume-normalized)
  lattice. Tighter values reject more spurious candidates; looser values
  tolerate more input noise but risk over-promotion to higher symmetry.
- `verify_stable::Bool=false` — opt-in stability check. When `true`, the
  algorithm re-runs at `tol/1000` and emits a `@warn` if the operation
  count differs between the two runs (i.e. the lattice is near a
  symmetry boundary). The returned group is unchanged.

This is the public entry point. It delegates to the internal
[`Spacey.pointGroup_robust`](@ref) which is tolerance-tunable and designed
for real-world noisy input. For other variants reachable via the qualified
name, see [`Spacey.pointGroup_fast`](@ref) (clean input, production speed)
and [`Spacey.pointGroup_simple`](@ref) (validation only, brute force).

# Examples
```jldoctest
julia> using LinearAlgebra

julia> u = [1.0, 0, 0]; v = [0, 1.0, 0]; w = [0, 0, 1.0];

julia> length(pointGroup(u, v, w))
48

julia> length(pointGroup(Matrix{Float64}(I, 3, 3)))
48
```
"""
pointGroup(u::AbstractVector, v::AbstractVector, w::AbstractVector;
           tol::Real=0.01, verify_stable::Bool=false, auto_reduce::Bool=true) =
    pointGroup_robust(u, v, w; tol=tol, verify_stable=verify_stable,
                                auto_reduce=auto_reduce)

pointGroup(A::AbstractMatrix; tol::Real=0.01, verify_stable::Bool=false,
                              auto_reduce::Bool=true) =
    pointGroup_robust(A[:,1], A[:,2], A[:,3]; tol=tol, verify_stable=verify_stable,
                                              auto_reduce=auto_reduce)

"""
    Spacey.aspectRatio(a1, a2, a3)

Return the lattice aspect ratio: longest / shortest basis vector after
Minkowski reduction. A useful diagnostic — high aspect ratios degrade the
numerical reliability of [`pointGroup`](@ref), and the underlying
`Spacey.pointGroup_robust` emits a `@warn` when the ratio exceeds 100.

Internal helper — not exported. Reach as `Spacey.aspectRatio(...)`.

# Examples
```jldoctest
julia> Spacey.aspectRatio([1.0, 0, 0], [0, 1.0, 0], [0, 0, 2.0])
2.0
```
"""
function aspectRatio(a1,a2,a3)
     a = minkReduce(a1,a2,a3)[1:3]
     return max(norm(a[1]),norm(a[2]),norm(a[3]))/min(norm(a[1]),norm(a[2]),norm(a[3]))
end

"""
    Spacey.aspectRatio(A)

Matrix-form wrapper around `Spacey.aspectRatio(a1, a2, a3)`: treats the
columns of `A` as the three basis vectors. Internal helper — not exported.
"""
function aspectRatio(A)
     return aspectRatio(A[:,1],A[:,2],A[:,3])
end

end