module Spacey
using MinkowskiReduction
using LinearAlgebra
using StatsBase
export pointGroup_fast, pointGroup_simple, threeDrotation,
       pointGroup, pointGroup_robust, snapToSymmetry_SVD, isagroup,
       snapToSymmetry_avg, aspectRatio,
       Crystal, isSpacegroupOp, fractional, cartesian, default_pos_tol,
       spacegroup

""" averageOverOps(vec,ops) 

Apply each operator in the list to the input vector. Take the average over the vectors that are approximately invariant
"""
function avgVecOverOps(vec,ops)
     cands = [iop*vec for iop ∈ ops]
     cands = filter(i->norm(i-vec)<.1*norm(vec),cands)
     avgVec = sum(cands)/length(cands)
     return avgVec
end

""" snapToSymmetry_avg(v1,v2,v3,ops)

Average three basis vectors over the operators

(The average is only over the resultant vectors that are ≈ to the original.)
"""
function snapToSymmetry_avg(v1,v2,v3,ops)
     w1 = avgVecOverOps(v1,ops)
     w2 = avgVecOverOps(v2,ops)
     w3 = avgVecOverOps(v3,ops)
     return w1,w2,w3
end

""" snapToSymmetry_avg(M,ops) 

Average three basis vectors (columns of M) over the operators.

(The average is only over the resultant vectors that are ≈ to the original.)
"""
function snapToSymmetry_avg(M,ops)
     res = snapToSymmetry_avg(M[:,1],M[:,2],M[:,3],ops)
     return [res[1] res[2] res[3]]
end

""" isagroup(members::Vector{<:AbstractMatrix{<:Integer}})

Determine whether `members` (matrices with *integer* elements) form a group
under matrix multiplication using *exact* equality.
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


""" isagroup(members::Vector{<:AbstractMatrix{<:AbstractFloat}}; atol = 1e-8, rtol = 1e-8)

Determine whether `members` (matrices with floating-point elements) form a
group under matrix multiplication, using `isapprox` with the provided
tolerances to handle finite-precision errors.
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

Numeric input is accepted as any `AbstractMatrix{<:Real}` / `AbstractVector{<:Real}`
and converted to `Float64` once at construction.
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
        r64 = coords === :cartesian ? inv(A64) * Float64.(r) : Float64.(r)
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
"""
fractional(c::Crystal) = c.r

"""
    cartesian(c::Crystal)

Return the 3×N matrix of atomic positions in Cartesian coordinates (same basis
and units as `c.A`).
"""
cartesian(c::Crystal) = c.A * c.r

"""
    default_pos_tol(c::Crystal)

Default position-matching tolerance for symmetry detection. Equal to
`0.01 · (V/N)^(1/3)` where `V = |det(A)|` and `N` is the number of atoms —
1% of the characteristic atom separation, unit-agnostic. See
`designDiscussions.md` for the rationale behind the 1% choice.
"""
default_pos_tol(c::Crystal) = 0.01 * (abs(det(c.A)) / size(c.r, 2))^(1/3)

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
threeDrotation(u,v,w,α,β,γ)

Rotate a basis by three angles to any orientation. Useful for building robust unit tests. 

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
Generate the symmetry operations of a lattice, defined by three 3-vectors
(This function is _very_ simple but not efficient. Used to test more efficient algorithms.)

```juliadoctest
julia> u = [1,0,0]; v = [.5,√3/2,0]; w = [0,0,√(8/3)];
julia> pointGroup_basic(u,v,w)
24-element Array{Array{Float64,2},1}:
[-1.0 0.0 0.0; -1.0 1.0 0.0; 0.0 0.0 -0.9999999999999999]
...
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
Generate the symmetry operations of a lattice, defined by three 3-vectors.
This function aims to be more efficient than `pointGroup_basic` and so is more complex

```juliadoctest
julia> u = [1,0,0]; v = [.5,√3/2,0]; w = [0,0,√(8/3)];
julia> pointGroup(u,v,w)
24-element Array{Array{Float64,2},1}:
[-1 -1 0; 0 1 0; 0 0 -1]
...
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

""" pointgroup_robust(a1, a2, a3)

Calculate the pointGroup using an epsilon based on input lattice 

The routine works in a similar fashion to the pointGroup_fast routine but
finite precision comparisons use an epsilon scaled to the input. The 
epsilon is quite large, 10% (may change) of the smallest scale of the 
input. With sufficient testing, this routine may become the defacto standard
for the Spacey package.

     ```juliadoctest
     julia>  u = [1,0,1e-4]; v = [.5,√3/2,-1e-5]; w = [0,1e-4,√(8/3)];
     julia> pointGroup_robust(u,v,w)
     24-element Array{Array{Float64,2},1}:
     [-1.0 0.0 0.0; -1.0 1.0 0.0; 0.0 0.0 -0.9999999999999999]
     ...
     ```
"""
function pointGroup_robust(u,v,w;tol=0.01, verify_stable::Bool=false)
# Mink reduction can change the basis even when the basis is already reduced (degenerate cases). So don't do it here. But do check that no reduction is needed.
if !(orthogonalityDefect(u,v,w)≈orthogonalityDefect(minkReduce(u,v,w)[1:3]...))
    error("Input basis for 'pointGroup' is not reduced. Use 'minkReduce' to pick the shortest basis vectors.")
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
result_Rc  = Rc[tp][1:maxl]
if verify_stable
    tight_tol = tol / 1000
    tight_ops, _ = pointGroup_robust(u, v, w; tol=tight_tol, verify_stable=false)
    if length(tight_ops) != length(result_ops)
        @warn "pointGroup_robust: group size depends on tolerance — lattice is near a symmetry boundary." tol group_at_tol=length(result_ops) tight_tol group_at_tight_tol=length(tight_ops)
    end
end
return result_ops, result_Rc
end

"""
    spacegroup(c::Crystal; ...)

Find the space-group operations of crystal `c`. Scheduled for Phase 2 per
`spacegroup_plan.md`; currently raises an error to prevent silent misuse.
"""
function spacegroup(c::Crystal)
    error("spacegroup: not yet implemented — scheduled for Phase 2 (see spacegroup_plan.md)")
end

""" snapToSymmetry_SVD(u,v,w,ops)

    Adjust input vectors and atomic basis to be an exact match to symmetry found. Adjust symmetries to be exact orthonormal transforms (to machine precision)

     In most applications, where robustness/accuracy is the most important consideration (rather than speed), one should probably always call the "robust" pointGroup finder and then follow up with a call to this routine. If the input is trustworthy (highly accurate), then calling this routine would be unnecessary.
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
iops,rops = pointGroup_robust(u,v,w)
return u,v,w,iops,rops
end

""" pointGroup(basisMatrix)

Calculate the point group of a lattice using a matrix of column vectors as the basis. The routine is a wrapper for the `pointGroup_robust` routine. See docstring for that routine for more details.
"""
function pointGroup(A;tol=0.1)
     return pointGroup_robust(A[:,1],A[:,2],A[:,3];tol=tol)
end

""" aspectRatio(a1,a2,a3)

Calculate the aspect ratio of a lattice.

The aspect ratio is the ratio of the longest to shortest lattice vector.
"""
function aspectRatio(a1,a2,a3)
     a = minkReduce(a1,a2,a3)[1:3]
     return max(norm(a[1]),norm(a[2]),norm(a[3]))/min(norm(a[1]),norm(a[2]),norm(a[3]))
end

""" aspectRatio(A)

Calculate the aspect ratio of a lattice using a matrix of column vectors as the basis. The routine is a wrapper for the `aspectRatio(a1,a2,a3)` routine.
"""
function aspectRatio(A)
     return aspectRatio(A[:,1],A[:,2],A[:,3])
end

end