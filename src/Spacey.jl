module Spacey
using MinkowskiReduction
using LinearAlgebra
using StatsBase
export pointGroup_fast, pointGroup_simple, threeDrotation, 
       pointGroup, pointGroup_robust, snapToSymmetry_SVD, isagroup,
       snapToSymmetry_avg

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

""" isagroup(members)
Check and see if a list of operators form a valid group. 
Assumes that the operators are integer matrices.
     
Checks that 1) each member is unique and 2) that the product of 
any two members is still in the list (closure).
"""
function isagroup(members)
     # Check that each member is unique
     if length(unique(members)) < length(members)
     return false
     end
     # Check that there is closure 
     for i ∈ members
     for j ∈ members
          if !(i*j ∈ members)
               return false
          end
     end
     end 
     return true
end

struct Crystal
     a1
     a2 
     a3
     r::Array{Float64,2}
     a::Array{Int}
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
function pointGroup_robust(u,v,w,tol=0.1)
# Mink reduction can change the basis even when the basis is already reduced (degenerate cases). So don't do it here. But do check that no reduction is needed.
if !(orthogonalityDefect(u,v,w)≈orthogonalityDefect(minkReduce(u,v,w)[1:3]...))
    error("Input basis for 'pointGroup' is not reduced. Use 'minkReduce' to pick the shortest basis vectors.")
end
A = [u v w] # Define a matrix with input vectors as columns
Ai = inv(A) 
norms=norm.([u,v,w]) # Compute the norms of the three input vectors
ε = tol*min(norms...) # Scale factor for comparisons (unit tests must decide correct rescaling)
vol = abs(u×v⋅w) # Volume of the parallelipiped formed by the basis vectors

# A list of all possible lattice vectors in a rotated basis. These are lattice points from the vertices of the 8 cells that have a corner at the origin. There are 27 of these (==3^3)
c = [A*[i,j,k] for i ∈ (-1,0,1) for j ∈ (-1,0,1) for k ∈ (-1,0,1)]
# Now keep only those vectors that have a norm matching one of the input vectors
# efficiency: Gather three groups, according to length. This limits the candidates even more
c1 = c[findall([isapprox(norms[1],norm(i),atol=ε) for i ∈ c])] # All vectors with first norm
c2 = c[findall([isapprox(norms[2],norm(i),atol=ε) for i ∈ c])] # All vectors with second norm
c3 = c[findall([isapprox(norms[3],norm(i),atol=ε) for i ∈ c])] # All vectors with third norm
# Construct all candidate bases, Rc (i.e., all combinations of c vectors), skip duplicate vectors.
Rc = [[i j k] for i ∈ c1 for j ∈ c2 if i !≈ j for k ∈ c3 if i ≈ k && j !≈ k]
Rc = Rc[findall([isapprox(abs(det(i)),vol,rtol=ε) for i in Rc])]  # Delete candidate bases with the wrong volume
Rc = [i*Ai for i ∈ Rc] # Compute the candidate rotations from the candidate bases

# This is the Uᵀ ̇U, where U transforms original basis to candidate basis
# If Tᵢ==identity then the Rc is orthogonal and is a symmetry of the lattice
T = [transpose(rc)*rc for rc ∈ Rc]
# Indices of candidate T's that match the identity
idx = findall([norm(t-I(3)) < ε for t ∈ T])
# Convert the transformations to integer matrices (formally they should be)
ops = [round.(Int,Ai*i*A) for i in Rc[idx]] # Need the 'Int' so integers are returned
return ops, Rc[idx]
end

""" spaceGroup(a1, a2, a3, r, ele) 

Calculate the spacegroup of a crystal

Relies on the `pointgroup_robust` function for point group operations
(need to replace the doctest)     
     ```juliadoctest
     julia>  u = [1,0,1e-4]; v = [.5,√3/2,-1e-5]; w = [0,1e-4,√(8/3)];
     julia> pointGroup_robust(u,v,w)
     24-element Array{Array{Float64,2},1}:
     [-1.0 0.0 0.0; -1.0 1.0 0.0; 0.0 0.0 -0.9999999999999999]
     ...
     ```
"""
function spacegroup(c::Crystal)
     return true
end

""" Adjust input vectors and atomic basis to be an exact match to symmetry
found. Adjust symmetries to be exact orthonormal transforms (to machine precision)

In most applications, where robustness/accuracy is the most important
consideration (rather than speed), one should probably always call the "robust"
pointGroup finder and then follow up with a call to this routine. If the input
is trustworthy (highly accurate), then calling this routine would be unnecessary.
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
ops,_ = pointGroup_robust(u,v,w)
if det([u v w]) < 0 
     u,v,w = v,u,w
end

return u,v,w,ops
end

""" pointGroup(basisMatrix)

Calculate the point group of a lattice using a matrix of column vectors as the basis. The routine is a wrapper for the `pointGroup_robust` routine. See docstring for that routine for more details.
"""
function pointGroup(A)
     return pointGroup_robust(A[:,1],A[:,2],A[:,3])
end
end 

