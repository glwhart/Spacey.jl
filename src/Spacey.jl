module Spacey
using MinkowskiReduction
using LinearAlgebra
export pointGroup_basic, threeDrotation, pointGroup

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
function pointGroup_basic(a1,a2,a3,debug=false)
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
[-1.0 0.0 0.0; -1.0 1.0 0.0; 0.0 0.0 -0.9999999999999999]
...
```
"""
function pointGroup(a1,a2,a3)
u,v,w = minkReduce(a1,a2,a3) # Always do this first, algorithm assumes reduced basis
A = [u v w] # Define a matrix with input vectors as columns
norms=norm.([u,v,w]) # Compute the norms of the three input vectors
vol = abs(u×v⋅w) # Volume of the parallelipiped formed by the basis vectors
Ai = inv(A)
AiAiT = Ai*transpose(Ai) # Use this for checking for orthogonality
# A list of all possible lattice vectors in a rotated basis 
# These are lattice points from the vertices of the 8 cells with a corner at the origin)
# There are 27 of these (==3^3)
c = [A*[i;j;k] for i ∈ (-1,0,1) for j ∈ (-1,0,1) for k ∈ (-1,0,1)]
# Now keep only those vectors that have a norm matching one of the input vectors
# efficiency: Maybe gather three groups, according to length. That would limit the sets even more
c = c[findall([any(norm(i).≈norms) for i ∈ c])] 
# Construct all possible bases, (i.e., all combinations of c vectors), skip duplicate vectors
R = [[i j k] for i ∈ c for j ∈ c if i ≠ j for k ∈ c if j ≠ k && i ≠ k]
R = R[findall([abs(det(r))≈vol for r in R])] # Delete candidate bases with the wrong volume
RT = [transpose(i) for i ∈ R]
# This is the Uᵀ ̇U, where U transforms original basis to candidate basis
# If Tᵢ==identity then the U was a symmetry of the lattice
T = [R[i]*AiAiT*RT[i] for i ∈ 1:length(R)]
# Indices of candidate T's that match the identity
idx = findall([t≈I(3) for t ∈ T].==true)
ops = [round.(Int,Ai*R[i]) for i in idx]
return ops
end

end

