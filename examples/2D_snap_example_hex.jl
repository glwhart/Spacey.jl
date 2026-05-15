using Plots: length
using LinearAlgebra, Plots, StatsBase
function cell_plot(lattice)
    grid = [[0;1] [1;1] [1;0] [1;-1] [0;-1] [-1;-1] [-1;0] [-1;1]]
    d = lattice*grid
    return plot(d[1,:],d[2,:],st=:scatter,legend=:none,aspect_ratio=1)
end


function rotate_and_shift_back(v,rot,lattice,latinv)
    """ Rotate a point using the candidate rotation and then translate it back as near as possible to the starting lattice point. 
    """
    Δv = rot*v-v
    lattShift = ceil.(latinv*Δv) # number of lattice shifts to return the point to the unit cell with an origin at the original point
    neighbPts = lattice*(lattShift.-[0 1 0 1; 0 0 1 1])
    dists = [norm(vec-Δv) for vec ∈ eachcol(neighbPts)]
    idxmin = findmin(dists)[2]
    return rot*v-neighbPts[:,idxmin]
end

#function snapToLattice(u,v)
"""
   Test ideas for the snap-to-symmetry algorithm
    (Explicitly 2D only right now)
"""
u, v = ([.0624;0.001],[-0.001,15.999])

scale = √det([u v])
A = [u v]/scale # Matrix of basis vectors, normalized
iA = inv(A) # Inverse of basis vectors, convert from lattice to Cartesian
l = norm.(eachcol(A)) # lengths of input lattice vectors
# Candidate lattice vectors, in lattice coordinates (vectors that might be rotated versions of the 
# starting basis vectors) 
cL = filter(x->any(isapprox.(norm(A*x),l; atol=0.2)),[[i;j] for i ∈ -1:1 for j ∈ -1:1])

# Convert the candidate lattice vectors to Cartesian coordinates
#cLV = [A*x for x ∈ cL]

# From the candidate vectors, construct all possible _bases_ and filter out invalid bases
# Candidate bases that preserve the original volume
### This det could be done in exact integer arithmetic? We're just looking for linear independence, don't want negative copies of the same vector
RL = filter(x->isapprox(abs(det(x)),1.0; atol=0.2),[[i j] for i ∈ cL for j ∈ cL if i ≠ j])
Ap = [A*k for k ∈ RL]
#

# Find the transform for turning A into Ap
R = [p*iA for p ∈ Ap]

# Is that transform orthogonal?
Ut =  [r*r' for r ∈ R]
diffs = [x.-I(2) for x ∈ Ut]
idx = findall(norm.(diffs).<1.0)
nops = length(idx)
# Keep the transforms that are orthogonal
R = R[idx]
RL = RL[idx]
Ap = [A*k for k ∈ RL]

# Apply all the successful transforms to the original basis
# If the transforms were *exactly* orthognal, the lengths of
# and angles between the transformed vectors would be unchanged. 
# This will only be approximately true. Get the average length 
# of each vector and the angles between the basis vectors
lengths = mean([[norm(i) for i ∈ eachcol(b)] for b ∈ Ap])
angles= mean([[acos(i'*j/norm(i)/norm(j)) for i ∈ eachcol(b) for j ∈ eachcol(b) if j<i] for b ∈ Ap])

B = diagm(lengths.^2)
n = length(lengths)
# Fill in the off-diagonal components in the B matrix
# get the "Cartesian indices" of the lower off-diagonal elements
offDiag = [(i,j) for i ∈ 1:n for j ∈ 1:n if j < i]
# for each index, assign the proper cos(angle)|a||b|==a⋅b 
for (i,idx) ∈ enumerate(offDiag)
    B[idx[1],idx[2]] = cos(angles[i])*lengths[idx[1]]*lengths[idx[2]]
    B[idx[2],idx[1]] = B[idx[1],idx[2]]
end
s = svd(B)
#Anew = s.U*diagm(sqrt.(s.S))*s.V'
Anew = diagm(sqrt.(s.S))*s.V'
T = A*inv(Anew)
t = svd(T)
Afinal = t.U*t.V'*Anew
#return Afinal*scale, nops, diffs
#end

# ##
# #u, v = ([1.0+0.05;0],[1.0/2-.01,√3/2.0-.03])
#u, v = ([1.03;-.03],[.51;.87])
# #u, v = ([.706;.708],[-.705,.709])
#Afinal,n,diffs = snapToLattice(u,v)

# # UBU = [u'*B*u for u ∈ RL]
# # BMUBU = [(B-rbr)[:] for ubu ∈ UBU]
# # t=transpose(hcat(BMRBR...))[:,[1;2;4]]
# # nullspace(t)
# # # Fill in upper off-diagonals using lower off-diagonal elements
# # B = Symmetric(B,:L)
# # s = svd(B)
# # A = s.U*diagm(sqrt.(s.S))
unew=Afinal[:,1]
vnew=Afinal[:,2]
println("Volume (before, after): ",det([u v])," ",det(Afinal))
println("norms:\n",norm(unew),"\n",norm(vnew),"\n",unew'*vnew/norm(unew)/norm(vnew))
println("Num ops: ",length(idx))