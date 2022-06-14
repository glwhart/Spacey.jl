using MinkowskiReduction
a1 = [1/16 - .0001, .001, .0001]
a2 = [-0.001, 16-.0001, -.0001]
a3 = [-.0001, .0001, 1.51]
u,v,w = minkReduce(a1,a2,a3) # Always do this first, algorithm assumes reduced basis
A = [u v w] # Define a matrix with input vectors as columns
Ai = inv(A) 
Aiᵀ = transpose(Ai)
norms=norm.([u,v,w]) # Compute the norms of the three input vectors
ε = 0.1min(norms...) # Scale factor for comparisons (unit tests must decide correct rescaling)
vol = abs(u×v⋅w) # Volume of the parallelipiped formed by the basis vectors
# A list of all possible lattice vectors in a rotated basis 
# These are lattice points from the vertices of the 8 cells with a corner at the origin)
# There are 27 of these (==3^3)
c = [A*[i,j,k] for i ∈ (-1,0,1) for j ∈ (-1,0,1) for k ∈ (-1,0,1)]
# Now keep only those vectors that have a norm matching one of the input vectors
# efficiency: Gather three groups, according to length. This limits the candidates even more
c1 = c[findall([isapprox(norms[1],norm(i),rtol=ε) for i ∈ c])] # All vectors with first norm
c2 = c[findall([isapprox(norms[2],norm(i),rtol=ε) for i ∈ c])] # All vectors with second norm
c3 = c[findall([isapprox(norms[3],norm(i),rtol=ε) for i ∈ c])] # All vectors with third norm
# Construct all possible bases, A′ (i.e., all combinations of c vectors), skip duplicate vectors
A′ = [[i j k] for i ∈ c1 for j ∈ c2 if !isapprox(i,j,rtol=ε) for k ∈ c3 if !isapprox(i,j,rtol=ε) && !isapprox(i,j,rtol=ε)]
A′ = A′[findall([isapprox(abs(det(i)),vol,rtol=ε) for i in A′])] # Delete candidate bases with the wrong volume
# The cross product is slightly (<1%) faster
#R = R[findall([abs(r[1]×r[2]⋅r[3])≈vol for r in R])] # Delete candidate bases with the wrong volume
A′ᵀ = [transpose(i) for i ∈ A′]
# This is the Uᵀ ̇U, where U transforms original basis to candidate basis
# If Tᵢ==identity then the U was a symmetry of the lattice
T = [Aiᵀ*A′ᵀ[i]*A′[i]*Ai for i ∈ 1:length(A′)]
# Indices of candidate T's that match the identity
idx = findall([all(norm.(t-I(3)) .< ε) for t ∈ T])
# Convert the transformations to integer matrices (formally they should be)
ops = [round.(Int,Ai*A′[i]) for i in idx] # Need the 'Int' so integers are returned
rops = [A′[i]*Ai for i in idx] 
