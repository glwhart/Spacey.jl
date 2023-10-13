using MinkowskiReduction
a1 = [1/16 - .0001, .0001, .0001]
a2 = [-0.001, 16-.0001, -.0001]
a3 = [-.0001, .0001, 1.51]
u,v,w = minkReduce(a1,a2,a3) # Always do this first, algorithm assumes reduced basis
#A = [a1 a2 a3]
a3 = [-.0000, .0001, 1.51]
a1 = [1+.01,0,0]; a2 = [0.,1-.01,0]; a3 = [0,0,1.5];
u,v,w=a1,a2,a3
#u,v,w = minkReduce(a1,a2,a3) # Always do this first, algorithm assumes reduced basis
#u = [1/16,0.,0.]; v=[0.,16.,0];w=[0.,0.,1.51];
A = [u v w] # Define a matrix with input vectors as columns
Ai = inv(A) 
Aiᵀ = transpose(Ai)
norms=norm.([u,v,w]) # Compute the norms of the three input vectors
#ε = 0.05min(norms...) # Scale factor for comparisons (unit tests must decide correct rescaling)
ε=.025
vol = abs(u×v⋅w) # Volume of the parallelipiped formed by the basis vectors
# A list of all possible lattice vectors in a rotated basis 
# These are lattice points from the vertices of the 8 cells with a corner at the origin)
# There are 27 of these (==3^3)
c = [A*[i,j,k] for i ∈ (-1,0,1) for j ∈ (-1,0,1) for k ∈ (-1,0,1)]
# Now keep only those vectors that have a norm matching one of the input vectors
# efficiency: Gather three groups, according to length. This limits the candidates even more
c1 = c[findall([isapprox(lens[1],norm(i),rtol=ε) for i ∈ c])] # All vectors with first norm
c2 = c[findall([isapprox(lens[2],norm(i),rtol=ε) for i ∈ c])] # All vectors with second norm
c3 = c[findall([isapprox(lens[3],norm(i),rtol=ε) for i ∈ c])] # All vectors with third norm
c1 = c[findall([isapprox(norms[1],norm(i),atol=ε) for i ∈ c])] # All vectors with first norm
c2 = c[findall([isapprox(norms[2],norm(i),atol=ε) for i ∈ c])] # All vectors with second norm
c3 = c[findall([isapprox(norms[3],norm(i),atol=ε) for i ∈ c])] # All vectors with third norm
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
norms = [norm(t-I(3)) for t ∈ T]
idx = findall([all(norm.(t-I(3)) .< ε) for t ∈ T])
idx = sortperm(norms)
# Convert the transformations to integer matrices (formally they should be)
ops = [round.(Int,Ai*A′[i]) for i in idx] # Need the 'Int' so integers are returned
rops = [A′[i]*Ai for i in idx] 

pointGroup_robust(a1,a2,a3)

isagroup(ops)

function makeGroup(ops)
    last = 1
    for i ∈ [2,3,4,6,8,12,16,24,48]
        println(i)
        if isagroup(ops[1:last])
            println("found groups")
            last = i 
        end
    end
    return ops[1:last]
end

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

A
ops
total = []
for i in [i*A[:,3] for i ∈ ops]
    if norm(i-A[:,3]) < .1
        println(round.(i,digits=3))
        append!(total,[i])   
    end
    end
sum(total)

b1,b3,b2=minkReduce(a1,a2,a3)
norm.([b1,b2,b3])
norm.([a1,a2,a3])

""" averageOverOps(vec,ops) 

Apply each operator in the list to the input vector. Take the average of the results than are approximately invariant
"""
function avgVecOverOps(vec,ops)
    cands = [iop*vec for iop ∈ ops]
    cands = filter(i->norm(i-vec)<.1*norm(vec),cands)
    avgVec = sum(cands)/length(cands)
    return avgVec
end

function snapToSymmetry_avg(v1,v2,v3,ops)
    w1 = avgVecOverOps(v1,ops)
    w2 = avgVecOverOps(v2,ops)
    w3 = avgVecOverOps(v3,ops)
    return w1,w2,w3
end

function snapToSymmetry_avg(M,ops)
    res = snapToSymmetry_avg(M[:,1],M[:,2],M[:,3],ops)
    return [res[1] res[2] res[3]]
end
#snapToSymmetry(u,v,w,ops)
A = [u v w] # Take the lattice basis as a matrix 
Ap = [A*k for k ∈ ops] # Convert the operators into Cartesian coordinates
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
     B[idx[2],idx[1]] = B[idx[1],idx[2]] # Symmetric matrix, copy elements across diagonal
end
s = svd(B)
Anew = diagm(sqrt.(s.S))*s.V'
T = A*inv(Anew)
t = svd(T)
Afinal = t.U*t.V'*Anew
u,v,w=[Afinal[:,i] for i ∈ 1:length(u)]
ops = pointGroup_robust(u,v,w)