using MinkowskiReduction
using LinearAlgebra

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

# Test case
a1 = [1/16 - .0001, .001, .0001]
a2 = [-0.001, 16-.0001, -.0001]
a3 = [-.0001, .0001, 1.51]
ops, rops = pointGroup_robust(minkReduce(a1,a2,a3)...)

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

using Spacey
A = [0.0 0.5 0.5; 0.5 0.0 0.5; 0.5 0.5 0.0]
pointGroup(A)