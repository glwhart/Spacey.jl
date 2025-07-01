using Spacey
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
using MinkowskiReduction

A = [0.0 0.5 0.5; 0.5 0.0 0.5; 0.5 0.5 0.0]
ε = 2e-3
@time for i ∈ 1:1000
    noise = (2*rand(3,3).-1)*ε
    Atemp = A + noise
    length(pointGroup_robust(minkReduce(eachcol(Atemp)...)[1:3]...)[1])!=48 && error("Point group is not 48")
end 

# Test case 2
begin
maxε = 1e-6
for ε ∈ logrange(2*maxε,1e-3,10)
ar = 1.0
maxε = 1e-6
for ε ∈ logrange(2*maxε,1e-3,10)
ar = 1.0
aspect_ratio = [1 0 0; 0 1 0; 0 0 ar]
A = aspect_ratio*[0.0 0.5 0.5; 0.5 0.0 0.5; 0.5 0.5 0.0]
for i ∈ 1:50
    noise = (2*rand(3,3).-1)*ε
    Atemp = A + noise
    if length(pointGroup_robust(minkReduce(eachcol(Atemp)...)[1:3]...)[1])!=16
        println("ε: ", ε)
        maxε = ε
        break
    end
    if maxε == ε; break; end
end 
println("Max ε: ", maxε)
end
end
end


## FCT test case
begin
p1=plot()
tols = logrange(5e-4,5e-1,10)        # Tolerance values to test
#colors = distinguishable_colors(length(tols))  # Generate clearly distinct colors
colors = palette(:viridis, length(tols))       # Use the viridis colour scheme
for (idx,tol) ∈ enumerate(tols)
     println("tol: ", round(tol,digits=5))
     A = [0.0 0.5 0.5; 0.5 0.0 0.5; 0.52 0.52 0.0]
    #A = [0.0 0.5 0.5; 0.5 0.0 0.5; 0.5 0.5 0.0]
    #A =[-1.0 1.0 1.0; 1.0 -1.0 1.0; 1.0 1.0 -1.0]
     a = 5e-30; Navg =100; Nsteps = 40;
     data = Vector{Float64}(undef,Nsteps)
     plim = logrange(5e-4,1e0,Nsteps)
     for (i,ε) ∈ enumerate(plim)
         data[i] = count([length(pointGroup_robust(minkReduce(eachcol((A+(2*rand(3,3).-1)*ε)*a)...)[1:3]...;tol=tol)[1])==48 for _ ∈ 1:Navg])/Navg
     end
p1=plot!(plim[findall(data.>0)],
         data[findall(data.>0)];
         yscale=:log10,
         xscale=:log10,
         xlabel="Noise level (ε/a)",
         ylabel="Success rate",
         title="BCC case, Unscaled tol ",
         label=string(round(tol,digits=5)),
         lw=3,
         color=colors[idx],
         legend=:bottomleft,
         #xticks=[.005,1e-2,1e-1,1e0,1e1,5e1],
         )
 end
 display(p1)
 end
# Presumably the tol setting in pointGroup_robust can be as much as 10% of the smallest lattice vector and we'll get lots of candidates an the symmetry finder will be slow but more robust.

for i ∈ -30:2:30
    for _ ∈ 1:200
        if !all([length(pointGroup(minkReduce(A+ (2*rand(3,3).-1)*0.02)*10.0^i;tol=1e-10)[1])==48]) error("Symmetry group is not 48") end
    end
end


length(pointGroup(minkReduce(A+ (2*rand(3,3).-1)*0.02)*1e-14)[1])==48
begin
Atest = (A+ (2*rand(3,3).-1)*0.02)*1e-15
pointGroup(minkReduce(Atest))[1]
end

begin
a = 5e-3; maxε = 0.0
for ε ∈ logrange(4e-3*a,6e-3*a,20)
    maxε = ε
    success = true
    for i ∈ 1:1000
        noise = (2*rand(3,3).-1)*ε*a
        Atemp = A*a + noise
        #@show Atemp
        nops = length(pointGroup_robust(minkReduce(eachcol(Atemp)...)[1:3]...)[1])
        if nops != 48
            println("Symmetry group is not 48. ε: ", round(ε/a,digits=5), "   nops: ", nops, "    noise: ", round(ε/a,digits=5))
            success = false
            break
        end
    end 
    if !success; maxε = ε; break; end
    println("ε/a ", round(ε/a,digits=5), " passed")
end
println("Max ε: ", round(maxε,digits=6))
end








