using LinearAlgebra, Plots
u,v = ([1;0],[.125;1.6]) # starting lattice vectors, peturbed from 
A = [u v] # Matrix of basis vectors
iA = inv(A) # Inverse of basis vectors, convert from lattice to Cartesian
B = iA*transpose(iA) # Useful to check for orthogonality

# Get the candidate lattice vectors, vectors that might be rotated versions of the starting basis
cLV = filter(x->isapprox(norm(x),1.0; atol=0.2)||isapprox(norm(x),1.6;atol=.1),[A*[i;j] for i ∈ -1:1 for j ∈ -1:1])

# Candidate lattice vectors, in lattice coordinates
cL = filter(x->isapprox(norm(A*x),1.0; atol=0.2)||isapprox(norm(A*x),1.6;atol=.1),[[i;j] for i ∈ -1:1 for j ∈ -1:1])

# Candidate bases that preserve the original volume
RL = filter(x->isapprox(abs(det(x)),1.0; atol=0.2),[[i j] for i ∈ cL for j ∈ cL if i ≠ j])
A′ = [A*k for k ∈ RL]
R = [p*iA for p ∈ A′]
# Candidate Rs

#Ut = [x*B*transpose(x) for x ∈ R]
Ut =  [r*transpose(r) for r ∈ R]
diffs = [x.-I(2) for x ∈ Ut]
idx = findall(norm.(diffs).<.3)
R = R[idx]
s = [svd(x) for x ∈ R]
Ro = [i.U*transpose(i.V) for i ∈ s]

#cR=filter(x->norm(x*B*transpose(x).-I(2))<√eps(1.0),RR)
#cR = [i*iA for i ∈ cR]
#L=filter(x->norm(A*x*B*transpose(A*x).-I(2))<√eps(1.0),RL)

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
# u
#ts = reduce(hcat,[A*ceil.(iA*(r*u-u)) for r ∈ Ro])
# pts = reduce(hcat,[rotate_and_shift_back(u,r,A,iA) for r ∈ Ro])
# scatter(pts[1,:],pts[2,:],aspect_ratio=1,legend=:none)
#plot!(A[1,:],A[2,:])

pts = reduce(hcat,[rotate_and_shift_back(v,r,A,iA) for r ∈ Ro])
pts = reduce(hcat,[r*A for r ∈ Ro])
scatter(pts[1,:],pts[2,:],aspect_ratio=1,legend=:none)


# git branch -m master main
# git fetch origin
# git branch -u origin/main main
# git remote set-head origin -aspect_ratio