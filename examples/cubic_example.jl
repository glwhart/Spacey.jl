using MinkowskiReduction
using Spacey
using LinearAlgebra
# distorted simple cubic
a1 = [1+.01,0,0]; a2 = [0.,1-.01,0]; a3 = [0,0,1-.001];
u,v,w = minkReduce(a1,a2,a3)
ops = pointGroup_robust(u,v,w)
a,b,c = snapToSymmetry(u,v,w,ops)
det([a b c])≈det([a1 a2 a3])

a1 = [1+.01,0,0]; a2 = [0.,1-.01,0]; a3 = [0,0,1.5];
u,v,w = minkReduce(a1,a2,a3)
ops = pointGroup_robust(u,v,w)
a,b,c,ops = snapToSymmetry(u,v,w,ops)
println("Num ops: ",length(ops))
norm(a)≈norm(b)
println(det([u v w]))
println(det([a1 a2 a3]))
println(det([a b c]))
A=[a1 a2 a3]
pg, rpg = pointGroup(A)