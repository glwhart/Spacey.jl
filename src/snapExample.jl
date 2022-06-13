using Spacey

u = [1/16 - .0001, .001, .0001]
v = [-0.001, 16-.0001, -.0001]
w = [-.0001, .0001, 1.51]

# This is appoximately and orthorhombic cell, should have 8 ops
# Fast gives 2 (it has a stringent ϵ)
ops=pointGroup_fast(u,v,w)

# Robust version gives 24 (it has a generous ϵ)
ops=pointGroup_robust(u,v,w)

snapToSymmetry(u,v,w,ops)
