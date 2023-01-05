using Spacey

u = [1/16 - .0001, .001, .0001]
v = [-0.001, 16-.0001, -.0001]
w = [-.0001, .0001, 1.51]

# This is appoximately an orthorhombic cell, should have 8 ops
# Fast gives 2 (it has a stringent ϵ)
ops=pointGroup_fast(u,v,w)

# Robust version gives 24 (it has a generous ϵ)
# Presumably 3 copies of the 8 correct operations
ops,rops=pointGroup_robust(u,v,w)

snapToSymmetry(u,v,w,ops)

u = [8.581165,     -0.067727,      0.063365]
v = [-0.068580,      8.576216,     -0.068580]
w = [0.063365,     -0.067727,      8.581165]

ops=pointGroup_fast(u,v,w)