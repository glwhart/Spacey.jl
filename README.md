# Spacey.jl
[![Runtests](https://github.com/glwhart/Spacey.jl/actions/workflows/Runtests.yml/badge.svg)](https://github.com/glwhart/Spacey.jl/actions/workflows/Runtests.yml)
[![codecov](https://codecov.io/gh/glwhart/Spacey.jl/branch/main/graph/badge.svg?token=OYAPFZI28I)](https://codecov.io/gh/glwhart/Spacey.jl)

The package provides functions related to finding the spacegroup of a crystal. (The spacegroup of a crystal is useful in the context density functional theory calculations, computational materials science, and machine learning for materials prediction.)

This package will eventually provide all of the functionality of the [enumlib package](https://github.com/msg-byu/enumlib) and some additional functionality, useful for machine learning applications and materials science. In this pre-release version, only a pointgroup finder is provided, as well as a "snap to symmetry" function for unit cells.

One of the unique features of the package is its robustness to finite-precision issues, a constant bane of other spacegroup codes.

A few simple examples are provided below. This package finds symmetries very efficiently* by relying on the assumption that the basis is as compact as possible. To assure the most compact basis, this code first reduces the input basis vectors to a "Minkowski-reduced" basis. (See [MinkowskiReduce.jl](https://github.com/glwhart/MinkowskiReduction.jl))

(* Normally, efficiency is not important, since finding the symmmetryies is a relatively quick computation. When the symmetries of tens of thousands of cases are needed in just a second or two, as in [autoGR](https://github.com/msg-byu/autoGR), efficiency becomes essential.)


## Example 1: Hexagonal lattice
```
julia> u,v,w = ([1, 0, 0], [0.5, √3/2, 0.0], [0.0, 0.0, √(8/3)])
([1, 0, 0], [0.5, 0.8660254037844386, 0.0], [0.0, 0.0, 1.632993161855452])

julia> pointGroup(u,v,w)
24-element Vector{Matrix{Int64}}:
 [-1 0 0; -1 1 0; 0 0 -1]
...
 [1 0 0; 1 -1 0; 0 0 1]
```

By default the output is is "direct coordinates", the transformation matrices showing integer combinations of the basis that yield equivalent, but rotated basis vectors.  

## Example 2: Rhombohedral lattice, arbitrary rotation of coordinate system

```
 julia> u = [1, 1, 2]; v = [1, 2, 1];  w = [2, 1, 1];

julia> u, v, w = threeDrotation(u, v, w, π / 3, π / 5, π / 7)
([1.0328466951537383, 1.855345731770484, 1.2210323173082058], [-0.2604424479658839, 2.2850084410500355, 0.8431525103014423], [0.7968127453125833, 2.3142904165695555, -0.09565206052008146])

julia> length(pointGroup(u, v, w)) == 12
true

julia> pointGroup(u, v, w)
12-element Vector{Matrix{Int64}}:
 [-1 0 0; -1 1 1; 0 0 -1]
 [-1 1 0; -1 0 -1; 0 0 1]
 [-1 0 0; 0 -1 0; 0 0 -1]
 [-1 1 0; 0 1 0; 0 0 1]
 [0 -1 -1; -1 0 -1; 0 0 1]
 [0 1 1; -1 1 1; 0 0 -1]
 [0 -1 -1; 1 -1 -1; 0 0 1]
 [0 1 1; 1 0 1; 0 0 -1]
 [1 -1 0; 0 -1 0; 0 0 -1]
 [1 0 0; 0 1 0; 0 0 1]
 [1 -1 0; 1 0 1; 0 0 -1]
 [1 0 0; 1 -1 -1; 0 0 1]
 ```

 ## Example 3: Slightly distorted cubic case, "snap back" to perfect cubic cell
The lattice vectors are not quite orthogonal, not quite all the same length. First `pointGroup_robust` gets the symmetries that _should_ be present for a perfect cubic cell. `snapToSymmetry` then finds the perfect cubic cell that is _as close as possible_ to the original cell.

```
julia> u = [1+.01,0,0]; v = [0.,1-.01,0]; w = [0,0,1-.001];

julia> ops = pointGroup_robust(u,v,w)
48-element Vector{Matrix{Int64}}:
 [-1 0 0; 0 -1 0; 0 0 -1]
 [-1 0 0; 0 -1 0; 0 0 1]
 [-1 0 0; 0 0 -1; 0 -1 0]
 [-1 0 0; 0 0 1; 0 -1 0]
 [-1 0 0; 0 0 -1; 0 1 0]
 [-1 0 0; 0 0 1; 0 1 0]
 [-1 0 0; 0 1 0; 0 0 -1]
 [-1 0 0; 0 1 0; 0 0 1]
 [0 -1 0; -1 0 0; 0 0 -1]
 ⋮
 [1 0 0; 0 -1 0; 0 0 -1]
 [1 0 0; 0 -1 0; 0 0 1]
 [1 0 0; 0 0 -1; 0 -1 0]
 [1 0 0; 0 0 1; 0 -1 0]
 [1 0 0; 0 0 -1; 0 1 0]
 [1 0 0; 0 0 1; 0 1 0]
 [1 0 0; 0 1 0; 0 0 -1]
 [1 0 0; 0 1 0; 0 0 1]
a,b,c,newops = snapToSymmetry(u,v,w,ops)
([0.9996332321644671, -1.1098158305294557e-16, 3.329447491588367e-16], [0.0, 0.9996332321644675, -3.330557677480862e-16], [0.0, -3.330557677480862e-16, 0.9996332321644676], [[-1 0 0; 0 -1 0; 0 0 -1], [-1 0 0; 0 -1 0; 0 0 1], [-1 0 0; 0 0 -1; 0 -1 0], [-1 0 0; 0 0 1; 0 -1 0], [-1 0 0; 0 0 -1; 0 1 0], [-1 0 0; 0 0 1; 0 1 0], [-1 0 0; 0 1 0; 0 0 -1], [-1 0 0; 0 1 0; 0 0 1], [0 -1 0; -1 0 0; 0 0 -1], [0 -1 0; -1 0 0; 0 0 1]  …  [0 1 0; 1 0 0; 0 0 -1], [0 1 0; 1 0 0; 0 0 1], [1 0 0; 0 -1 0; 0 0 -1], [1 0 0; 0 -1 0; 0 0 1], [1 0 0; 0 0 -1; 0 -1 0], [1 0 0; 0 0 1; 0 -1 0], [1 0 0; 0 0 -1; 0 1 0], [1 0 0; 0 0 1; 0 1 0], [1 0 0; 0 1 0; 0 0 -1], [1 0 0; 0 1 0; 0 0 1]])
```
