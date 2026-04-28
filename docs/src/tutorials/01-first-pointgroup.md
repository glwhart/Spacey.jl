# Find the point group of a cubic lattice

Welcome! In this 5 min. tutorial you'll find the **point group** of a 3D Bravais lattice — the set of integer rotations that map the lattice to itself. We'll work through two short examples: a simple cubic lattice (which has the maximum 48 symmetries) and a tetragonal lattice (which has fewer, because it's stretched along one axis).

By the end you'll know how to call `pointGroup` on a lattice you build by hand and read off the result. This tutorial assumes you have Spacey installed and a working Julia REPL. If `using Spacey` errors with "package not found," install it first:

```julia
julia> using Pkg

julia> Pkg.add(url="https://github.com/glwhart/Spacey.jl")
```

## 1. Set up

Start a Julia session and load Spacey, plus `LinearAlgebra` for the identity-matrix shorthand we'll use in a moment:

```jldoctest tut1
julia> using Spacey, LinearAlgebra
```

(The semicolon-free first line just imports the modules — it has no output to print.)

## 2. Define a simple cubic lattice

The simple cubic lattice has three perpendicular basis vectors of equal length. The cleanest way to spell this in Julia is the 3×3 identity matrix:

```jldoctest tut1
julia> A = Matrix{Float64}(I, 3, 3)
3×3 Matrix{Float64}:
 1.0  0.0  0.0
 0.0  1.0  0.0
 0.0  0.0  1.0
```

Each *column* of `A` is a basis vector — `(1, 0, 0)`, `(0, 1, 0)`, `(0, 0, 1)`. Spacey's column-as-basis-vector convention matches most crystallography codes, but it's worth being explicit because the row-vs-column choice is a common confusion.

If you prefer to specify the three vectors individually, that works too:

```jldoctest tut1
julia> u = [1.0, 0, 0]; v = [0, 1.0, 0]; w = [0, 0, 1.0];
```

Both forms work as input to `pointGroup`; we'll use both below.

## 3. Call `pointGroup`

`pointGroup` returns a tuple `(LG, G)`:

- `LG` is a `Vector{Matrix{Int}}` — each entry is one symmetry operation, expressed as an integer matrix in the lattice basis.
- `G` is the same operations, expressed as Cartesian rotations (`Vector{Matrix{Float64}}`).

For just the *order* of the group (the count of operations) we want `length(LG)`:

```jldoctest tut1
julia> LG, G = pointGroup(A);

julia> length(LG)
48

julia> length(pointGroup(u, v, w)[1])    # same answer using the three-vector form
48
```

**Forty-eight.** That's the maximum possible point group for a 3D Bravais lattice — the cubic point group, written O_h in [Schönflies notation](https://en.wikipedia.org/wiki/Sch%C3%B6nflies_notation). Every cubic lattice (simple cubic, FCC, BCC) has 48 lattice symmetries; the differences between them show up only when you add atoms (the next tutorial).

You can confirm the identity is one of the 48 ops:

```jldoctest tut1
julia> Matrix{Int}(I, 3, 3) in LG
true
```

Note that `LG[1]` is *not* guaranteed to be the identity. `pointGroup` returns the 48 ops in an unspecified order, so the identity might be at any index. (`spacegroup` is different — it does sort the identity to index 1.) If you need to retrieve a specific op, look it up by what the matrix is rather than by where it sits in the list — for the identity, the test `Matrix{Int}(I, 3, 3) in LG` above; for any other op, `findfirst(==(target), LG)`.

## 4. Try a tetragonal lattice

A *tetragonal* lattice is a cube stretched along one axis: still two equal-length perpendicular vectors in the plane, but the third has a different length. Set the c-axis to 1.5:

```jldoctest tut1
julia> A_tet = [1.0  0    0;
                0    1.0  0;
                0    0    1.5];

julia> length(pointGroup(A_tet)[1])
16
```

**Sixteen ops.** The `c ≠ a` distinction breaks every symmetry that mixes the c-axis with the a- or b-axes — what remains is the tetragonal point group [D_4h](https://en.wikipedia.org/wiki/Crystallographic_point_group), of order 16. As you stretch the cell further (`c = 2`, `c = 10`, …), the count *stays at 16*; it's only the moment `c` becomes distinguishable from `a` (within tolerance) that the count drops from 48.

To see the boundary in action, try `c = 1 + 1e-15`:

```jldoctest tut1
julia> A_almost_cubic = [1.0  0    0;
                         0    1.0  0;
                         0    0    1.0+1e-15];

julia> length(pointGroup(A_almost_cubic)[1])
48
```

At a 1e-15 distortion (well below `tol`), Spacey still sees this as cubic. That's the right answer for "is this lattice within float-precision of a cubic one." Try the same construction with `1 + 1e-3` and the answer drops to 16 — the structural distortion is now larger than the default `tol = 0.01`, so the algorithm rejects the cubic-only operations. This tolerance behavior is the subject of the [Detect tolerance-dependent answers](../how-to/detect-tolerance-dependence.md) how-to and the [Over-promotion](../explanation/over-promotion.md) explanation.

## What you learned

- `pointGroup(A)` accepts a 3×3 matrix or three basis vectors.
- It returns `(LG, G)` — integer-matrix and Cartesian forms of the same operations.
- A simple cubic lattice has 48 symmetries; a tetragonal lattice has 16.
- The number is intrinsic to the *lattice geometry* — atom positions don't enter at this level.

## Next steps

- Tutorial: [Find the space group of NaCl](02-first-spacegroup.md) — adds atoms to the cell and walks through four common materials.
- How-to: [Find a point group](../how-to/find-pointgroup.md) — tolerance handling, the `_fast` and `_simple` variants, rotation invariance.
- Explanation: [Why Minkowski reduction](../explanation/why-minkowski.md) — the theorem behind the algorithm.
