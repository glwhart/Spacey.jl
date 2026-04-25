# Classify the Bravais system

[`crystal_system`](../reference/crystals.md) returns the Bravais-system symbol of a lattice — one of seven values:

| Symbol | Crystal system | Holohedry order |
|---|---|---|
| `:triclinic`    | triclinic    | 2  |
| `:monoclinic`   | monoclinic   | 4  |
| `:orthorhombic` | orthorhombic | 8  |
| `:trigonal`     | trigonal     | 12 |
| `:tetragonal`   | tetragonal   | 16 |
| `:hexagonal`    | hexagonal    | 24 |
| `:cubic`        | cubic        | 48 |

Identification uses the order of the lattice's point group (its holohedry), which uniquely determines the system. This is the **7-class** (Layer 1) classification — for the full **14-class** Bravais lattice classification (which adds centering: P/I/F/C/R), see [Crystal system vs full Bravais](../explanation/crystal-system-vs-bravais.md).

## From a 3×3 lattice matrix

```jldoctest
julia> using Spacey, LinearAlgebra

julia> crystal_system(Matrix{Float64}(I, 3, 3))
:cubic

julia> crystal_system([1.0 0 0; 0 1 0; 0 0 1.5])
:tetragonal

julia> crystal_system([1.0 0 0; 0 1.2 0; 0 0 1.5])
:orthorhombic

julia> crystal_system([1.0 0.3 0.2; 0 1.1 0.4; 0 0 0.9])
:triclinic
```

## From a `Crystal`

Same answer — `crystal_system(c)` reports the lattice's symmetry, not the crystal's full space group.

```jldoctest
julia> using Spacey, LinearAlgebra

julia> A = Matrix{Float64}(I, 3, 3);

julia> r = reshape([0.0, 0.0, 0.0], 3, 1);

julia> c = Crystal(A, r, [:X]; coords=:fractional);

julia> crystal_system(c)
:cubic
```

## Tune `lattice_tol`

`crystal_system` accepts the same `lattice_tol` keyword as [`pointGroup_robust`](../reference/point-groups.md), which it calls internally. Default is `0.01`. Loosen for noisy lattices; tighten if you suspect over-promotion (see [Detect tolerance-dependent answers](detect-tolerance-dependence.md)).

```jldoctest
julia> using Spacey

julia> A_noisy = [1.0 0 0; 0 1.0+1e-4 0; 0 0 1.0];   # near-cubic, slightly tetragonal

julia> crystal_system(A_noisy; lattice_tol=1e-3)   # at this tol, looks cubic
:cubic

julia> crystal_system(A_noisy; lattice_tol=1e-6)   # tight: actually tetragonal
:tetragonal
```

This is the canonical over-promotion footgun — a tolerance that's looser than the structural distortion silently classifies a tetragonal lattice as cubic. The size of the lattice distortion that causes this is `~lattice_tol` in normalized units.

## See also

- Reference: [`crystal_system`](../reference/crystals.md)
- Explanation: [Crystal system vs full Bravais](../explanation/crystal-system-vs-bravais.md), [Over-promotion](../explanation/over-promotion.md)
