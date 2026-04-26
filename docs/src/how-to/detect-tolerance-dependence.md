# Detect tolerance-dependent answers

Both [`pointGroup`](../reference/point-groups.md) and [`spacegroup`](../reference/space-groups.md) accept a `verify_stable=true` keyword. When set, the algorithm runs twice — once at the requested tolerance, once at 1/1000 of it — and emits a `@warn` if the operation count differs between the two runs. A difference means the answer depends on the tolerance: the lattice (or atomic positions) is near a symmetry boundary.

This is the canonical guard against silent over-promotion. Use it when you cannot personally vouch for the input's noise level.

## When to use it

- **Always**, on first analysis of unfamiliar input — experimental refinements, new structures from a database, output of someone else's relaxation pipeline.
- **Periodically**, in a CI loop over a corpus of structures, as a regression detector.
- **Skip** when the input is exact / synthetic and you control its construction.

The cost is one extra invocation at tighter tolerance — the same algorithmic cost as the first call. Worth it.

## Pattern: point group

```jldoctest
julia> using Spacey

julia> u = [1.0, 0, 0]; v = [0, 1.0, 0]; w = [0, 0, 1.0];   # exact cubic

julia> length(pointGroup(u, v, w; verify_stable=true)[1])   # silent: 48 at all tols
48
```

When the lattice is near a symmetry boundary (e.g. tetragonal that's *almost* cubic at a loose tolerance), the warning fires:

```julia
julia> ε = 1e-3;

julia> u = [1.0, 0, 0]; v = [0, 1.0, 0]; w = [0, 0, 1+ε];   # tetragonal, 16 ops at tight tol

julia> # at loose tol, this looks cubic (48 ops) — the warning catches the disagreement
julia> pointGroup(u, v, w; tol=0.1, verify_stable=true);
┌ Warning: pointGroup: result depends on `tol` — lattice is near a symmetry boundary.
│   tol = 0.1
│   tight_tol = 0.0001
│   nops_at_tol = 48
│   nops_at_tight_tol = 16
└ @ Spacey ...
```

The returned group is the one at the *requested* `tol` — `verify_stable` does not change the answer, only flags it.

## Pattern: space group

The same flag exists on `spacegroup`, where it re-runs at `pos_tol / 1000`. This catches **position over-promotion** — the canonical case is a ferroelectric like BaTiO₃ where atoms are slightly displaced from the cubic prototype:

```jldoctest
julia> using Spacey, LinearAlgebra

julia> A = Matrix{Float64}(I, 3, 3);

julia> c = Crystal(A, reshape([0.0, 0.0, 0.0], 3, 1), [:X]; coords=:fractional);

julia> length(spacegroup(c; verify_stable=true))   # silent: 48 ops, exact cubic
48
```

When atomic displacements are below `pos_tol`, Spacey treats the crystal as the higher-symmetry parent and `verify_stable` flags the disagreement:

```julia
# A cubic lattice with Ti displaced from body-center by 0.05 Å (BaTiO₃-like).
# True space group is P4mm (8 ops); at loose pos_tol it looks cubic (48).
julia> spacegroup(c_BaTiO3; pos_tol=0.1, verify_stable=true);
┌ Warning: spacegroup: operation count depends on pos_tol — crystal is near a position-symmetry boundary.
│   pos_tol = 0.1
│   ops_at_pos_tol = 48
│   tight_pos_tol = 0.0001
│   ops_at_tight_pos_tol = 8
└ @ Spacey ...
```

## What to do when the warning fires

1. **Don't ignore it.** Silent over-promotion is the canonical class of bug in symmetry analysis — the wrong space group propagates into every downstream calculation.
2. **Investigate the input.** What's the noise level you actually expect? Is the structure genuinely near a parent symmetry (a real ferroelectric, an Andersen-class transition), or is the noise dominating signal?
3. **Tighten the tolerance** to one that's smaller than the structural distortion you want to detect, and re-run without `verify_stable`. If the answer is now stable, that's the truth.
4. **If you cannot tighten** (genuine experimental noise that exceeds the displacement of interest), the question is unanswerable from this input alone. Consider [snapping to symmetry](snap-to-symmetry.md) and re-analyzing, or accept the higher-symmetry answer with the caveat documented in your output.

## Cost

`verify_stable` doubles the runtime of the call. For point-group analysis on a single lattice this is microseconds → microseconds; for space-group analysis on a 100-atom crystal it goes from O(N³) to 2·O(N³). For corpus-scale processing (thousands of structures), keep it on for the first pass and turn it off only after the corpus is verified clean.

## See also

- Reference: [`pointGroup`](../reference/point-groups.md), [`spacegroup`](../reference/space-groups.md)
- How-to: [Find a point group](find-pointgroup.md), [Find a space group](find-spacegroup.md), [Handle noisy real-world data](handle-noisy-data.md)
- Explanation: [Over-promotion](../explanation/over-promotion.md)
