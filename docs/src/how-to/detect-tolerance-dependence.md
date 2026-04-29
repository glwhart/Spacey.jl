# Detect tolerance-dependent answers

Both [`pointGroup`](../reference/point-groups.md) and [`spacegroup`](../reference/space-groups.md) accept a `verify_stable=true` keyword. When set, the algorithm runs twice — once at the requested tolerance, once again at 1/1000 of it — and emits a `@warn` if the operation count differs between the two runs. A difference means the answer depends on the tolerance: the lattice (or atomic positions) is near a symmetry boundary.

This setting provides a guard against silent over-promotion (finding a symmetry that is too high). Use it when you cannot personally vouch for the input's noise level.

## When to use it

- **Always**, on first analysis of unfamiliar input — experimental refinements, new structures from a database, output of someone else's relaxation pipeline.
- **Periodically**, in a CI loop over a corpus of structures, as a regression detector.
- **Skip it** when the input is exact / synthetic and you control its construction — for example, when you build a known Bravais prototype to validate downstream code.

The cost is one extra invocation at tighter tolerance — the same algorithmic cost as the first call. Unless the routine is being called for tens of thousands of structures (say in automatic k-point generation) then the marginal cost of being careful is likely worth it.

## Pattern: point group

```jldoctest
julia> using Spacey

julia> u = [1.0, 0, 0]; v = [0, 1.0, 0]; w = [0, 0, 1.0];   # exact cubic

julia> length(pointGroup(u, v, w; verify_stable=true))   # silent: 48 at all tols
48
```

When the lattice is near a symmetry boundary (e.g., tetragonal that's *almost* cubic at a loose tolerance), the warning fires:

```julia
julia> ε = 1e-3;

julia> u = [1.0, 0, 0]; v = [0, 1.0, 0]; w = [0, 0, 1+ε];   # tetragonal, 16 ops at tight tol

julia> # at loose tol, this looks cubic (48 ops) — the warning catches the disagreement
julia> pointGroup(u, v, w; tol=0.1, verify_stable=true);
┌ Warning: pointGroup_robust: group size depends on tolerance — lattice is near a symmetry boundary.
│   tol = 0.1
│   group_at_tol = 48
│   tight_tol = 0.0001
│   group_at_tight_tol = 16
└ @ Spacey ...
```

The returned group is the one at the *requested* `tol` — `verify_stable` does not change the answer, only flags the fact that the answer is sensitive to the tolerance.

## Pattern: space group

The same flag exists on `spacegroup`, where it re-runs at `pos_tol / 1000`. This catches **position over-promotion** — when atomic displacements are below `pos_tol`, Spacey treats the crystal as the higher-symmetry parent and the answer flips depending on how `pos_tol` is set.

The canonical case is a ferroelectric like BaTiO₃: the cubic perovskite Pm3̄m parent has Ba at the corner, Ti at the body center, and three O at the face centers. Below the Curie point the Ti shifts by ~0.05 Å along z, dropping the symmetry to P4mm (8 ops). At loose `pos_tol`, Spacey reports the parent Pm3̄m (48 ops):

```julia
using Spacey, LinearAlgebra

A = 4.0 * Matrix{Float64}(I, 3, 3)              # cubic lattice, 4 Å
ε = 0.05                                         # Ti displacement in Å
δ = ε / 4.0                                      # in fractional coords
r = [0.0  0.5    0.5  0.5  0.0;                  # Ba  Ti      O    O    O
     0.0  0.5    0.5  0.0  0.5;
     0.0  0.5-δ  0.0  0.5  0.5]
c_BaTiO3 = Crystal(A, r, [:Ba, :Ti, :O, :O, :O]; coords=:fractional)

# At loose pos_tol (0.1 Å > 0.05 Å), the displaced Ti maps onto the cubic
# body-center image — Spacey reports the parent cubic group, silently.
length(spacegroup(c_BaTiO3; pos_tol=0.1))   # 48 ops (over-promoted)

# verify_stable re-runs at pos_tol/1000 = 1e-4 Å, where the displacement
# IS resolved, sees only 8 ops, and warns about the disagreement.
spacegroup(c_BaTiO3; pos_tol=0.1, verify_stable=true);
# ┌ Warning: spacegroup: operation count depends on pos_tol — crystal is near a position-symmetry boundary.
# │   pos_tol = 0.1
# │   ops_at_pos_tol = 48
# │   tight_pos_tol = 0.0001
# │   ops_at_tight_pos_tol = 8
# └ @ Spacey ...
```

For an exact, non-degenerate structure, `verify_stable` stays silent — the operation count doesn't change as `pos_tol` shrinks. Zincblende GaAs in its primitive 2-atom cell (Ga at the origin, As at ¼,¼,¼) is a clean example:

```jldoctest
julia> using Spacey

julia> A = [0.0  0.5  0.5;        # FCC primitive lattice vectors
            0.5  0.0  0.5;
            0.5  0.5  0.0];

julia> r = [0.0  0.25;             # Ga at origin, As at (¼,¼,¼)
            0.0  0.25;
            0.0  0.25];

julia> c = Crystal(A, r, [:Ga, :As]; coords=:fractional);

julia> length(spacegroup(c; verify_stable=true))   # 24 ops (T_d), no warning
24
```

## What to do when the warning fires

1. **Don't ignore it.** Silent over-promotion is the canonical class of bug in symmetry analysis — the wrong space group propagates into downstream calculations.
2. **Investigate the input.** What's the noise level you actually expect? Is the structure genuinely near a parent symmetry (a real ferroelectric, a displacive transition near its high-temperature parent, an octahedral tilt in a perovskite), or is noise causing a broken symmetry when it shouldn't?
3. **Tighten the tolerance** to one that's smaller than the structural distortion you want to detect, and re-run with `verify_stable`. If the answer is now stable, that's the truth.
4. **If you cannot tighten** (genuine experimental noise that exceeds the displacement of interest), the question is unanswerable from this input alone — it's a science question. Consider [snapping to symmetry](snap-to-symmetry.md) and re-analyzing, or accept the higher-symmetry answer with the caveat documented in your output.

## Cost

`verify_stable` effectively doubles the runtime of the call. For a 100-atom crystal, milliseconds become a few milliseconds; for a 1000-atom crystal, seconds become a few seconds. For corpus-scale processing (thousands of structures), keep it on for the first pass and turn it off only after the corpus is verified clean.

## See also

- Reference: [`pointGroup`](../reference/point-groups.md), [`spacegroup`](../reference/space-groups.md)
- How-to: [Find a point group](find-pointgroup.md), [Find a space group](find-spacegroup.md), [Handle noisy real-world data](handle-noisy-data.md)
- Explanation: [Over-promotion](../explanation/over-promotion.md)
