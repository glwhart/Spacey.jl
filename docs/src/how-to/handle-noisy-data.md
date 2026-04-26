# Handle noisy real-world data

Real lattices and atomic positions from experiment or simulation are never exact. This page is a working scientist's guide to choosing tolerances for [`pointGroup`](../reference/point-groups.md) and [`spacegroup`](../reference/space-groups.md) given the noise level in your input.

## Two tolerances, two scales

Spacey has two independent tolerances, each operating on a different quantity:

| Keyword | Used by | Operates on | Default |
|---|---|---|---|
| `tol`         | `pointGroup`, `spacegroup` (as `lattice_tol`) | (volume-normalized) lattice vectors | `0.01` |
| `pos_tol`     | `spacegroup`, `isSpacegroupOp`                | atomic positions in the user's lattice basis | `default_pos_tol(c)` ≈ 1% of `(V/N)^(1/3)` |

**Choose them independently.** Lattice noise (drift in unit-cell parameters) and position noise (drift in atomic coordinates) come from different sources and have different magnitudes. A DFT-relaxed structure typically has lattice parameters good to 4–5 significant figures but atomic positions good to only 3–4.

## Lattice tolerance: `tol`

`tol` is **relative** — it applies to the lattice after volume normalization, so it's unitless and the same `tol = 0.01` works for a 1-Å unit cell or a 100-Å one. Roughly, it's "what fractional deviation in normalized lattice vectors should still count as a symmetry."

Recommended starting points:

| Source of input | Suggested `tol` |
|---|---|
| Synthetic / hand-built / post-snap | `1e-6` |
| DFT-relaxed (well-converged) | `1e-4` |
| DFT-relaxed (loose convergence) | `1e-3` |
| Experimental refinement, good single-crystal | `1e-3` to `1e-2` |
| Experimental, powder or low quality | `1e-2` (the default) |

**Looser tolerance accepts more candidate operations and risks over-promotion** — silently reporting a higher symmetry than is actually present. The default `0.01` is intentionally on the tight side for that reason; the fix when you see a result that's "obviously wrong because it's missing the real symmetry" is to loosen, not the other way around.

```jldoctest
julia> using Spacey

julia> u = [1.0, 0, 0]; v = [0, 1.0, 0]; w = [0, 0, 1.0];   # exact cubic

julia> length(pointGroup(u, v, w; tol=1e-6)[1])             # tight, finds 48
48
```

## Position tolerance: `pos_tol`

`pos_tol` is **absolute**, in the same units as your lattice matrix `c.A`. The default is `0.01 · (V/N)^(1/3)` — 1% of the characteristic atom separation, computed from cell volume `V` and atom count `N`.

This formula scales correctly under any unit choice (Ångström, Bohr, nanometres) because both `V` and `(V/N)^(1/3)` rescale together.

Recommended starting points (assuming `c.A` is in Ångström):

| Crystal context | Suggested `pos_tol` |
|---|---|
| Synthetic / hand-built (positions exact) | `1e-6` Å |
| DFT-relaxed (positions good to ~0.001 Å) | `default_pos_tol(c)` (typically `~0.02–0.05` Å) |
| Experimental refinement | match to the reported atomic-displacement uncertainty |
| Material with known small distortion you want to detect (ferroelectrics, Jahn–Teller) | **smaller than the distortion** |

That last row is the trap. If you set `pos_tol` larger than the off-center displacement of the polar atom (e.g. Ti in BaTiO₃ at ~0.05 Å), Spacey reports the higher-symmetry parent (cubic, 48 ops) instead of the true ferroelectric (P4mm, 8 ops). It does so **silently**, with no warning, because at that tolerance the displaced atom does map to a parent-position image. Use [`verify_stable`](detect-tolerance-dependence.md) to detect this case automatically.

## Recipe: choose tolerances given known noise

1. **Estimate the noise** in your input. For DFT, look at the convergence threshold of the relaxation step. For experiment, look at the standard deviations on lattice parameters and atomic positions reported by the refinement. Pick a *single* noise scale per quantity (don't get fancy with anisotropic noise per axis at this level).

2. **Set `tol`** to roughly your lattice noise / `(V)^(1/3)` (i.e. relative to the characteristic length scale). For a 0.001 Å drift on a 4 Å cell, `tol ≈ 0.001/4 = 2.5e-4`. Round up to the nearest order of magnitude (`1e-3`) for safety.

3. **Set `pos_tol`** to roughly your position noise scale, in absolute units. For a 0.005 Å drift, `pos_tol = 0.005`. If you don't know, use `default_pos_tol(c)` and accept its ~1% heuristic.

4. **Run with `verify_stable=true`.** This re-runs at 1/1000 of your chosen tolerance and warns if the answer differs. See [Detect tolerance-dependent answers](detect-tolerance-dependence.md).

5. **If the warning fires**, you cannot resolve the symmetry from this input alone. See the discussion on the verify-stable page.

## Concrete example

A DFT-relaxed crystal with lattice noise ~5e-4 Å on a ~4 Å cell, position noise ~1e-3 Å:

```julia
using Spacey, LinearAlgebra

A = ...   # 3×3 lattice from your DFT output
r = ...   # 3×N fractional positions
c = Crystal(A, r, types; coords=:fractional)

ops = spacegroup(c;
                 lattice_tol = 1e-3,
                 pos_tol     = 1e-3,
                 verify_stable = true)
length(ops)
```

If `verify_stable` is silent at these tolerances, the result is trustworthy. If it warns, tighten and re-investigate — see [Detect tolerance-dependent answers](detect-tolerance-dependence.md).

## What `default_pos_tol` actually computes

```jldoctest
julia> using Spacey, LinearAlgebra

julia> A = Matrix{Float64}(I, 3, 3);                   # V = 1, N = 1

julia> c = Crystal(A, reshape([0.0, 0, 0], 3, 1), [:X]; coords=:fractional);

julia> default_pos_tol(c)
0.01
```

For a typical 4-Å cubic perovskite with 5 atoms per cell, `V/N = 64/5 ≈ 12.8 Å³`, `(V/N)^(1/3) ≈ 2.34 Å`, and `default_pos_tol ≈ 0.023 Å` — comfortably tighter than the BaTiO₃ Ti displacement (0.048 Å) and tighter still than typical Jahn–Teller distortions (0.05–0.2 Å). The default is engineered to **not** silently over-promote those classes of materials. See `designDiscussions.md` for the full reasoning behind the 1% factor.

## See also

- Reference: [`default_pos_tol`](../reference/crystals.md), [`pointGroup`](../reference/point-groups.md), [`spacegroup`](../reference/space-groups.md)
- How-to: [Detect tolerance-dependent answers](detect-tolerance-dependence.md), [Snap a noisy lattice to symmetry](snap-to-symmetry.md)
- Explanation: [Tolerances](../explanation/tolerances.md), [Over-promotion](../explanation/over-promotion.md)
