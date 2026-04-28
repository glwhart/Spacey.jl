# Algorithm overview

This page walks through what Spacey actually *does* when you call `pointGroup` or `spacegroup`. It's the conceptual pipeline — not the code line-by-line, but enough that you can understand each step's purpose and what could go wrong.

The point-group algorithm is the foundation; `spacegroup` builds on it. We'll cover the point-group pipeline first, then say what changes for crystals.

## The point-group pipeline

Given a 3×3 lattice matrix `A` (columns are basis vectors), `pointGroup` returns the integer-matrix and Cartesian forms of every symmetry operation of the lattice, in five steps:

1. **Minkowski-reduce** the input basis. Replaces `A` with an equivalent basis whose vectors are as short as possible. See [Why Minkowski reduction](why-minkowski.md) for the theorem this depends on.
2. **Generate** the 27 candidate vectors `A · [i, j, k]` for `i, j, k ∈ {-1, 0, 1}`. The 27-neighbor result guarantees that every lattice symmetry maps each basis vector into this set.
3. **Filter** candidate symmetry operations through three cheap tests:
    - **Norm match**: the candidate basis must have the same vector lengths as the original (in some order).
    - **Volume conservation**: the candidate basis must have the same volume.
    - **Integer-matrix test**: `U = inv(A) · B` (where `B` is the candidate basis) must be an integer matrix to within `tol`.
4. **Close** under multiplication: from the operations that survived step 3, find the largest subset that forms a group. Spacey tries the canonical group sizes `[48, 24, 16, 12, 8, 4, 2]` in decreasing order and accepts the largest one that closes.
5. **(Optional) Verify stability**: if `verify_stable=true`, re-run the entire pipeline at `tol/1000` and emit a `@warn` if the resulting group order differs.

The rest of this page is what each step is doing and why.

## Step 1: Minkowski reduction

Calls into [MinkowskiReduction.jl](https://github.com/glwhart/MinkowskiReduction.jl). The output basis satisfies the 12 Minkowski conditions: each basis vector is the shortest possible given the previous ones, and angles are bounded in `[60°, 120°]`. This is what makes the search in step 2 *provably finite and complete* — see [Why Minkowski reduction](why-minkowski.md).

If you call `pointGroup` directly with a non-reduced basis, the function errors out and asks you to reduce first. (`Spacey.pointGroup_robust` is the same: explicit Mink-reduction is required.) The matrix-form wrapper handles this transparently: it Mink-reduces the input columns automatically.

## Step 2: 27-candidate generation

Each basis vector image must be a length-preserved lattice vector. The 27-neighbor theorem says these are exactly the integer combinations `A · [i, j, k]` with `i, j, k ∈ {-1, 0, 1}` — twenty-seven vectors, no more.

In code:

```julia
c = [A * [i, j, k] for i in (-1, 0, 1) for j in (-1, 0, 1) for k in (-1, 0, 1)]
```

A candidate symmetry operation chooses three of these vectors as the new basis. Naively that's `27³ = 19,683` candidate triples. The next step prunes that to dozens.

## Step 3: Three-stage filtering

Most candidate triples are obviously not symmetries. Three filters cut through the space cheaply:

- **Norm match.** A symmetry sends each basis vector to a vector of the *same length*. Group the 27 candidates by their norm; only triples that match the original lengths (in some order) survive. For a cubic lattice with all three basis vectors of equal length, this leaves a few hundred triples; for a triclinic lattice with three distinct lengths, much fewer.
- **Volume conservation.** A symmetry preserves the determinant: `|det(B)| = |det(A)|`. Checking this kills another large fraction.
- **Integer-matrix test.** This is the heart of the algorithm. The candidate `B` represents a symmetry operation iff `U = A⁻¹ · B` has integer entries (the symmetry is then a lattice transformation, expressible exactly in integer coordinates). For a *clean* input we can require `U` to be exactly integer; for a *noisy* input we require `U` to be close to integer, with closeness controlled by `tol`. Specifically: each entry must be within `tol` of the nearest integer.

Before any of these filters runs, Spacey **volume-normalizes the input**: divides the basis vectors by `∛|det(A)|` so that the cell has unit volume. This pulls every comparison onto a canonical scale — `tol = 0.01` means the same thing whether your lattice is in Ångström or in Bohr. See [Tolerances](tolerances.md) for the full discussion.

After all three filters, you have a set of candidate operations — typically a few hundred for clean cubic input, fewer for lower-symmetry lattices. Most of them are real symmetries; some may be near-misses that passed the integer-matrix test by a small margin.

## Step 4: Group closure

The set of candidates from step 3 is *almost* the point group, but it can contain spurious near-misses (an operation that passed the integer test but doesn't actually compose properly with the others). The way to clean this up is to require that the operations form a *group* — closed under multiplication.

A point group of a 3D Bravais lattice has order `48`, `24`, `16`, `12`, `8`, `4`, or `2` (per the seven crystal systems: cubic / hexagonal / tetragonal / trigonal / orthorhombic / monoclinic / triclinic). Spacey tries each size from largest to smallest, looking for the *largest* subset of candidates that closes under composition. The first size that succeeds is reported.

This step is what makes Spacey robust to spurious near-misses on its own: a spurious operation that passes the integer-matrix test but doesn't compose cleanly with the others is automatically dropped during closure. The flip side, also worth understanding: if multiple spurious operations *do* close into a higher-symmetry group, Spacey reports that group. That's the over-promotion failure mode — see [Over-promotion](over-promotion.md).

## Step 5: Optional stability verification

`verify_stable=true` re-runs the full pipeline at a tolerance 1000× tighter and warns if the resulting group order is different. The idea: a *real* symmetry should still be a symmetry at three orders of magnitude tighter tolerance; if it isn't, the answer was tolerance-dependent. This catches near-boundary cases (a tetragonal lattice that *almost* looks cubic) cheaply — at the cost of one extra invocation per call. See [Detect tolerance-dependent answers](../how-to/detect-tolerance-dependence.md) for the user-facing recipe.

## What changes for `spacegroup`

`spacegroup(c::Crystal)` does steps 1–5 on `c.A` to get the lattice point group, then for each candidate rotation `R` it enumerates candidate translations `τ` by looking at differences of probe-atom positions. Each `(R, τ)` is verified by `isSpacegroupOp` (does the op map each atom to an atom of the same type, modulo the lattice?) and the surviving operations are returned. The `verify_stable` flag works analogously, re-running at `pos_tol/1000`.

## How Spacey compares to other tools

| Code | Approach | Strengths | Trade-offs |
|---|---|---|---|
| **Spacey.jl** | Mink-reduce + 27-candidate enumeration + group closure | Provably complete; transparent tolerance; no naming/Niggli machinery to disagree with | No standard-setting names |
| **spglib** | Iterative tolerance + atom-distance matching | Universally adopted; provides space-group names | Tolerance is absolute (Å), heuristic adjustments; no over-promotion warning |
| **AFLOW-SYM** | Iterative tolerance + Niggli + table lookup | Full Bravais classification | Niggli boundary cases; coupled to large AFLOW corpus |
| **FINDSYM** (ISOTROPY) | Tolerance-driven testing of all 230 space groups | Rich post-processing (group–subgroup, distortions) | Slower for batch processing |

Spacey's deliberate scope is *just* the symmetry operations themselves, with explicit tolerance and stability checking. For names, Niggli reduction, or standard settings, use one of the other tools downstream.

## See also

- Explanation: [Why Minkowski reduction](why-minkowski.md), [Tolerances](tolerances.md), [Over-promotion](over-promotion.md), [Validation strategy](validation-strategy.md)
- Reference: [`pointGroup`](../reference/point-groups.md), [`spacegroup`](../reference/space-groups.md)
- How-to: [Find a point group](../how-to/find-pointgroup.md), [Find a space group](../how-to/find-spacegroup.md)
