# Spacey.jl

Spacey.jl finds the **point group** of a Bravais lattice and the **space group** of a crystal, with algorithms designed to be robust to finite-precision floating-point error in the input.

This documentation follows the [Diátaxis](https://diataxis.fr/) framework — four sections, each with a different purpose. Pick the one that matches what you need right now.

| If you're …                              | Go to …                              | Start with                                 |
|---                                       |---                                   |---                                         |
| **New to Spacey** — want to learn        | [Tutorials](tutorials/index.md)      | [Find the point group of a cubic lattice](tutorials/01-first-pointgroup.md) |
| **Solving a specific task** — need steps | [How-to guides](how-to/index.md)     | [Find a space group](how-to/find-spacegroup.md) or [Find a point group](how-to/find-pointgroup.md) |
| **Looking up a function** — need the API | [Reference](reference/index.md)      | [Crystals](reference/crystals.md), [Point groups](reference/point-groups.md), [Space groups](reference/space-groups.md) |
| **Understanding the why** — want depth   | [Explanation](explanation/index.md)  | [Why Minkowski reduction](explanation/why-minkowski.md), [Algorithm overview](explanation/algorithm-overview.md) |

## 30-second smoke test

```julia
julia> using Pkg; Pkg.add(url="https://github.com/glwhart/Spacey.jl")

julia> using Spacey

julia> LG = pointGroup([1.0, 0, 0], [0, 1.0, 0], [0, 0, 1.0]);

julia> length(LG)   # 48 — the cubic point group
48
```

If that runs and prints `48`, the package is loaded and working.

## What's distinctive about Spacey

- **Provably complete symmetry search.** Spacey's algorithm is exhaustive over a finite, small set of integer matrices (27 candidate vectors per basis vector, never more) thanks to the Minkowski-reduction theorem — see [Why Minkowski reduction](explanation/why-minkowski.md). Other tools use either iterative tightening (spglib, AFLOW-SYM) or hard-coded large search windows (VASP); both have failure modes that Spacey's approach avoids.
- **Two-tolerance design.** A relative `tol` for lattice geometry and an absolute `pos_tol` for atomic positions. They scale with different physical noise sources and shouldn't be conflated. See [Tolerances](explanation/tolerances.md).
- **`verify_stable` for tolerance-dependent answers.** Pass `verify_stable=true` to either `pointGroup` or `spacegroup` and Spacey re-runs at 1/1000 the requested tolerance, warning if the answer changes. Catches over-promotion (silently reporting higher symmetry than the structure has) without committing the user to a tighter answer than they asked for. See [Over-promotion](explanation/over-promotion.md).
- **Validated against the full AFLOW corpus.** 1095 prototypes from three published papers, checked on two independent invariants (op count + crystal system). See [Validation strategy](explanation/validation-strategy.md).

## What's *not* in scope

Spacey is a focused tool. It does *not* provide:

- Space-group **names** (`Pm3̄m`, `Fd3̄m`, …) or other Hermann–Mauguin symbol output. Use [spglib](https://spglib.readthedocs.io/), [FINDSYM](https://iso.byu.edu/), or one of the others linked in the [README](https://github.com/glwhart/Spacey.jl#readme) for that.
- Niggli cell reduction or standard-setting transformations.
- Wyckoff-position analysis.
- Magnetic space groups.

Spacey provides the symmetry *operations* themselves, in a form that downstream tools can consume directly.
