# Crystal system vs full Bravais

Spacey ships a `crystal_system` function that returns one of seven labels — `:triclinic`, `:monoclinic`, `:orthorhombic`, `:trigonal`, `:tetragonal`, `:hexagonal`, `:cubic`. It does **not** ship a function that returns the 14 Bravais lattices (the seven systems further split by primitive / centered / face-centered / body-centered variants). This page explains why we stopped at seven.

## Two layers of lattice classification

There are two natural levels of granularity for "what kind of lattice is this?":

| Layer | Output | Cost | Coverage |
|---|---|---|---|
| 1 | Crystal **system** (7 classes: tri / mono / ortho / tet / trig / hex / cubic) | Trivial — wrapper over `pointGroup` | Catches lattice over-promotion (`a/b ≈ 1`) and lattice mis-identification within a system |
| 2 | Full **Bravais lattice** (14 classes: P / C / I / F variants of each system) | Non-trivial — needs Niggli or metric-tensor classification with edge cases | Catches centering mis-identification on top of Layer 1 |

`crystal_system(A)` is Layer 1.

## Layer 1 is one lookup table

Each crystal system has a unique *holohedry* — the point group of the lattice itself, with no atomic decoration — and each holohedry has a unique order:

| Order | System | Holohedry |
|---|---|---|
| 2 | triclinic | C_i |
| 4 | monoclinic | C_2h |
| 8 | orthorhombic | D_2h |
| 12 | trigonal | D_3d |
| 16 | tetragonal | D_4h |
| 24 | hexagonal | D_6h |
| 48 | cubic | O_h |

Because the orders are all distinct, getting the system from a lattice is two lines of Julia: run `pointGroup` on the Minkowski-reduced basis, look up the order in the table above. That's exactly what `crystal_system(A)` does. Same tolerance behaviour and same failure modes as `pointGroup` — including over-promotion at boundary cases (a tetragonal lattice with `c/a` very close to 1 will return `:cubic` at default `tol`).

## Layer 2 is a different problem

Going from "system" (7 classes) to "Bravais lattice" (14 classes) means distinguishing primitive (P), base-centered (C), body-centered (I), face-centered (F), and rhombohedral (R) variants *within* a system. For primitive cells (which is what most input formats — POSCAR, CIF — supply), the centering shows up as a geometric pattern in the basis vectors:

- BCC primitive has equal-length vectors with arccos(1/3) inter-vector angles.
- FCC primitive has equal-length vectors with 60° angles.
- C-centered orthorhombic primitive has two specific length-pair relationships.
- … and so on for each of the 14 cases.

The decision tree is roughly thirty branches over the metric tensor `G = AᵀA`, with several edge cases at boundaries between Bravais types. Modern implementations (FINDSYM, spglib, AFLOW-SYM) all do this via Niggli reduction followed by a lookup table from the canonical Niggli parameters. Each handles the boundary cases differently, with documented failure modes around the edges.

That's a few days of careful work, useful but not yet motivated for Spacey's current use cases. Spacey's downstream consumers (cluster expansion, k-point grid generation) consume the symmetry operations directly and don't need the Bravais label — they get all the information they need from `pointGroup` and `spacegroup`. The Layer 1 label is sufficient for the validation use case (see below).

## Why ship Layer 1 at all

Layer 1's cost-to-benefit ratio is excellent for one specific reason: it's a *cheap second invariant for testing*.

The full AFLOW prototype-corpus tests cross-check Spacey's output against the space group encoded in the prototype label. The label encodes both the space group *number* (via the third underscore-separated segment) and the Pearson-symbol crystal system (via the first one or two letters of the symbol):

```
a*  →  triclinic    m*  →  monoclinic   o*  →  orthorhombic
t*  →  tetragonal   c*  →  cubic
hP  →  hexagonal    hR  →  trigonal
```

So for each AFLOW prototype, two independent invariants are checkable: the **op count** and the **crystal system**. Different bugs surface differently in these two tests. Some Spacey deviations pass op-count but fail crystal-system (lattice over-promotion that happens to round to a "right-looking" op count). Others pass crystal-system but fail op-count (atomic-position misidentification on a correctly-labelled lattice). Together, the two invariants catch more deviations than either alone — and the AFLOW corpus run uses both.

Without Layer 1 we'd have only one invariant, and miss the second class of deviations. With Layer 2 we'd catch even more, but the marginal gain over Layer 1 is small (the Pearson-letter centering rarely disagrees with the space-group number when the system already does), and the cost is much higher.

## Why Layer 2 is deferred

Three reasons, in order of importance:

1. **No motivating use case in Spacey's downstream consumers.** Cluster expansion, k-point generation, structural enumeration — none of them need the centering label. They consume operations.
2. **Marginal validation benefit over Layer 1.** The AFLOW-corpus deviation set didn't ask for finer granularity than what Layer 1 already provides.
3. **Boundary-case handling is non-trivial.** The Niggli-classification path has documented edge cases (a hexagonal cell with `c/a` matching the FCC primitive ratio, etc.) that other codes still occasionally get wrong. Spacey's policy is "don't ship machinery we can't fully test" — and Layer 2 testing requires at minimum the curation of structures that exercise each Bravais transition cleanly.

If the use case ever appears (a downstream consumer that needs the centering, or a benchmarking exercise that asks for the 14-way classification), Layer 2 can be added without breaking the existing API: a new function `bravais(A)` returning a 14-class label, with `crystal_system(A)` continuing to return the 7-class one. The 7-class API is permanent regardless.

## See also

- Reference: [`crystal_system`](../reference/crystals.md), [`pointGroup`](../reference/point-groups.md)
- How-to: [Classify the Bravais system](../how-to/classify-bravais.md)
- Explanation: [Validation strategy](validation-strategy.md) — the AFLOW corpus tests that motivate Layer 1
