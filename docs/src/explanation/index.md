# Explanation

Background and rationale: *why* Spacey decides what it decides. Read these when you want to understand the algorithm, the trade-offs, and the failure modes — not just call the API.

- [Why Minkowski reduction](why-minkowski.md) — the structural theorem that makes Spacey's symmetry search provably finite. The 27-neighbor result, with the bounded-neighborhood figure from Hart-Jorgensen-Morgan-Forcade (2019).
- [Algorithm overview](algorithm-overview.md) — the end-to-end pipeline of `pointGroup` and `spacegroup`: Minkowski reduce → 27 candidates → 3-stage filter → group closure → optional `verify_stable`. Includes a comparison table against spglib, AFLOW-SYM, and FINDSYM.
- [Tolerances](tolerances.md) — why two tolerances (`tol`, `pos_tol`), how the defaults were calibrated, and how Spacey's two-tolerance design contrasts with the single-tolerance designs of other tools.
- [Over-promotion](over-promotion.md) — the canonical Spacey failure mode (silent reporting of higher symmetry than the structure actually has) and how `verify_stable` detects it. The `(tol, ε)` boundary regime, walked through with the tetragonal-near-cubic case.
- [Canonicalizing τ](canonicalizing-tau.md) — why `SpacegroupOp` translations snap to small rationals at construction. The round-vs-snap episode (`1/3 + 1/3 ≠ 2/3` in 10-digit decimals, silently breaking trigonal-screw closure tests).
- [Crystal system vs full Bravais](crystal-system-vs-bravais.md) — why Spacey ships the 7-class `crystal_system` function but defers the 14-class Bravais classification.
- [Validation strategy](validation-strategy.md) — how the test suite cross-validates three independent algorithm variants against the AFLOW corpus on two complementary invariants (op count + crystal system), and what the `@test_broken` deviation catalog represents.

If you want to *do* something specific, the [How-to guides](../how-to/index.md) are a better fit. For just the API surface, see the [Reference](../reference/index.md).
