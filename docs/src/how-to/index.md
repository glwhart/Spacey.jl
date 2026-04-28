# How-to guides

Recipe-style guides that assume you already know what you want and just need the steps. Each guide is short, focused, and runnable.

**Setting up your input**

- [Construct a crystal](construct-a-crystal.md) — `Crystal` constructor patterns: required `coords` kwarg, matrix vs. three-vector input, type labels, common errors.

**Finding symmetries**

- [Find a point group](find-pointgroup.md) — `pointGroup`, the four internal variants reachable as `Spacey.<name>`, rotation invariance, tolerance basics.
- [Find a space group](find-spacegroup.md) — build a `Crystal`, call `spacegroup`, read `R` and `τ`, sanity-check by op count.
- [Classify the Bravais system](classify-bravais.md) — `crystal_system` from a matrix or `Crystal`, with the 7-output table and the over-promotion footgun.

**Working with the output**

- [Compose and apply operations](compose-and-apply-ops.md) — `SpacegroupOp` cookbook: apply, compose, invert, equality, `toCartesian`.
- [Snap a noisy lattice to symmetry](snap-to-symmetry.md) — `snapToSymmetry_SVD` recipe with a 1%-noise near-cubic example.

**Handling noise and uncertainty**

- [Handle noisy real-world data](handle-noisy-data.md) — translating experimental noise into appropriate `tol` and `pos_tol` values; the recommended-defaults table by input source.
- [Detect tolerance-dependent answers](detect-tolerance-dependence.md) — `verify_stable=true` on `pointGroup` and `spacegroup`; what to do when the warning fires.

If this is your first time, start with the [Tutorials](../tutorials/index.md). For *why* the API is shaped this way, see the [Explanation](../explanation/index.md) section.
