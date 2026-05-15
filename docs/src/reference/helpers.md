# Helpers

## Public

```@docs
is_group
is_equiv_lattice
is_derivative
```

## Internals

These helpers are not exported. They are reachable as `Spacey.<name>` and are documented here for the small audience that wants to drive Spacey's diagnostic plumbing directly (test scaffolding, lattice quality reporting, faster-but-less-robust symmetry snapping). They are not part of the API contract — signatures may change without a major-version bump.

```@docs
Spacey.aspect_ratio
Spacey.rotate_basis_3d
Spacey.snap_to_symmetry_avg
```

## Index

```@index
Pages = ["helpers.md"]
```

See also:

- Reference: [`snap_to_symmetry_svd`](snap.md) (the public, recommended snap-to-symmetry path; `Spacey.snap_to_symmetry_avg` is its faster but less-robust alternative)
- How-to: [Snap a noisy lattice to symmetry](../how-to/snap-to-symmetry.md)
- Explanation: [Validation strategy](../explanation/validation-strategy.md) (the `_simple` / `_fast` / `_robust` cross-validation rationale)
