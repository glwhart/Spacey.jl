# Point groups

The public point-group API is one function: [`pointGroup`](#Spacey.pointGroup). It accepts either three basis vectors or a 3×3 matrix and returns the tuple `(LG, G)` of integer-matrix lattice operations and Cartesian rotations.

## Public API

```@docs
pointGroup
```

## Internals

These three variants are kept inside the package and reachable as `Spacey.<name>`. They are not exported and not part of the API contract — signatures may change without a major-version bump. `pointGroup` is a thin wrapper around `Spacey.pointGroup_robust` with the same defaults.

```@docs
Spacey.pointGroup_robust
Spacey.pointGroup_fast
Spacey.pointGroup_simple
```

## Index

```@index
Pages = ["point-groups.md"]
```

See also: [Find a point group](../how-to/find-pointgroup.md), [Why Minkowski reduction](../explanation/why-minkowski.md), [Over-promotion](../explanation/over-promotion.md).
