# Point groups

Spacey provides four point-group finders. All four agree on exact (non-noisy) input for every Bravais lattice; the differences are in tolerance handling, speed, and intended use case. See [Find a point group](../how-to/find-pointgroup.md) for guidance on which to call.

## Robust finder (recommended for real input)

```@docs
pointGroup_robust
pointGroup
```

## Fast and simple finders (clean / validation use)

```@docs
pointGroup_fast
pointGroup_simple
```

## Index

```@index
Pages = ["point-groups.md"]
```

See also: [Why Minkowski reduction](../explanation/why-minkowski.md), [Over-promotion](../explanation/over-promotion.md).
