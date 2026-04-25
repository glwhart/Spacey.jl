# Crystals

The `Crystal` type bundles a lattice with atomic positions and type labels. All space-group routines accept a `Crystal`. Positions are stored internally in fractional coordinates regardless of the form used at construction.

## Construction and accessors

```@docs
Crystal
fractional
cartesian
```

## Tolerances and classification

```@docs
default_pos_tol
crystal_system
```

## Index

```@index
Pages = ["crystals.md"]
```

See also: [Construct a crystal](../how-to/construct-a-crystal.md), [Classify the Bravais system](../how-to/classify-bravais.md).
