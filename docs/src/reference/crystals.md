# Crystals

The `Crystal` type bundles a lattice with atomic positions and type labels. All space-group routines accept a `Crystal`. Positions are stored internally in fractional coordinates regardless of the form used at construction.

## Construction and accessors

```@docs
Crystal
fractional
cartesian
read_poscar
```

## Tolerances and classification

```@docs
default_pos_tol
crystal_system
```

## Primitivization

```@docs
is_primitive
make_primitive
```

## Index

```@index
Pages = ["crystals.md"]
```

See also:

- How-to: [Construct a crystal](../how-to/construct-a-crystal.md), [Classify the Bravais system](../how-to/classify-bravais.md), [Handle noisy real-world data](../how-to/handle-noisy-data.md)
- Explanation: [Tolerances](../explanation/tolerances.md) (the `default_pos_tol` formula), [Crystal system vs full Bravais](../explanation/crystal-system-vs-bravais.md)
