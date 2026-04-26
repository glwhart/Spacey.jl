# Space groups

The space-group API. `spacegroup(c)` is the entry point; everything else either supports it or operates on the operations it returns.

## Finding the space group

```@docs
spacegroup
isSpacegroupOp
```

## The `SpacegroupOp` type

A single space-group operation `r ↦ R·r + τ`, in lattice (fractional) coordinates. Composition (`*`), inversion (`inv`), application (`op(r)`), and equality are all overloaded — see the type docstring for details. The translation `τ` is canonicalized to `[0, 1)` with snap-to-rational at construction; this is what makes default `==` and `hash` agree with the periodic-boundary semantics a user expects.

```@docs
SpacegroupOp
toCartesian
```

## Index

```@index
Pages = ["space-groups.md"]
```

See also: [Find a space group](../how-to/find-spacegroup.md), [Compose and apply operations](../how-to/compose-and-apply-ops.md), [Canonicalizing τ](../explanation/canonicalizing-tau.md).
