# Snapping to symmetry

Take a noisy lattice (and the symmetry operations found by [`pointgroup`](@ref)) and produce a snapped basis whose symmetry is exact.

```@docs
snap_to_symmetry_svd
```

`snap_to_symmetry_svd` is the recommended path. A faster but less-robust alternative, `Spacey.snap_to_symmetry_avg`, is documented under [Helpers → Internals](helpers.md#Internals).

## Index

```@index
Pages = ["snap.md"]
```

See also:

- How-to: [Snap a noisy lattice to symmetry](../how-to/snap-to-symmetry.md)
- Explanation: [Over-promotion](../explanation/over-promotion.md) (the failure mode that `snap_to_symmetry_svd` deliberately commits to, when invoked)
