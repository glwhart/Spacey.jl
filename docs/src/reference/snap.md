# Snapping to symmetry

Take a noisy lattice (and the symmetry operations found by [`pointGroup_robust`](@ref)) and produce a snapped basis whose symmetry is exact.

```@docs
snapToSymmetry_SVD
```

`snapToSymmetry_SVD` is the recommended path. A faster but less-robust alternative, `Spacey.snapToSymmetry_avg`, is documented under [Helpers → Internals](helpers.md#Internals).

## Index

```@index
Pages = ["snap.md"]
```

See also: [Snap a noisy lattice to symmetry](../how-to/snap-to-symmetry.md).
