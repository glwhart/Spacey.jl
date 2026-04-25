# Helpers

## Public

```@docs
isagroup
```

## Internals

These helpers are not exported. They are reachable as `Spacey.<name>` and are documented here for the small audience that wants to drive Spacey's diagnostic plumbing directly (test scaffolding, lattice quality reporting, faster-but-less-robust symmetry snapping). They are not part of the API contract — signatures may change without a major-version bump.

```@docs
Spacey.aspectRatio
Spacey.threeDrotation
Spacey.snapToSymmetry_avg
```

## Index

```@index
Pages = ["helpers.md"]
```
