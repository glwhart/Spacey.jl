# Changelog

## v0.9.0 — 2026-05-15

**Breaking:** package-wide rename to snake_case for consistency with the
Julia style guide. No deprecation aliases — code that worked on v0.8.0
must be updated to use the new names below.

### Renames (public API)

| Old (v0.8)             | New (v0.9)               |
|------------------------|--------------------------|
| `pointGroup`           | `pointgroup`             |
| `snapToSymmetry_SVD`   | `snap_to_symmetry_svd`   |
| `isagroup`             | `is_group`               |
| `isSpacegroupOp`       | `is_spacegroup_op`       |
| `toCartesian`          | `to_cartesian`           |

Kept (already conformant): `Crystal`, `SpacegroupOp`, `spacegroup`,
`fractional`, `cartesian`, `default_pos_tol`, `crystal_system`,
`is_equiv_lattice`, `is_derivative`, `is_primitive`, `make_primitive`,
`read_poscar`.

### Renames (internal, reachable via `Spacey.<name>`)

| Old (v0.8)             | New (v0.9)               |
|------------------------|--------------------------|
| `Spacey.pointGroup_robust` | `Spacey.pointgroup_robust` |
| `Spacey.pointGroup_fast`   | `Spacey.pointgroup_fast`   |
| `Spacey.pointGroup_simple` | `Spacey.pointgroup_simple` |
| `Spacey.snapToSymmetry_avg`| `Spacey.snap_to_symmetry_avg` |
| `Spacey.avgVecOverOps`     | `Spacey.avg_vec_over_ops` |
| `Spacey.aspectRatio`       | `Spacey.aspect_ratio`     |
| `Spacey.threeDrotation`    | `Spacey.rotate_basis_3d`  |

### Other changes

- **Typed exceptions:** argument-validation errors now throw
  `ArgumentError` instead of `ErrorException`. Affects `Crystal`,
  `is_spacegroup_op`, `read_poscar`, `inv(::SpacegroupOp)`, and
  `pointgroup_robust`'s `auto_reduce=false` path. Internal-invariant
  failures still throw `ErrorException` (deliberate distinction).
- **Keyword-only debug flag:** `pointgroup_simple`'s `debug` argument is
  now a keyword (`pointgroup_simple(a1, a2, a3; debug=true)`) rather than
  positional.
- **Style sweep:** indentation, operator/comma spacing, and iteration
  variable naming brought into line with the Julia style guide. Behavior
  unchanged.
- **Examples relocated:** the previously orphaned `src/debugging.jl`,
  `src/cubic_example.jl`, `src/snapExample.jl`, and the three
  `2D_snap_example_*.jl` files at the repo root have been moved into
  `examples/`. They are *not* `include`d by the module — they remain
  scratch / illustrative material, unchanged from v0.8.0.

### How to fix your code

A find-and-replace using the tables above is sufficient for almost all
callers — the rename is mechanical and the signatures of every renamed
function are identical to v0.8.0.
