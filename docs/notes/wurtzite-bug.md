# Wurtzite over-symmetric `spacegroup` when called with uniform types

Surfaced during Enumlib chunk 6.5b (Regime-C) cross-validation against the Fortran enumlib corpus on 2026-05-19.

## TL;DR

For shifted wurtzite (hexagonal lattice, 4-atom basis: 2 cation + 2 anion), `spacegroup(c::Crystal)` returns **|G| = 24** when `c.types == [1, 1, 1, 1]` (uniform) and **|G| = 12** when `c.types == [1, 1, 2, 2]` (cation vs. anion distinguished). The 24-op group is mathematically correct for the position set with no type information, but Enumlib's downstream multilattice machinery (R50.2a `dset_perms`) wants the labeled symmetry (12 ops for P6₃mc), and currently calls into `spacegroup` with `ones(Int, ndset)` because the active-vs-inactive distinction isn't known at `ParentLattice` construction time.

The over-symmetric group bleeds into Enumlib's `_filter_perm_group_by_mask` step: that filter drops parent ops mapping different-`allowed_labels` sublattices to each other, but the `getPermG` `unique!()` step has already collapsed multiple distinct ops to a single site-permutation by then, so legitimate ops can't be told apart from illegitimate ones. Result: wurtzite enumeration over-deduplicates (julia=3 vs. Fortran=4 at n=1, julia=10 vs. Fortran=42 at n=2, julia=58 vs. Fortran=260 at n=3).

## Minimal reproducer

```julia
using Spacey
A = [1.0  -0.5                 0.0;
     0.0   0.866025403784439   0.0;
     0.0   0.0                 1.6329931618554521]   # column-vector hexagonal basis
r = hcat([0.0, 0.0, 0.0],
         [1.0/3.0, 2.0/3.0, 0.5],
         [0.0, 0.0, 0.375],
         [1.0/3.0, 2.0/3.0, 0.875])                   # shifted wurtzite dset

uniform = Crystal(A, r, [1, 1, 1, 1]; coords = :fractional)
labeled = Crystal(A, r, [1, 1, 2, 2]; coords = :fractional)

length(spacegroup(uniform))   # → 24
length(spacegroup(labeled))   # → 12  (matches P6₃mc)
```

Standard wurtzite (cation at Wyckoff 2b: `(1/3, 2/3, 0)` and `(2/3, 1/3, 1/2)`) shows the same pattern.

## Why the 24-op answer is mathematically correct (and the question for Spacey)

For uniform types, `spacegroup` finds every isometry of the lattice that maps the position set to itself. With 4 indistinguishable atoms in the shifted-wurtzite dset, the position set has 24 such isometries — twice the P6₃mc op count. The extra 12 are the ops that swap cation positions with anion positions; with uniform types, those are valid (positions are interchangeable).

This is *correct behavior for the input* — `spacegroup` does not infer chemical species. The question is whether `spacegroup` (or a higher-level helper) should expose a way to compute the **multilattice symmetry group induced by an equivalence relation on dset positions** — i.e., "treat positions {1, 2} as one class and {3, 4} as another, give me the ops that preserve this partition."

## Where Enumlib trips

`Enumlib/src/types/parent_lattice.jl:79-84` constructs the parent's space group via:

```julia
crystal = Spacey.Crystal(Af, r, ones(Int, length(ds)); coords = :fractional)
ops = Spacey.spacegroup(crystal)
```

The `ones(Int, length(ds))` was a deliberate choice: at `ParentLattice` construction the `Sites` object (which carries `allowed_labels`) hasn't been seen yet, and the comment notes this gives the "multilattice space group" rather than the Bravais point group. For Regime A (single dset position) and Regime B (uniform `allowed_labels`), it's the right answer. For Regime C (heterogeneous `allowed_labels`), the over-symmetric group's extra ops survive `_filter_perm_group_by_mask` for some HNFs because `getPermG`'s `unique!()` has already coalesced them with legitimate ops on the supercell site set.

## Fix candidates

In rough order of blast radius:

1. **Enumlib-side workaround (smallest).** Pass `Sites` into `ParentLattice` (or compute `space_group` lazily at `enumerate(...)` time, with `Sites` in scope). Use type IDs that encode the position's `allowed_labels` equivalence class: e.g. positions with the same `allowed_labels` get the same type, distinct `allowed_labels` get distinct types. Spacey then returns the right (label-aware) group.

2. **Spacey-side helper.** A new entry point like `spacegroup(c::Crystal, equivalence::Vector{Int})` where `equivalence` is a per-atom class label distinct from `types` (so the caller doesn't have to abuse `types`). Conceptually equivalent to (1) but with the partition lifted into Spacey's API.

3. **Enumlib-side: filter at the parent-space-group level rather than per-supercell.** Drop ops from `parent.space_group` whose `dset_perm` maps a position to a different `allowed_labels` class. Then `getPermG` only sees label-preserving ops, no over-symmetry to deal with. Equivalent in effect to (1); ties the parent-level group to a specific `Sites`, which may or may not be desirable.

(2) is the cleanest API-level fix but the biggest change. (1) is probably the easiest Enumlib-side patch — it doesn't need Spacey to change at all, just a `ParentLattice` API that accepts an explicit per-position equivalence class.

## Corpus data (Enumlib repo)

Per-volume cross-validation numbers and the full per-row `(case, volume, hnf_degen, orbit_size, count)` corpus are checked in at `Enumlib.jl/test/data/chunk6.5_fortran_corpus.csv`. Wurtzite rows are the ones that fail; the other four cases (perovskite, zincblende, half-Heusler, full-Heusler with relabeled inactive Y=Z={2}) all pass.
