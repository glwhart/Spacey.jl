# symlib research notes for Spacey.jl

A code review of [symlib](https://github.com/msg-byu/symlib) — the Fortran symmetry-finding library that predates Spacey.jl — capturing what could usefully migrate, what shouldn't, what algorithmic differences to learn from, and what advantages Spacey already has that the symlib README could point at.

The recommendations are prioritized by impact on Spacey's primary downstream consumer, [Enumlib.jl](https://github.com/glwhart/Enumlib.jl), which currently uses Spacey only for `pointGroup(A)` (3D lattice point group from a basis matrix) and gets its HNF/SNF infrastructure from [NormalForms.jl](https://github.com/glwhart/NormalForms.jl).

---

## What's in symlib

`src/` contains roughly 4,600 lines of Fortran across these modules:

| Module | Purpose |
|---|---|
| `symmetry.f90` (875 lines) | Point-group + space-group finders, `make_primitive`, `bring_into_cell`, `find_site_equivalencies` |
| `vector_matrix_utilities.f90` (578) | Minkowski reduction (`minkowski_reduce_basis`, `reduce_C_in_ABC`, `gaussian_reduce_two_vectors`, `minkowski_conditions_check`), determinant, cross product, orthogonality defect |
| `rational_mathematics.f90` (631) | `HermiteNormalForm`, `SmithNormalForm` (and a `_li` long-int variant), GCD utilities |
| `compare_structures.f90` | `compare_arbitrary_structures`, `is_equiv_lattice`, `is_derivative`, `is_lattice_point`, `mapping_translation_exists` |
| `numerical_utilities.f90` (343) | `equal` (with relative + absolute tol pair) for ranks 0–3 |
| `group_theory.f90` (~70) | `grouper` — generate a permutation group from a list of generators |
| `combinatorics.f90`, `itertools.f90`, `classes.f90` | Imported from the polya project; combinatorial enumeration helpers |
| `utilities.f90`, `num_types.f90` | Reallocation helpers + double-precision kind |

Drivers: `symdriver.f90` (`symcheck.str`-format input), `symdriverPOSCAR.f90`, `spacegroup.f90`, `inverse_check.f90`, `symcheck.jl` (a Julia harness that drives the Fortran `sym.x` binary on MTP-style config files).

Test infrastructure: every public routine has a `*.xml` fortpy spec listing input cases (often 1–1000 per routine) and reference output files. The actual test data (input fixtures + reference outputs) lives in `../tests/<module_name>/` directories; the fortpy framework runs the routine on each case and diffs against the reference.

---

## Where Spacey already wins (for the symlib README pointer)

These are the items worth listing in symlib's README to redirect users who don't have a hard Fortran requirement:

1. **Native Julia, no FFI.** Spacey ships as a standard `Pkg.add` install; symlib needs gfortran or ifort and a Makefile. For mixed-language projects that can absorb Fortran this is fine, but Julia-only consumers (Enumlib.jl, autoGR, downstream ML pipelines) take a meaningful hit.
2. **Two-tolerance design with calibrated defaults.** Symlib uses a single `eps` (default `1e-10`) plus an internal absolute tolerance (`atol = 5e-4`, hardcoded). Spacey separates lattice tolerance (`tol`, relative, volume-normalized) from position tolerance (`pos_tol`, absolute, formula-driven `0.01 · (V/N)^(1/3)`). The two scales correspond to physically distinct noise sources — see `docs/src/explanation/tolerances.md` and `designDiscussions.md`.
3. **`verify_stable` opt-in over-promotion detection.** Spacey re-runs at 1/1000 of the user's tolerance and warns if the answer changes. Symlib has no equivalent; an over-promoted answer is silently returned.
4. **Provably-complete symmetry search via Minkowski reduction.** Spacey's 27-candidate search on a Mink-reduced basis is theoretically complete with a small fixed search space. Symlib's `get_lattice_pointGroup` builds a sphere-of-lattice-points search up to the max basis vector length (`n1, n2, n3` ceilings computed from `max_norm · |b_j × b_k|/V`), which is correct but not size-bounded ahead of time and doesn't depend on prior reduction. (See "Algorithmic comparisons" below.)
5. **3,400-test suite with the full 1095-prototype AFLOW corpus.** Symlib's README claims it is "~10% unit tested" with fortpy. Spacey runs the AFLOW Library of Crystallographic Prototypes (Mehl 2017, Hicks 2019, Hicks 2021) on every push, against two independent invariants (op count + crystal system) and uses `@test_broken` as an explicit deviation catalog.
6. **`SpacegroupOp` type with operator overloads.** Composition (`*`), inverse (`inv`), callable (`op(r)`), exact equality with mod-1 canonicalization. Symlib returns separate `rots(:,:,:)` and `shifts(:,:)` arrays; downstream code re-implements composition by hand.
7. **τ canonicalization via snap-to-rational** (denominator ≤ 12). Eliminates Float64 drift in space-group composition tests; symlib doesn't canonicalize τ.
8. **Diátaxis documentation tree.** Tutorials, how-to guides, reference, explanation; published at `https://glwhart.github.io/Spacey.jl/stable`. Symlib has function-header comments and a brief README.
9. **Active maintenance.** Spacey targets Julia 1.11; symlib's last revision was 2.0.4 with a comment that "many small changes ... were made so long ago I don't remember why".

A short paragraph in symlib's README pointing at Spacey for new Julia users — keeping symlib for the existing Fortran-dependent consumers (UNCLE, parts of enum4) — would let users self-select. Wording suggestion saved at the end of this file.

---

## Algorithmic comparisons

### Point-group finding

| | symlib `get_lattice_pointGroup` | Spacey `pointGroup` |
|---|---|---|
| Pre-condition | Accepts arbitrary basis | Requires Minkowski-reduced (vector form errors out otherwise; matrix wrapper auto-reduces) |
| Search space | Sphere of lattice points up to `max_norm`; counted via `n_i = ⌈max_norm · |b_j × b_k|/V⌉` | Fixed 27-vector neighborhood `{-1, 0, 1}³` |
| Filter pipeline | Length match → unique-vector check → volume → orthogonality (`R Rᵀ = I`) | Length match → volume → integer-matrix test (`inv(A) · B ≈ Z`) → group closure |
| Tolerance | Single `eps` (default 1e-10), augmented internal `atol = 5e-4` | Volume-normalized `tol` (default 0.01) for the integer test; volume normalization makes `tol` unit-independent |
| Output form | `lattpg_op(:,:,:)` — Cartesian rotation matrices (real) | Tuple `(LG, G)` with both lattice-coordinate integer matrices and Cartesian rotations (the `G` half is on the v1.0 chopping block — see "v1.0 task list" below) |
| Theoretical foundation | Implicit (sphere bound is a sufficient condition; not minimal) | Explicit (27-neighbor theorem on Mink-reduced basis is provably complete and minimal) |

**The 27-neighbor approach is the cleaner story** — provably complete, smaller candidate set, the standard textbook foundation. Symlib's sphere-search is correct (the basis vectors of any symmetry image have length `≤ max_norm`) but searches more candidates than necessary. For typical inputs this is invisible overhead; for high-aspect-ratio bases the sphere search can balloon.

### Space-group finding

Both codes follow the same recipe:

1. Find the lattice point group (above).
2. For each candidate rotation `R`, enumerate candidate translations `τ` by taking differences of probe-atom positions of the same type.
3. Verify each `(R, τ)` by checking it maps the full atomic set onto itself modulo the lattice.

The differences are surface-level:

- Symlib stores `(rots, shifts)` as separate arrays; Spacey returns a `Vector{SpacegroupOp}` with the rotation/translation paired.
- Symlib's `get_spaceGroup` accepts non-primitive cells if `make_primitive` was called first; Spacey requires the user to supply the cell they want analyzed (no auto-primitivization).
- Symlib uses the first atom as the probe; Spacey uses the atom type with the smallest count (smaller candidate-`τ` set per `R`, fewer `isSpacegroupOp` calls).
- Symlib has no equivalent of Spacey's identity-at-index-1 guarantee or `verify_stable` flag.
- Spacey canonicalizes `τ` to `[0, 1)` via snap-to-rational (denominators ≤ 12) at `SpacegroupOp` construction. Symlib does not, leading to closure-test issues on trigonal-screw structures unless tolerances are tuned.

### Minkowski reduction

Symlib has its own (internal) Minkowski reduction in `vector_matrix_utilities.f90` (`minkowski_reduce_basis`, `reduce_C_in_ABC`, `gaussian_reduce_two_vectors`). Spacey delegates to MinkowskiReduction.jl, which is a separate package shared with autoGR. The algorithms are equivalent (same loop structure, same termination condition); the choice is just whether the reduction logic lives inside the symmetry library or as a separate dependency.

### Tolerance handling

Symlib uses one tolerance (with a hardcoded `atol = 5e-4` companion in several routines). The `equal` function in `numerical_utilities.f90` accepts a relative + absolute pair and applies `|a - b| < rtol · max(|a|, |b|) + atol`. This is a sound design but uniformly applied to lattice and atom-position comparisons — symlib doesn't separate the two scales the way Spacey does.

Symlib has no over-promotion detector. The `check_spaceGroup` routine verifies group closure on the *output* of `get_spaceGroup`, which catches algorithmic bugs but not tolerance-driven over-promotion.

---

## Symlib functions worth porting to Spacey (gap analysis)

Prioritized by impact on Enumlib.jl use cases.


### High priority: belongs in Spacey

#### `make_primitive(c::Crystal) → Crystal`

Symlib's `make_primitive` (lines 324–499 of `symmetry.f90`) takes a possibly non-primitive crystal and returns a primitive equivalent: it finds atom-position-preserving fractional translations within the cell, treats those as candidate primitive lattice vectors, and verifies that the resulting basis describes the same lattice in primitive form. The interface returns the reduced atom list, new lattice vectors, and an optional `removed_` index of dropped atoms.

**Why it matters for Enumlib.** Enumlib generates supercells from a primitive cell, so the forward direction is "small → big". But users supply *cells* that may or may not be primitive — and Enumlib currently has no clean way to verify or reduce them. Today, callers must trust their input; a `make_primitive` predicate (and `is_primitive` boolean wrapper) closes that loop. Also useful for general Spacey users importing CIF / POSCAR data, which often comes in conventional rather than primitive form.

**Code-sharing constraint** (per user direction). The "find fractional self-translations of the structure" logic is the *same operation* that Spacey's `spacegroup` already does to enumerate candidate translations `τ` per rotation `R` (with `R = I` it reduces to the make-primitive case). Implementation must factor that common logic into an internal helper (something like `_find_self_translations(c::Crystal; pos_tol)::Vector{Vector{Float64}}`) and have both `make_primitive` and `spacegroup` call it. **No duplication.**

**Implementation cost.** ~50 lines of Julia for the public `make_primitive` + `is_primitive` once the shared helper is factored. The trickier piece is the "find a primitive triplet from {existing basis vectors} ∪ {found internal translations}" search — a triple-nested loop already present in the Fortran. The volume-conservation test (`abs(det(B_new)) ≈ V_orig / k` where `k` is the centering multiplicity) is direct.

**API sketch:**
```julia
make_primitive(c::Crystal; pos_tol = default_pos_tol(c)) → (Crystal, removed_indices)
is_primitive(c::Crystal; pos_tol = default_pos_tol(c)) → Bool
```

#### `is_equiv_lattice(A, B)` and `is_derivative(parent, child)`

Two-line predicates from `compare_structures.f90`, operating on bare basis matrices (not Crystals). They belong in Spacey because the test is purely geometric (integer-matrix relationship between two lattices) — no atoms, no symmetry operations applied.

```julia
is_equiv_lattice(A, B; tol) = let S = inv(A) * B
    isapprox(abs(det(S)), 1.0; atol=tol) && isapprox(S, round.(S); atol=tol)
end

is_derivative(parent, child; tol) = let S = inv(parent) * child
    isapprox(S, round.(S); atol=tol)  # integer relationship, no volume constraint
end
```

**Why it matters for Enumlib.** These are *exactly* the integer-matrix tests that `basesAreEquiv` and `getFixingLatticeOps` already do inline in Enumlib.jl (`src/LatticeColoringEnumeration.jl`). Putting them in Spacey gives Enumlib a single source of truth for "are two lattices the same up to a unimodular transformation" / "is one a sublattice of the other" — currently each downstream consumer reimplements the test with its own `1e-6` constant.

**Implementation cost.** Trivial; ~10 lines each plus tests. Public exports (no `Spacey.` qualification needed).


##### Unit-test plan

Each predicate is small enough that the test sets are the load-bearing part of the design:

```julia
@testset "is_equiv_lattice" begin
    I3 = Matrix{Float64}(I, 3, 3)

    # 1. Identity / reflexive case
    @test is_equiv_lattice(I3, I3)

    # 2. Unimodular transforms (det = ±1, integer entries) preserve the lattice
    for M in [
        [1 1 0; 0 1 0; 0 0 1],            # shear, det = 1
        [-1 0 0; 0 1 0; 0 0 1],           # reflection, det = -1
        [0 1 0; 1 0 0; 0 0 1],            # axis swap, det = -1
        RandUnimodMat3(8),                # random unimodular
    ]
        @test is_equiv_lattice(I3, I3 * M)
        @test is_equiv_lattice(I3 * M, I3)   # symmetric
    end

    # 3. Different volume → NOT equivalent
    @test !is_equiv_lattice(I3, 2 * I3)
    @test !is_equiv_lattice(I3, [2 0 0; 0 1 0; 0 0 1])

    # 4. Same volume but non-integer transform → NOT equivalent
    A_skew = [1.0 0.5 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]   # det = 1, not integer-related to I
    @test !is_equiv_lattice(I3, A_skew)

    # 5. Tolerance respected: noisy near-equivalent
    A_noisy = I3 + 1e-12 * randn(3, 3)
    @test  is_equiv_lattice(I3, A_noisy; tol = 1e-8)   # noise below tol
    @test !is_equiv_lattice(I3, A_noisy; tol = 1e-15)  # noise above tol

    # 6. Self-equivalence under each Bravais point group
    #    (the LG operations applied to the basis must produce equivalent lattices)
    for prototype_basis in BRAVAIS_HAND_PANEL  # the 14-Bravais list already in test/runtests.jl
        LG, _ = pointGroup(prototype_basis)
        for op in LG
            @test is_equiv_lattice(prototype_basis, prototype_basis * op)
        end
    end
end

@testset "is_derivative" begin
    I3 = Matrix{Float64}(I, 3, 3)

    # 1. Equivalent lattices are trivially derivatives (det = 1 case)
    @test is_derivative(I3, I3)
    @test is_derivative(I3, I3 * [1 1 0; 0 1 0; 0 0 1])

    # 2. Integer supercells are derivatives
    for hnf in [[2 0 0; 0 1 0; 0 0 1],     # 2x in one direction
                [2 0 0; 0 2 0; 0 0 2],     # 8x cubic
                [3 1 0; 0 2 0; 0 0 1]]     # general HNF, det = 6
        @test is_derivative(I3, I3 * hnf)
    end

    # 3. Non-integer transform → NOT a derivative
    @test !is_derivative(I3, [1.0 0.5 0; 0 1 0; 0 0 1])

    # 4. The derivative relationship is one-way: parent IS derivative of itself,
    #    but parent is NOT a derivative of a 2× supercell (the inverse transform
    #    has half-integer entries)
    parent = I3
    super = parent * [2 0 0; 0 1 0; 0 0 1]
    @test  is_derivative(parent, super)   # parent's lattice contains super's
    @test !is_derivative(super, parent)   # super's lattice does NOT contain parent's

    # 5. Cross-check against the existing HNF generator (if present)
    #    Every output of getAllHNFs(n) applied to a parent must be a derivative
    for n in (2, 3, 4)
        for hnf in getAllHNFs(n)            # uses Enumlib's existing HNF enumerator
            @test is_derivative(I3, I3 * hnf)
        end
    end
end
```

The Bravais panel cross-check (item 6 of `is_equiv_lattice`) and the HNF cross-check (item 5 of `is_derivative`) are the two most valuable tests — they wire the new predicates into Spacey's existing test invariants and Enumlib's existing primitives so a regression in either camp is immediately visible.


### Medium priority: probably belongs in Enumlib, not Spacey (deferred)

#### `compare_crystals(c1, c2; pos_tol)` → defer to Enumlib

Symlib's `compare_arbitrary_structures` returns one of six status codes:

| Status | Meaning |
|---|---|
| 0 | Structures are equivalent |
| 1 | Inequivalent (but lattices match, atoms are on the same lattice, atom counts match, type counts match) |
| 2 | Underlying lattices are different |
| 3 | Atoms of `str2` do not lie on `str1`'s lattice |
| 4 | Different number of atoms |
| 5 | Different number of types |

The whole idea of "are these two structures equivalent" is central to the enumeration problem (deduplicating generated configurations) but separate from symmetry-finding (which is what Spacey is for). Spacey provides operations; Enumlib decides what to do with them.

**Decision: belongs in Enumlib.jl, not Spacey.** Saved as a project-memory note for the future Enumlib work; symlib's `compare_structures.f90` is the reference implementation when that lands.


### Not needed in Spacey

- **`bring_into_cell`.** The `Crystal` constructor folds positions to `[0, 1)` automatically (the `mod 1` change earlier this session). For users who already have a `Crystal`, this routine has no work to do. For users with raw Cartesian points outside the `Crystal` abstraction, two lines of inline `inv(A) * r` + `mod.(_, 1.0)` is enough — not worth a public utility.
- **`HermiteNormalForm` / `SmithNormalForm`.** Already in [NormalForms.jl](https://github.com/glwhart/NormalForms.jl) and [SmithNormalForm.jl](https://github.com/wildart/SmithNormalForm.jl). No reason to absorb into Spacey.
- **`grouper`.** Generates a permutation group from generators. Useful but Enumlib already implements its own `getTransGroup`/`getPermG` for the specific case it needs.
- **`get_spaceGroup_atomTypes`.** UNCLE-specific label-coupling logic. Not generally useful.
- **`find_site_equivalencies`.** Cluster-expansion-specific (only used in UNCLE per symlib's HISTORY). May belong in [JuCE.jl](https://github.com/glwhart/JuCE.jl) eventually — saved as a project-memory note.
- **`check_spaceGroup`.** Validates that an op set is a group. Spacey's own `isagroup` already covers the integer-matrix case; the floating-point-tolerance variant is a small generalization but not a high-impact gap.
- **`put_pointGroup_in_latticeCoords`.** Spacey already returns both forms (`(LG, G)`) — and once `G` is dropped (see v1.0 task list), the lattice-coord form will be the default and the Cartesian form will be a separate helper.
- **The fortpy XML test infrastructure.** Spacey's `runtests.jl` + Documenter doctest model is structurally cleaner; the symlib XML specs are tightly coupled to fortpy and don't transfer well.

---


## Test corpus

Symlib's tests are fortpy-driven: each routine has an XML spec declaring up to 1000 input cases, and the framework runs the routine on each input file and diffs against a reference output file. The test data lives outside the source tree (`../tests/<module>/<routine>_*.in.<N>`) and isn't in the local clone, but the *patterns* are visible in the XML.

### What to actually adopt

A close audit of MinkowskiReduction.jl's existing test suite shows it **already has extensive random-basis testing** — about 290 unimodular-scrambled cases across the cubic / FCC / BCC / 7-Bravais panels, plus aspect-ratio + noise sweeps up to AR ~ 10⁸, plus permutation/sign invariance, idempotence, and transform-matrix sanity sweeps. Adding "1000 unimodular-scrambled cases" would be largely redundant.

**The marginal addition that's not redundant: pure random Gaussian inputs.** All current MinkowskiReduction.jl random tests start from a structured Bravais lattice and apply a `RandUnimodMat3` (`|det| = 1`) transform. Truly unstructured random Float64 inputs — `randn(3, 3)` with arbitrary determinant and aspect ratio — exercise a regime not currently covered. **Action:** add a small (~100-case) pure-random panel in MinkowskiReduction.jl/test/runtests.jl. Lives in MinkowskiReduction.jl, not Spacey.

### What's *missing* from symlib's tests that Spacey has

Three test patterns Spacey uses that don't have a direct symlib analog:

#### 1. Near-boundary heatmap diagnostics

`test/nearMissBoundary.jl` (point group) and `test/nearMissBoundaryCrystal.jl` (BaTiO₃-style space group) sweep a `(ε, tol)` 2D grid and record the operation count Spacey returns at each cell. ε is the structural distortion (how far the input is from a higher-symmetry parent — e.g. tetragonal `c = 1 + ε` for `pointGroup`, Ti displacement ε for `spacegroup`); tol is Spacey's tolerance parameter. The output is an ASCII heatmap with cells marked according to whether Spacey reports the right (lower) symmetry, the over-promoted (higher) symmetry, and whether `verify_stable` would have caught the over-promotion. The diagonal `tol ≈ ε` failure region is visually striking.

These diagnostics are *not part of CI* (they take seconds to minutes to run as full sweeps), but **specific cells are pinned by unit tests** — the over-promoted-and-warning-fires case, the stable-and-no-warning case, etc. The diagnostic file itself doubles as a regression checkpoint: re-running it lets a developer see at a glance whether the failure region has shifted.

Symlib doesn't have anything equivalent — its tests are individual input-output regression checks, not 2D parameter sweeps.

#### 2. AFLOW corpus dual-invariant validation

Each AFLOW prototype label encodes **two independent invariants**:

- **Operation count** — derivable from the space-group number (third underscore-separated segment of the prototype label, `1`–`230`); the order is a known integer per the ITA tables.
- **Crystal system** — derivable from the Pearson-symbol prefix (first one or two letters: `a` triclinic, `m` monoclinic, `o` orthorhombic, `t` tetragonal, `c` cubic, `hP` hexagonal, `hR` trigonal).

Spacey's test suite checks **both** invariants for every prototype. The same AFLOW structure is fed to `spacegroup(c)` for op-count comparison and to `crystal_system(c)` for system comparison — these are independent code paths in Spacey, and different bugs surface in each. Some Spacey deviations pass op-count but fail crystal-system (lattice over-promotion that happens to round to a "right-looking" op count); others pass crystal-system but fail op-count (atomic-position mis-identification on a correctly-labeled lattice). The two invariants together catch a strictly larger deviation set than either alone.

This is the conceptual win behind shipping `crystal_system` cheaply (it's ~10 lines and re-uses `pointGroup`). Without that second invariant we'd have only one test signal per prototype and would miss a real class of deviations.

Documented in `docs/src/explanation/validation-strategy.md`.

#### 3. Cross-validation between `_simple` / `_fast` / `_robust`

Spacey ships three internal point-group finders (none exported directly; reachable as `Spacey.<name>`):

- `pointGroup_simple` — brute-force iteration over all 19,683 candidate matrices in `{-1,0,1}⁹`, filtered by `T = UᵀU ≈ I` with default `isapprox`. Slow but transparent.
- `pointGroup_fast` — the same pipeline with norm/volume pre-filters. Fast, strict tolerance, no `tol` knob.
- `pointGroup_robust` — same pipeline plus `tol`-controlled integer-matrix test and group-closure post-step. The version `pointGroup` actually delegates to.

The test suite **cross-checks all three on every clean (non-noisy) Bravais lattice** — they must produce the same operation count for each of the 14 hand-built prototypes. If any of the three drifts away from the others, that's an immediate signal: either a bug in the variant, or a regression in shared code.

This is `pointGroup_simple`'s only purpose. It exists to be a reference truth — slow but correctness-by-inspection — that the optimized paths must match. Symlib has only one point-group implementation and consequently no cross-validation of this form.

---


## v1.0 task list (in priority order)

**Status: all six items shipped in v0.8.0 (2026-04-29).** Items 1–5 landed in
Spacey.jl; item 6 (MinkowskiReduction.jl Gaussian stress panel) landed in that
package's repo earlier in the same review pass. The v0.8.0 → v1.0.0 bump is
held back deliberately to season the new APIs (predicates, `make_primitive`,
`read_poscar`) under real downstream use before committing to long-term API
stability. Reasons captured in the Spacey.jl release notes for v0.8.0.

User-confirmed scope, with effort estimates:

1. **Allow non-Mink-reduced input to the three-vector form of `pointGroup`** (~30 min). The matrix-form wrapper already auto-reduces; the three-vector form errors out. Symmetric handling removes the footgun. Add a kwarg `auto_reduce=true` (default) with the current behavior preserved as `auto_reduce=false` for users who want to assert their input is already reduced. **— Done, commit 95b9959.**

2. **`is_equiv_lattice` / `is_derivative` predicates in Spacey** (~1 hour). Lattice-only predicates (no Crystal, no atoms). Test plan above. **— Done, commit 74a5c0e. `is_derivative` docstring strengthened in cf52f62 to lead with the directional framing (the asymmetry is load-bearing, not a footgun).**

3. **`make_primitive` / `is_primitive`** (~½ day, including refactoring `spacegroup` to share the "find fractional self-translations" helper). Closes the most visible functionality gap relative to symlib. Hard requirement: no duplicated code with `spacegroup`'s τ-enumeration step. **— Done, commit 546cb70. Three internal helpers (`Spacey._probe_atoms`, `Spacey._find_translations_for_rotation`, `Spacey._find_self_translations`) factored from `spacegroup`'s τ-loop; both `spacegroup` and `make_primitive` go through them.**

4. **POSCAR reader** (~½ day). Function that reads a POSCAR-format file and returns a `Crystal`. Lots of users would benefit (the symlib drivers all use POSCAR-like input, and downstream ML-training-data workflows consume POSCAR by default). Could live in Spacey or as a small companion package; user preference is to ship in Spacey for accessibility. **— Done, commit fb20213. Supports VASP 4 / VASP 5+, Direct/Cartesian coords, scaling factor (positive or negative-as-target-volume), `Selective dynamics` line.**

5. **Drop `G` from `pointGroup`'s return** (substantial breaking change, ~½–1 day for the API + docs sweep + downstream Enumlib update). The `(LG, G)` tuple has been a footgun — users see two outputs, are unsure which to use, and downstream code shows a mix of `LG, _ = pointGroup(A)` (correct) and `LG, G = pointGroup(A)` with `G` then discarded. New API: `pointGroup(A)` returns just `LG`; a separate helper `to_cartesian(LG, A) → G` for users who genuinely want Cartesian rotations. Migrating Enumlib.jl is part of the same change. **— Done first (so items 1–4 could land on the cleaner API). Spacey commit 4f2f029, Enumlib commit daa99af.**

This was previously raised as a v0.8 cleanup in `phase2_plan.md` §6.2; promoted to a v1.0 task because the LG/G confusion compounds the longer it lives.

6. **Pure-random Gaussian stress panel in MinkowskiReduction.jl** (~30 min — see "Test corpus" above). Not Spacey-side but related; the symlib review surfaced it. **— Done earlier in the same review pass; in MinkowskiReduction.jl repo.**

### Explicitly out of scope for v1.0

- `compare_crystals` (deferred to Enumlib; saved to project memory)
- `find_site_equivalencies` (may belong in JuCE.jl; saved to project memory)
- `bring_into_cell` (redundant given the `Crystal` constructor's mod-1 fold)
- `tools/symcheck.jl`-style ML-data validator (user reports no current need)
- `regime` convenience kwarg (deferred; was already on the `designDiscussions.md` "future polish" list)
- `verify_stable` regression on additional symlib-HISTORY edge cases (low priority)

---


## Symlib README addition (proposed wording)

A paragraph that could be added near the top of the symlib README to point Julia and Python users at Spacey:

> **For Julia users (and Python users via the [`juliacall`](https://juliapy.github.io/PythonCall.jl/stable/juliacall/) bridge):** A successor library, [Spacey.jl](https://github.com/glwhart/Spacey.jl), provides the same point-group and space-group functionality natively in Julia. Spacey ships as a standard `Pkg.add` install (no Fortran compiler required), uses a provably-complete 27-neighbor symmetry search on a Minkowski-reduced basis (smaller candidate set than symlib's sphere search), separates lattice and atomic-position tolerances explicitly, and includes an opt-in `verify_stable` flag that re-runs at tighter tolerance and warns on disagreement (catching silent over-promotion that symlib does not). Documentation: <https://glwhart.github.io/Spacey.jl/stable>. The two codes will continue to coexist — symlib remains the right choice for Fortran consumers (UNCLE, the legacy enum4 path) — but new Julia projects, including [Enumlib.jl](https://github.com/glwhart/Enumlib.jl) and [autoGR](https://github.com/msg-byu/autoGR), should prefer Spacey.

The "easily called from Python" claim was softened to "via the `juliacall` bridge" — accurate without overstating compared to a direct-pip-install C library like spglib. Tweak to taste; the key claims are fact-checked against `git log` and `docs/src/explanation/algorithm-overview.md`.
