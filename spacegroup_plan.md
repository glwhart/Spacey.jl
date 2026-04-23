# Space Group Plan for Spacey.jl

This document plans the extension of Spacey.jl from lattice point groups to crystal space groups. Scope is deliberately narrow: **find the symmetry operations `(R, τ)` of a crystal**. Classification into the 230 ITA space-group types and standard-setting transforms are **explicitly out of scope** for this plan; they are deferred to a later effort.

This is a living planning document. Decisions and resolutions will be annotated in place as we work through each phase — items are not deleted when resolved.

---

## 1. Scope

### In scope
- A usable `Crystal` type (lattice + atomic basis + atom types).
- Algorithm to find all space-group operations of a crystal: every pair `(R, τ)` where `R` is a point-group operation of the lattice and `τ` is a translation such that `(R, τ)` maps the crystal to itself.
- Tolerance-aware variant for noisy atomic positions.
- Unit tests against well-known crystals with known operation counts.

### Out of scope (explicit)
- Identifying which of the 230 ITA space groups the returned operation set corresponds to.
- Transforming to standard (ITA) setting.
- Space-group symbol generation (Hermann–Mauguin, Schoenflies).
- Structural analysis beyond detecting the operations (Wyckoff positions, site symmetries, etc.).

### Relationship to existing point-group code
- `pointGroup_robust` applied to the lattice vectors gives the candidate `R` pool. Only a subset may survive as space-group operations depending on the atomic basis.
- The tolerance-based robustness discussion in `research.md` §2.1 carries over directly: tolerance choice at boundary cases (e.g. an atom near but not at a special position) produces the same over-/under-report structure we saw for lattices near symmetry boundaries.

---

## 2. Algorithm

### 2.1 Core idea

For a crystal with lattice `A = [a₁ a₂ a₃]`, atomic positions `{rᵢ}`, and atom types `{tᵢ}`:

1. Compute the point group of the lattice using `pointGroup_robust` — call it `{R_k}`.
2. For each `R_k`:
   1. Apply `R_k` to every atom: `r'ᵢ = R_k · rᵢ`.
   2. `(R_k, τ)` is a symmetry iff there exists a single translation `τ` such that for every atom `i`, there is an atom `j` of the same type with `r'ᵢ + τ ≡ r_j (mod lattice)`.
   3. Candidate `τ` values are differences `r_j − r'_{i₀}` for a single probe atom `i₀` and every atom `j` of the same type as `i₀`. Only finitely many to test (≤ N).
   4. For each candidate `τ`, verify the match holds for every atom.
3. Collect all surviving `(R_k, τ)` pairs.

### 2.2 Complexity
- `|{R_k}| ≤ 48`.
- ≤ N candidate τ per R (N = number of atoms of the probe type).
- Each τ tested against N atoms.
- Total: `O(48 · N²)` position comparisons. Small for typical crystals (N ~ 1–100).

### 2.3 Subtleties
- **Periodic equivalence.** All position comparisons are modulo the lattice. Fractional coordinates make this a `mod 1` comparison.
- **Atom types.** The mapping `i → j` must preserve type.
- **Origin choice.** `τ` depends on the coordinate origin; `R` does not. The *set* of `(R, τ)` pairs modulo origin choice is the crystal's space group.
- **Fractional translations (screw axes, glide planes).** These are the non-trivial `τ` values; the algorithm handles them naturally because it enumerates candidate τ from the atoms themselves.
- **Centering translations.** Pure translations (R = I) with non-zero τ signal a non-primitive (conventional) cell. These should appear as separate ops.

### 2.4 Known reference implementations
- **Spglib** (Togo): lattice-guided R search + operation-matching, with tolerance handling.
- **FINDSYM**: Niggli reduction + table lookup + operation matching.
- **AFLOW-SYM**: multi-pass with position averaging and multiple tolerance levels.

Our plan is the simplest form of this family, minus classification.

---

## 3. Phased implementation plan

### Phase 1 — Infrastructure
- Revise `Crystal` struct: concrete types, validated inner constructor, documented layout convention.
- Add `fractional(c::Crystal)` / `cartesian(c::Crystal)` conversion helpers.
- Add `isSpacegroupOp(R, τ, c::Crystal; tol)` — returns `true` if `(R, τ)` maps `c` to itself within tolerance.
- **Exit criterion:** unit tests for `isSpacegroupOp` pass on trivial cases (identity op, pure lattice translations, obviously-not-symmetric perturbations).

### Phase 2 — Core algorithm
- Implement `spacegroup(c::Crystal; tol)` returning a list of `(R, τ)` pairs.
- Internally: call `pointGroup_robust` on the lattice, enumerate candidate τ per R via atom-difference probing, verify each candidate against all atoms with `isSpacegroupOp`.
- Replace the `spacegroup(c)::Bool = true` stub in the same commit.
- **Exit criterion:** simple cubic lattice with one atom at the origin returns 48 ops, all with `τ = 0`.

### Phase 3 — Tests against known crystals
Add `@testset "Space group"` with, at minimum:

| Crystal | Lattice | Basis | Expected ops | Why in the list |
|---|---|---|---|---|
| Simple cubic, 1 atom | sc | `{(0,0,0)}` | 48 | Trivial, matches point group |
| CsCl | sc | Cs `(0,0,0)`, Cl `(½,½,½)` | 48 (Pm3̄m) | Two-atom primitive, τ = 0 only |
| NaCl | fcc | Na `(0,0,0)`, Cl `(½,½,½)` | 192 (Fm3̄m) | Centering translations must appear |
| Diamond | fcc | C `(0,0,0)`, C `(¼,¼,¼)` | 192 (Fd3̄m) | Non-symmorphic: glide planes present |
| HCP | hexagonal | 2 atoms at Wyckoff `2c` | 24 (P6₃/mmc) | Screw axis test |
| Quartz (α) | trigonal | Si ×3, O ×6 | 6 (P3₁21) | Non-symmorphic, low symmetry, helical |

- **Exit criterion:** all return the expected operation count.

### Phase 4 — Robustness
- Position-matching tolerance distinct from lattice tolerance.
- Noise sweep: random position noise, verify operation count is stable.
- Near-boundary crystals: e.g., an atom near but not at a special position; tetragonal ferroelectric with a slightly-displaced cation.
- `verify_stable`-analog for space groups (re-run at tighter `pos_tol`, warn on disagreement).
- **Exit criterion:** new test panel mirroring the near-boundary point-group tests (`test/nearMissBoundary.jl`-style diagnostic).

### Phase 5 *(deferred, out of current scope)*
Classification into 230 ITA types, standard-setting transforms, space-group symbols, Wyckoff position analysis.

---

## 4. Open design choices

These need resolution before Phase 1 begins. Initial recommendations given; final decisions to be annotated in place.

### 4.1 Coordinate convention for atomic positions

Note: this question is *only* about atomic positions `r`. The point-group code is basis-agnostic (operations come back as integer matrices in lattice coordinates; no "fractional vs Cartesian" distinction applies). Lattice vectors `a₁, a₂, a₃` are stored as 3-vectors in whatever basis the user supplied.

**Options:**
- **A. Fractional (lattice coordinates).** `r` stored as `(f₁, f₂, f₃)` with `r = f₁a₁ + f₂a₂ + f₃a₃`. Entries typically in [0, 1). Periodic-equivalence check is `mod 1`. Matches ITA convention. Pairs naturally with integer-matrix R and rational τ.
- **B. Cartesian.** `r` stored in the same basis and units as `a₁, a₂, a₃`. Every periodic-equivalence check needs `inv(A) · r` to reduce mod lattice.
- **C. Accept either at construction, normalize to one canonical form internally.** Constructor takes a `coords=:fractional` / `:cartesian` keyword; stored form is always fractional.
- **D. Store both, lazy conversion.** Sync risk, wasted memory; rejected.

**Recommendation:** **C**, with internal canonical form = fractional. Users coming from CIF get fractional for free; users coming from XYZ-style tools can hand in Cartesian and be converted once at construction. Algorithm code stays simple (fractional-only inner loop). Explicit keyword at the call site prevents silent-wrong-answer errors from passing the wrong form without a flag.

**Decision (2026-04-23):** **C with the `coords` keyword *required* at construction — no default.**

```julia
Crystal(A, r, types; coords=:fractional)   # explicit
Crystal(A, r, types; coords=:cartesian)    # explicit
Crystal(A, r, types)                       # ERROR: coords kwarg required
```

Rationale: the `(0.5, 0.5, 0.5)` body-centre ambiguity is unresolvable without explicit declaration (a range check on fractional input catches Cartesian-with-large-values but not small-lattice Cartesian). Since the failure mode is a silent wrong space group — no crash, no warning — maximising user safety is worth one extra keyword per call. Internal storage is always fractional regardless of how the user supplied it.

### 4.2 Crystal struct shape

Current (in `src/Spacey.jl`):
```julia
struct Crystal
    a1
    a2 
    a3
    r::Array{Float64,2}
    a::Array{Int}
end
```

**Issues:**
- `a1, a2, a3` untyped.
- Whether columns or rows of `r` are atoms is not specified.
- No validation that `size(r, 2) == length(a)`.

**Proposed revision:**
```julia
struct Crystal
    A::Matrix{Float64}      # 3×3 lattice, columns = a1, a2, a3
    r::Matrix{Float64}      # 3×N fractional positions, columns = atoms
    types::Vector{Int}      # length N
    function Crystal(A, r, types)
        size(A) == (3, 3) || error("A must be 3×3")
        size(r, 1) == 3 || error("r must have 3 rows (fractional coords)")
        size(r, 2) == length(types) || error("r columns must match types length")
        new(A, r, types)
    end
end
```

**Sub-question:** `types::Vector{Int}` or `Vector{Symbol}` (`:Fe`, `:O`)? Int is simpler and matches the current stub; Symbol is more expressive and harder to misuse. We could also support both via a type parameter: `struct Crystal{T}; ... types::Vector{T}; end`.

**Recommendation:** Parameterized `Crystal{T}` with `T <: Union{Int, Symbol, String}`. Zero-cost generality.

**Decision (2026-04-23):**

1. **Numeric storage:** `A::Matrix{Float64}`, `r::Matrix{Float64}` (strict). The constructor accepts `AbstractMatrix{<:Real}` and converts once, so users can pass `Int`, `Rational`, `Float32`, etc. without friction — but the library has exactly one numeric form to reason about internally. BigFloat for extreme aspect ratios is left for a future parameterised twin; not unlocked by just relaxing `Crystal`.

2. **Types field:** parameterised `Crystal{T}` with `types::Vector{T}`, no subtype restriction. The algorithm only ever compares types with `==`, which works for any `T`. Users pick whatever they like — `Int` (compact), `Symbol` (`:Fe`, `:O`), `String` (`"Fe"`), or a custom type. Zero cost to the library.

3. **Three-vector convenience constructor** (resolved 2026-04-23):

   ```julia
   Crystal(a1, a2, a3, r, types; coords) = Crystal(hcat(a1, a2, a3), r, types; coords=coords)
   ```

   One-line delegation; matches Spacey's existing `pointGroup_robust(u, v, w)` convention. The matrix form stays canonical and owns all validation; the vector form just `hcat`s. `a1, a2, a3` become the **columns** of `A`.

4. **Final struct shape:**

   ```julia
   struct Crystal{T}
       A::Matrix{Float64}          # 3×3 lattice, columns = a1, a2, a3
       r::Matrix{Float64}          # 3×N fractional positions, columns = atoms
       types::Vector{T}            # length N, any type comparable with ==
       function Crystal{T}(A::AbstractMatrix{<:Real}, r::AbstractMatrix{<:Real},
                           types::AbstractVector{T}; coords::Symbol) where T
           size(A) == (3, 3) || error("A must be 3×3")
           size(r, 1) == 3 || error("r must have 3 rows")
           size(r, 2) == length(types) || error("r columns must match types length")
           coords ∈ (:fractional, :cartesian) || error("coords must be :fractional or :cartesian")
           A64 = Float64.(A)
           r64 = coords === :cartesian ? inv(A64) * Float64.(r) : Float64.(r)
           new{T}(A64, r64, collect(types))
       end
   end
   Crystal(A, r, types::AbstractVector{T}; coords) where T = Crystal{T}(A, r, types; coords=coords)
   Crystal(a1, a2, a3, r, types; coords) = Crystal(hcat(a1, a2, a3), r, types; coords=coords)
   ```

   (Final details — exact signatures, whether to wrap positions into `[0, 1)` at construction, etc. — to be settled during Phase 1 implementation.)

### 4.3 Tolerance model

Point-group `tol` operates on normalized (unit-volume) lattice vectors. For positions we need a separate tolerance — atoms don't live at the lattice-vector scale.

**Options:**
- **A. Single `tol`.** Applied to both. Simple API but conflates two scales.
- **B. Separate `lattice_tol`, `pos_tol`.** Explicit but verbose.
- **C. `tol` auto-propagated.** Silent coupling; hard to debug.

**Recommendation:** **B**, with `pos_tol` defaulting to `lattice_tol` if unspecified.

**Decision (2026-04-23):**

1. **Two tolerances:** separate `lattice_tol` and `pos_tol` kwargs. They measure different things at different scales; conflating them would be confusing to document and debug.
2. **`lattice_tol` default:** `0.01`, matching existing `pointGroup_robust`. Operates on normalised (unit-volume) lattice vectors.
3. **`pos_tol` default:** `0.01 · (V/N)^(1/3)` where `V = |det(A)|` and `N` is the number of atoms. Unit-agnostic (scales with whatever units the user chose) and self-calibrating (scales with atom density). Exposed as a public helper `default_pos_tol(c::Crystal)` so users can inspect the value before calling `spacegroup`.
4. **α = 0.01 specifically.** The alternative α = 0.05 was considered and rejected because it silently over-promotes BaTiO₃-class ferroelectric distortions (Ti displaced by ~0.05 Å from cubic centre, just below a 0.1 Å tolerance). The asymmetry that decides the default: over-promotion is silent, under-promotion is obvious to the user; the default should err toward the recoverable failure mode. Full reasoning, counter-examples for both directions, and a regime-vs-α table are preserved in `designDiscussions.md`.

Signatures sketch:
```julia
spacegroup(c::Crystal; lattice_tol=0.01, pos_tol=default_pos_tol(c))
default_pos_tol(c::Crystal) = 0.01 * (abs(det(c.A)) / size(c.r, 2))^(1/3)
```

### 4.4 Return type

**Options:**
- **A. `Vector{Tuple{Matrix{Int}, Vector{Float64}}}`.** Straightforward, untyped.
- **B. `Vector{NamedTuple{(:R, :τ)}}`.** Self-documenting.
- **C. A dedicated `SpacegroupOp` struct**, enabling methods like `op1 * op2`, `inv(op)`, `apply(op, r)`.

**Recommendation:** **C**. Small investment, major downstream benefit (composition, verification, application to test positions all become natural).

**Decision (2026-04-23):** **C — dedicated `SpacegroupOp` struct.** Method overloads (`*`, `inv`, callable `op(r)`, mod-1 `==`) fall out cleanly and will be needed in the test suite almost immediately. ~30 lines of declaration + methods. Full reasoning, comparison table, and gotchas (integer-inverse rounding, mod-1 hashing, mutability) preserved in `designDiscussions.md`.

Struct skeleton:
```julia
struct SpacegroupOp
    R::Matrix{Int}
    τ::Vector{Float64}
end
```
Return type: `Vector{SpacegroupOp}`.

### 4.5 Cartesian vs. lattice representation of the returned ops

Point-group code returns both `LG` (lattice) and `G` (Cartesian). For space groups, the natural choice is fractional (R integer, τ rational).

**Recommendation:** Fractional only for the main return. Provide `toCartesian(op, A)` helper.

**Decision (2026-04-23):** Fractional only. `SpacegroupOp.R` stays integer-matrix in lattice coords; `SpacegroupOp.τ` stays fractional. A helper function converts to Cartesian on demand:

```julia
toCartesian(op::SpacegroupOp, A::Matrix{Float64}) =
    (A * op.R * inv(A), A * op.τ)   # returns (R_cart, τ_cart)
```

Rationale: the fractional form is exact (integer R, rational τ when snapped), matches ITA convention, and is what the algorithm uses internally. Cartesian is derivable and use-case-specific (e.g. visualisation, comparison with external tools that speak Cartesian); keeping it out of the stored op avoids the sync-drift and numerical-noise issues of storing both forms.

### 4.6 Stub transition

Current `spacegroup(c::Crystal)` returns `true` unconditionally.

**Options:**
- **A.** Replace with `error("not yet implemented")` in Phase 1; replace with real implementation in Phase 2.
- **B.** Leave stub until Phase 2 lands.
- **C.** Delete stub now, re-add in Phase 2.

**Recommendation:** **A**. Clearest signal of WIP status; prevents silent misuse in the intervening commits.

**Decision (2026-04-23):** **A.** Phase 1 replaces the current `spacegroup(c::Crystal) = true` stub with `error("spacegroup: not yet implemented — scheduled for Phase 2")`. Phase 2 replaces it with the real implementation in the same commit that lands the algorithm.

---

## 5. Risks and open questions

- **Position tolerance at symmetry special positions.** An atom near but not at a high-symmetry position produces over-/under-report ambiguity analogous to the lattice near-boundary case. A `verify_stable`-style tool for crystals is likely needed in Phase 4.
- **Non-primitive cells.** If the user supplies a conventional (non-primitive) cell, the lattice point group is too small — conventional cell misses the centering translations. We should either auto-primitize or require primitive input. `plan.md` §3.7 discusses the lattice-side version.
- **Rational `τ` recognition.** Numerically, `τ = 1/2` may come out as `0.5000000001`. For presentation and equality-testing we probably want to round/snap to small rationals. The `snapToSymmetry_SVD` machinery is a partial precedent.
- **Multi-type crystals with one-of-a-kind atoms.** A single-atom probe works only if that atom has at least one same-type counterpart per symmetry op. For crystals with only one of some type (e.g. a defect), the probe atom must be chosen from a type that has ≥ 2 atoms; handle the edge case where every type has exactly 1 atom.
- **Numerical identity vs. almost-identity.** `R = I` with `τ ≈ 0` must distinguish from `R = I` with `τ = 1/2`. This is tolerance-sensitive; handled naturally by the enumeration, but worth a dedicated test.

---

## 6. Next step

Resolve the open design choices in §4 (each marked _pending_), then begin Phase 1.
