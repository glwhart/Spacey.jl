# Phase 2 Implementation Plan

Working plan for Phase 2 of the space-group extension (per `spacegroup_plan.md`). Scope: implement `spacegroup(c::Crystal; ...)` returning a `Vector{SpacegroupOp}` of the symmetry operations of crystal `c`, plus the `SpacegroupOp` struct with its core method overloads.

**Phase 2 exit criterion (per `spacegroup_plan.md`):** simple cubic, one atom at origin returns 48 ops, all with `τ = 0`.

**What this single test exercises.** A simple cubic lattice (`A = I`) with a single atom at the origin is the minimum-complexity crystal that touches every part of the algorithm. Its space group is **Pm3̄m (type #221), order 48** — the maximum 3D point group (cubic Oh, 48 rotations) combined with one atom at the origin. Every rotation preserves the atom, so `τ = 0` for all 48 ops. Passing this confirms four things simultaneously: (a) the loop over point-group candidates reaches all 48; (b) the `minkReduce` + basis-transform round-trip leaves an already-reduced lattice unchanged without drift; (c) `isSpacegroupOp` correctly validates every candidate; (d) we don't duplicate or miss ops. A bug in any of the four would show up as a wrong count or wrong τ-vector.

**This is a 4-eye-check document.** Add comments inline (prefixed `> GWH:` or similar) and I'll revise before touching code. Questions for you are marked **❓**.

---

## 1. Files to touch

- **`src/Spacey.jl`**
  - Add `struct SpacegroupOp` with field-equality, multiplication, inverse, callable, `one`, `show`.
  - Add `toCartesian(op, A)` helper.
  - Add `Base.hash(op::SpacegroupOp)` **(open question — see §6.1)**.
  - Implement `spacegroup(c::Crystal; lattice_tol, pos_tol)` — replaces the stub from Phase 1.
  - Extend `export` list: `SpacegroupOp`, `toCartesian`.
- **`test/runtests.jl`**
  - New `@testset "SpacegroupOp methods"` — tests struct overloads directly.
  - New `@testset "spacegroup: Phase 2 core cases"` — exit criterion + a few more sanity cases (simple tetragonal, simple orthorhombic, triclinic, CsCl).
- **`CLAUDE.md`**
  - Remove the "`spacegroup(A)` returns `true` — not implemented" line.
  - Remove "`Crystal` struct is defined but unused" line.
  - Add a one-paragraph description of the new API.
- **`spacegroup_plan.md`**
  - Annotate Phase 2 section with the commit hash once landed.

No other files touched.

---

## 2. `SpacegroupOp` struct and methods

Per §4.4/§4.5 decisions, copied here for reference and with concrete implementation sketches. I believe this is mostly direct transcription of the `designDiscussions.md` skeleton, but double-check the subtle bits.

```julia
struct SpacegroupOp
    R::Matrix{Int}         # 3×3 integer rotation in lattice coords
    τ::Vector{Float64}     # length-3 fractional translation, canonicalised to [0, 1)
    SpacegroupOp(R, τ) = new(R, mod.(τ, 1.0))
end
```

**Inner constructor canonicalises `τ`** (per §6.1 decision). Any input `τ` is folded into `[0, 1)` via `mod.(τ, 1.0)` before storage. Consequence: ops constructed with `τ = [1.0, 0, 0]`, `τ = [0.0, 0, 0]`, or `τ = [2.5, 0, 0]` all produce the same stored representation. This lets Julia's default field-by-field `==` and `hash` behave correctly without a custom method (see §2.1 for what this removes).

### 2.1 Methods

```julia
# Composition: (R1, τ1) ∘ (R2, τ2) = (R1·R2, R1·τ2 + τ1)
# Convention: op1 * op2 means "apply op2 first, then op1" (function-composition
# semantics, so (f ∘ g)(x) = f(g(x))). Derivation:
#     op2: r        ↦ R2·r + τ2
#     op1: R2·r+τ2  ↦ R1·(R2·r + τ2) + τ1  =  R1·R2·r + R1·τ2 + τ1
# Hence (op1 ∘ op2) = (R1·R2, R1·τ2 + τ1).
Base.:*(a::SpacegroupOp, b::SpacegroupOp) =
    SpacegroupOp(a.R * b.R, a.R * b.τ + a.τ)

# Inverse: (R, τ)⁻¹ = (R⁻¹, -R⁻¹·τ)
# |det R| = 1 for any valid lattice rotation, so R⁻¹ is integer
function Base.inv(op::SpacegroupOp)
    Rinv_f = inv(Float64.(op.R))
    Rinv = round.(Int, Rinv_f)
    # Sanity: verify rounding was exact
    maximum(abs, Rinv_f .- Rinv) < 1e-8 || error("inverse of R is not integer — det(R) ≠ ±1?")
    SpacegroupOp(Rinv, mod.(-Rinv * op.τ, 1.0))
end

# Apply to a fractional position vector (callable struct)
(op::SpacegroupOp)(r::AbstractVector) = mod.(op.R * r + op.τ, 1.0)

# Identity op
Base.one(::Type{SpacegroupOp}) = SpacegroupOp(Matrix{Int}(I, 3, 3), zeros(3))

# Equality: no custom method needed — with τ canonicalised to [0,1) at
# construction, Julia's default field-by-field == (and hash) do the right thing.
# Trade-off: two ops whose τ differ by floating-point epsilon (e.g. 0.5 vs
# 0.5 + 1e-15) would compare unequal. In Phase 2 we never construct such
# pairs (τ comes from atom differences or integer-matrix transforms, both
# exact), so this is acceptable. If a downstream need arises, add a
# tolerance-based `≈`/`isapprox` method alongside `==`.

# REPL printing — never called directly by user code. Julia invokes it
# automatically whenever it needs to display a SpacegroupOp: REPL echo,
# println, string interpolation, vector printing, etc. Defining `show`
# is the standard Julia way to control how a custom struct appears.
Base.show(io::IO, op::SpacegroupOp) =
    print(io, "SpacegroupOp(R = ", op.R, ", τ = ", op.τ, ")")

# Cartesian conversion (returns a tuple, not a SpacegroupOp)
toCartesian(op::SpacegroupOp, A::AbstractMatrix) =
    (A * op.R * inv(A), A * op.τ)
```

### 2.2 Design notes worth flagging

- ~~**`==` does not compare hash-equivalently.**~~ **Resolved:** canonicalising `τ` to `[0, 1)` at construction (§6.1 decision) makes default `==` and `hash` consistent. The only case this doesn't handle is floating-point-epsilon differences in `τ` — which Phase 2 code paths don't produce. Flagged as a possible later addition if downstream code needs a tolerance-based `isapprox`.
- **`inv` requires `det(R) = ±1`.** True for any rotation in a lattice point group, but we verify with an assertion for safety. A `SpacegroupOp` with `det(R) = 2` would not satisfy this — but we never construct such in `spacegroup()`, since point-group ops have `|det| = 1`. **❓ Should we also validate this at `SpacegroupOp` construction?** (I lean no — cost of a check per construction, and we construct only internally.) GH: No, lets not. **Decided: no construction-time det check.**
- **`toCartesian` returns a `Tuple`, not a `SpacegroupOp`.** Because the Cartesian `R` is a `Matrix{Float64}`, not `Matrix{Int}`. A `SpacegroupOp` is defined to hold integer `R` (lattice coords); Cartesian output is deliberately a different shape.

---

## 3. `spacegroup` algorithm — full walk-through

The math is the Phase 2 sketch from `spacegroup_plan.md` §2.1. The subtle part is the **basis-transformation dance** required because `pointGroup_robust` expects a Minkowski-reduced input, but the atoms are given in the user's (possibly non-reduced) basis.

### 3.1 Why we need to dance

`pointGroup_robust(u, v, w; tol)` errors if the input is not Mink-reduced (line ~216 of `src/Spacey.jl`). The user's `c.A` is whatever they passed in — generally not guaranteed reduced. So we must:

1. Compute `A_red = minkReduce(c.A)`.
2. Transform atoms to the reduced basis: `r_red = U_or * c.r`, where `U_or = A_red⁻¹ · c.A` is the (integer!) change-of-basis.
3. Find ops in the reduced basis: `pointGroup_robust(A_red_cols...)` → `{R_red}`; enumerate `τ_red` per R.
4. Transform ops back to the user's basis before returning: `R_orig = U_ro · R_red · U_or`, `τ_orig = U_ro · τ_red`, where `U_ro = inv(U_or)`.

The matrices `U_or` and `U_ro` are integer unimodular (`|det| = 1`) because `minkReduce` only does integer-valued basis operations (permutations and integer-coefficient subtractions, see `MinkowskiReduction.jl` source).

The returned ops are in the **user's original basis** — so `op(c.r[:, i])` gives an orbit position in the user's fractional coordinates, consistent with what they passed in.

### 3.2 Pseudocode

```julia
function spacegroup(c::Crystal; lattice_tol::Real=0.01,
                                pos_tol::Real=default_pos_tol(c))
    # --- 3.2.1 Minkowski-reduce the lattice ------------------------------
    # GH: Why not use the matrix version of minkReduce?  [Using it — cleaner.]
    A_red = minkReduce(c.A)
    u_red, v_red, w_red = eachcol(A_red)

    # --- 3.2.2 Change-of-basis integer matrices --------------------------
    # A_orig · U_ro = A_red    (so U_ro: reduced-coord → original-coord)
    # A_red  · U_or = A_orig   (U_or = inv(U_ro): original-coord → reduced-coord)
    U_ro = round.(Int, inv(c.A) * A_red)
    U_or = round.(Int, inv(A_red) * c.A)
    # Sanity: unimodular + round-trip identity
    abs(det(U_ro)) == 1 || error("change-of-basis not unimodular: det = $(det(U_ro))")
    U_ro * U_or == I(3) || error("change-of-basis round-trip failed")
    norm(c.A * U_ro - A_red) < 1e-8 * opnorm(c.A) ||
        error("reduced basis does not agree with integer transform of original")

    # --- 3.2.3 Transform atoms to reduced basis --------------------------
    r_red = mod.(U_or * c.r, 1.0)
    c_red = Crystal(A_red, r_red, c.types; coords=:fractional)

    # --- 3.2.4 Point group of the lattice (reduced basis) ----------------
    LG_red, _ = pointGroup_robust(u_red, v_red, w_red; tol=lattice_tol)

    # --- 3.2.5 Per-R, enumerate candidate τ via probe-atom differences ---
    # Pick the atom type with the FEWEST atoms as the probe (smallest
    # candidate set). Ties broken arbitrarily.
    types_unique = unique(c.types)
    counts = [count(==(t), c.types) for t in types_unique]
    probe_type = types_unique[argmin(counts)]
    probe_indices = findall(==(probe_type), c.types)
    i0 = probe_indices[1]

    # --- 3.2.6 Collect surviving (R_red, τ_red) --------------------------
    ops_red = Tuple{Matrix{Int}, Vector{Float64}}[]
    for R in LG_red
        image_i0 = R * c_red.r[:, i0]
        for j in probe_indices
            τ = mod.(c_red.r[:, j] .- image_i0, 1.0)
            if isSpacegroupOp(R, τ, c_red; tol=pos_tol)
                push!(ops_red, (Matrix{Int}(R), τ))
            end
        end
    end

    # --- 3.2.7 Transform ops back to user's basis ------------------------
    ops_out = [SpacegroupOp(U_ro * R_red * U_or, U_ro * τ_red)
               for (R_red, τ_red) in ops_red]

    # --- 3.2.8 Sort identity to the front (§6.3 decision) ----------------
    # The identity (R=I, τ=0 mod 1) is guaranteed to be present because it
    # is always a symmetry of any crystal. Put it at index 1 so users can
    # rely on `spacegroup(c)[1]` being the identity.
    e = one(SpacegroupOp)
    id_idx = findfirst(==(e), ops_out)
    id_idx === nothing && error("identity missing from computed space group — bug")
    id_idx != 1 && (ops_out[1], ops_out[id_idx] = ops_out[id_idx], ops_out[1])
    return ops_out
end
```

### 3.3 Complexity

- Minkowski reduction: O(1) (bounded-iteration, trivial for small 3D).
- Change-of-basis: O(1).
- Point group: O(1) (27-neighbor candidate search).
- **Per-R candidate τ enumeration: O(N_probe) candidates, each tested with `isSpacegroupOp` which is O(N²).**
- **Total: O(|LG| · N_probe · N²) ≤ O(48 · N · N²) = O(N³).**

For typical crystals (N ≤ 100), this is < 10⁷ ops — trivial. For larger crystals (N ~ 1000+), the `isSpacegroupOp` inner loop becomes the bottleneck — Phase 4 territory.

### 3.4 Correctness argument

Three claims need to hold for the algorithm to be complete and sound:

1. **Completeness of R_red.** `pointGroup_robust` returns all point-group operations of the *lattice*. Every space-group op's `R` is in this set (a space-group op preserves the lattice, so its rotation preserves it too). ✓

2. **Completeness of τ candidates.** If `(R, τ)` is a symmetry, then the probe atom `i0` maps to some atom `j` of the same type: `R·r_{i0} + τ ≡ r_j (mod 1)`. So `τ ≡ r_j − R·r_{i0} (mod 1)` for some `j` of the same type. Our enumeration covers exactly this set. ✓

3. **Verification via `isSpacegroupOp`.** A candidate `(R, τ)` is accepted iff every atom maps to an atom of the same type under the op, modulo the lattice, within `pos_tol`. This is the definition of a space-group op. ✓

The only thing not guaranteed by construction: that the `claimed` injectivity guard inside `isSpacegroupOp` handles edge cases correctly for loose tolerances — but this is the Phase 1 function we already tested.

---

## 4. Test plan

### 4.1 `@testset "SpacegroupOp methods"` — direct struct tests

These exercise the new type without invoking the finder. Handwritten examples.

```julia
@testset "SpacegroupOp methods" begin
    I3 = Matrix{Int}(I, 3, 3)

    # Identity
    e = one(SpacegroupOp)
    @test e.R == I3
    @test e.τ == zeros(3)

    # Composition
    R1 = [0 -1 0; 1 0 0; 0 0 1]  # 4-fold about z
    op1 = SpacegroupOp(R1, [0.5, 0.0, 0.0])
    op2 = SpacegroupOp(R1, [0.0, 0.5, 0.0])
    op12 = op1 * op2
    @test op12.R == R1 * R1
    @test op12.τ ≈ R1 * [0.0, 0.5, 0.0] + [0.5, 0.0, 0.0]

    # Inverse: op * inv(op) == e, mod 1
    op = SpacegroupOp(R1, [0.25, 0.0, 0.0])
    @test (op * inv(op)) == one(SpacegroupOp)
    @test (inv(op) * op) == one(SpacegroupOp)

    # Callable
    r = [0.1, 0.2, 0.3]
    @test op(r) ≈ mod.(R1 * r + op.τ, 1.0)

    # Equality mod-1 on τ
    opA = SpacegroupOp(I3, [0.0, 0.0, 0.0])
    opB = SpacegroupOp(I3, [1.0, 0.0, 0.0])   # equivalent mod 1
    opC = SpacegroupOp(I3, [0.5, 0.0, 0.0])
    @test opA == opB
    @test opA != opC

    # toCartesian round-trip on identity
    A = [1.0 0.5 0.0; 0.0 1.0 0.0; 0.0 0.0 1.5]
    Rc, τc = toCartesian(one(SpacegroupOp), A)
    @test Rc ≈ Matrix{Float64}(I, 3, 3)
    @test τc ≈ zeros(3)
end
```

### 4.2 `@testset "spacegroup: Phase 2 core cases"` — finder tests

```julia
@testset "spacegroup: Phase 2 core cases" begin
    # 4.2.1 Exit criterion: simple cubic, 1 atom at origin → 48 ops, all τ=0
    A_cubic = Matrix{Float64}(I, 3, 3)
    c_sc = Crystal(A_cubic, reshape([0.0, 0.0, 0.0], 3, 1), [:X]; coords=:fractional)
    ops = spacegroup(c_sc)
    @test length(ops) == 48
    @test all(op.τ ≈ zeros(3) for op in ops)
    @test all(abs(det(op.R)) == 1 for op in ops)

    # 4.2.2 Simple cubic, 1 atom shifted to body-centre (0.5, 0.5, 0.5):
    #       still 48 ops (origin shift preserves the abstract space group,
    #       though τ values change).
    c_sc_shift = Crystal(A_cubic, reshape([0.5, 0.5, 0.5], 3, 1), [:X]; coords=:fractional)
    @test length(spacegroup(c_sc_shift)) == 48

    # 4.2.3 Simple tetragonal (a=1, c=1.5), 1 atom at origin → 16 ops
    A_tet = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.5]
    c_tet = Crystal(A_tet, reshape([0.0, 0.0, 0.0], 3, 1), [:X]; coords=:fractional)
    @test length(spacegroup(c_tet)) == 16

    # 4.2.4 Simple orthorhombic (a, b, c all different), 1 atom → 8 ops
    A_ortho = [1.0 0.0 0.0; 0.0 1.2 0.0; 0.0 0.0 1.5]
    c_ortho = Crystal(A_ortho, reshape([0.0, 0.0, 0.0], 3, 1), [:X]; coords=:fractional)
    @test length(spacegroup(c_ortho)) == 8

    # 4.2.5 Triclinic (general), 1 atom → 2 ops (E, i)
    A_tri = [1.0 0.3 0.2; 0.0 1.1 0.4; 0.0 0.0 0.9]
    c_tri = Crystal(A_tri, reshape([0.0, 0.0, 0.0], 3, 1), [:X]; coords=:fractional)
    @test length(spacegroup(c_tri)) == 2

    # 4.2.6 CsCl: Cs at (0,0,0), Cl at (½,½,½) → 48 ops, all τ=0 (Pm3̄m)
    c_cscl = Crystal(A_cubic, [0.0 0.5; 0.0 0.5; 0.0 0.5], [:Cs, :Cl];
                     coords=:fractional)
    ops_cscl = spacegroup(c_cscl)
    @test length(ops_cscl) == 48
    @test all(op.τ ≈ zeros(3) for op in ops_cscl)

    # 4.2.7 Cross-check: every returned op passes isSpacegroupOp
    for op in ops
        @test isSpacegroupOp(op.R, op.τ, c_sc)
    end
    for op in ops_cscl
        @test isSpacegroupOp(op.R, op.τ, c_cscl)
    end

    # 4.2.8 Closure: op1 * op2 is also in the group (spot check)
    ops48 = spacegroup(c_sc)
    for _ in 1:5
        a, b = rand(ops48), rand(ops48)
        @test (a * b) ∈ ops48
    end

    # 4.2.9 Identity is always at index 1 (§6.3 decision)
    @test ops48[1] == one(SpacegroupOp)
    @test ops_cscl[1] == one(SpacegroupOp)
    @test spacegroup(c_tri)[1] == one(SpacegroupOp)

    # 4.2.10 Defensive regression test: an atom placed exactly on the
    # [0, 1) boundary is folded cleanly (see §7 item 5).
    c_boundary = Crystal(A_cubic,
        reshape([1.0, 0.0, 0.0], 3, 1),      # at x = 1.0 exactly, equivalent to x = 0
        [:X]; coords=:fractional)
    # ^ Note: the Crystal constructor doesn't currently fold, so this
    # tests that spacegroup is robust to the equivalent-position case.
    # Expected: same as c_sc (48 ops).
    @test length(spacegroup(c_boundary)) == 48
end
```

### 4.3 Tests I'm explicitly NOT adding in Phase 2

Deferred to **Phase 3** (per `spacegroup_plan.md` §3 Phase 3):

- NaCl (conventional FCC cell, 192 ops) — tests centering translations
- Diamond (Fd3̄m, 192 ops, non-symmorphic) — tests glide planes
- HCP (P6₃/mmc, 24 ops) — tests screw axes
- Quartz (P3₁21, 6 ops) — non-symmorphic low-symmetry

These all work with the same algorithm; they're in Phase 3 because that phase is about **testing against known crystal-structure tables**, and those crystals are the canonical benchmarks.

Deferred to **Phase 4** (robustness):
- Near-boundary crystals with noise.
- `verify_stable`-analog.
- Tolerance-sweep diagnostics.

---

## 5. CLAUDE.md update

Remove two stale lines:

```
- `spacegroup(A)` returns `true` — not implemented, reserved for future work
- `Crystal` struct is defined but unused
```

Replace with a short paragraph under "Exported API" (table):

| Function | Purpose |
|---|---|
| `Crystal(A, r, types; coords)` | Construct a crystal (fractional or Cartesian positions). |
| `spacegroup(c; lattice_tol, pos_tol)` | Find all `(R, τ)` symmetry operations. |
| `isSpacegroupOp(R, τ, c; tol)` | Check if `(R, τ)` is a symmetry of `c`. |
| `fractional(c)`, `cartesian(c)` | Atomic positions in either basis. |
| `default_pos_tol(c)` | Scale-invariant position tolerance default. |
| `toCartesian(op, A)` | Convert a `SpacegroupOp` to Cartesian form. |

Plus a one-liner: "Space-group ops are returned as `Vector{SpacegroupOp}` with integer `R` and fractional `τ`, both in the user's original basis."

---

## 6. Open questions

### 6.1 ❓ Should `SpacegroupOp` implement `Base.hash`?

Currently, two ops with `τ = [0.0, 0, 0]` and `τ = [1.0, 0, 0]` are `==` but would hash differently under Julia's default (field-by-field) hash. This means `Set{SpacegroupOp}` and `Dict{SpacegroupOp, _}` could have duplicate entries under the `==` semantics.

Fix: canonicalise `τ` in the constructor (`mod.(τ, 1.0)`), then default hash is consistent with `==`.
GH: Yes, do this.

**Lean:** canonicalise at construction. Cheap, closes a footgun. The only case this changes is if someone constructs an op with `τ = [2.5, 0, 0]` — we'd store `τ = [0.5, 0, 0]`. Any reasonable workflow treats these as the same op anyway. GH: just fine

❓ **Do you agree?** GH: Yes

**Decision (2026-04-23):** Canonicalise `τ` to `[0, 1)` via `mod.(τ, 1.0)` in an inner constructor. Consequence: we **drop the custom `==` method** and rely on Julia's default field-by-field `==` (and hash) — which is now consistent with the canonicalised representation. See §2.1 comment block for the trade-off (floating-point-epsilon differences in `τ` are treated as unequal).

### 6.2 ❓ Should `spacegroup` also return the Cartesian form for compatibility with `pointGroup_robust`'s `(LG, G)` return?

`pointGroup_robust` returns `(LG, G)` — integer-matrix + Cartesian-rotation. Our `spacegroup` returns only `Vector{SpacegroupOp}` (fractional/lattice-coord form), per §4.5.

Alternatives:
- **A.** Keep §4.5 decision: fractional-only; `toCartesian` helper is the escape hatch.
- **B.** Return both, matching `pointGroup_robust`'s `(LG, G)` interface style, as `(Vector{SpacegroupOp}, Vector{Tuple{Matrix{Float64}, Vector{Float64}}})`.

**Lean:** A. It's what §4.5 resolved. `toCartesian` handles the rare need for Cartesian form. But I want to flag that `pointGroup_robust`'s two-return style is a minor inconsistency in the package API — worth a note.

**Decision (2026-04-23):** Stick with **A** for Phase 2. `spacegroup` returns `Vector{SpacegroupOp}` (fractional only); `toCartesian` is the escape hatch.

**Follow-up question raised by GH:** *should `pointGroup_robust` itself be changed to return one type plus a Cartesian helper, mirroring this style?* Resolution: yes in the long run, but not in Phase 2 — see "Future cleanup" below.

**Future cleanup (deferred, not Phase 2):** refactor `pointGroup_robust` to return only `LG` (integer-matrix form in lattice coords) with a `toCartesian`-style helper for the `G` case. This is a breaking API change that touches every existing call site (internal Spacey, `runtests.jl`, `nearMissBoundary.jl`, any user code). The right place for it is a future "API consistency" pass (v0.8 or similar) bundled with the other small cleanups from `plan.md` §2.1–2.6. Noted here so it doesn't get lost.

### 6.3 ❓ What does `spacegroup` do when the crystal has minimal symmetry, and what order should the returned ops be in?

**Minimum symmetry of a Bravais lattice.** Every Bravais lattice — including triclinic — has inversion as a point-group symmetry. So a triclinic crystal with one atom at a generic position `r` has **2** space-group ops, not 1:
- `(I, 0)`
- `(-I, 2r mod 1)` — pure inversion about the origin maps `r → -r`, but `(-I, τ)` with `τ = 2r mod 1` maps `r → -r + 2r = r`, so it is a symmetry.

This is what the §4.2.5 test expects (`@test length(spacegroup(c_tri)) == 2`). For a crystal with truly **point group 1 (only identity)**, you need a chiral/asymmetric atomic arrangement that destroys inversion — e.g. multiple different-type atoms at generic positions with no inversion-related pairs. The algorithm handles this naturally: only `R = I` passes, only `τ = 0` passes.

(Earlier draft of this section incorrectly claimed a 1-atom triclinic crystal had only `{E}` — corrected per GH.)

**Order of the returned Vector.** Initial lean was unspecified order. Per GH ask, identity is sorted to index 1 instead — the cost is one linear scan post-collection, and `ops[1]` becoming a reliable identity reference is a useful sanity check for callers.

**Decision (2026-04-23):** put `(I, 0)` at index 1; leave the remaining order unspecified.

### 6.4 ❓ Should `spacegroup` accept a `verify_stable` kwarg?

The point-group finder has `verify_stable` (re-runs at tighter tol, warns on disagreement). An analogous flag for `spacegroup` would re-run at tighter `pos_tol` and warn if the operation count changes.

**Lean:** defer to Phase 4 (robustness). Phase 2 is about the core algorithm; Phase 4 adds the stability check for crystals specifically.

❓ **Agree?** GH: defer to phase 4

**Decision (2026-04-23):** Defer to Phase 4.

### 6.5 ❓ Probe-atom choice — min-count type or first atom?

I propose picking the atom *type* with the fewest atoms (smallest τ candidate set per R, so fewer isSpacegroupOp calls). Tie-broken by first-appearance.

Alternative: always use the first atom (simpler code).

**Lean:** min-count. The speedup is real when one type has 1 atom and another has 100; the added complexity is ~3 lines.

❓ **Agree?** GH: Agree

**Decision (2026-04-23):** Min-count probe type, tie broken by first-appearance.

---

## 7. Things that could go wrong (pre-mortem)

Anticipated failure modes during implementation:

1. **Integer rounding of `U_ro` / `U_or` fails for non-well-conditioned lattices.** If `inv(c.A) * A_red` has entries far from integer (say `0.5000001`), the `round.(Int, …)` produces the wrong integer. Sanity checks on `det(U_ro) == ±1` and `U_ro * U_or == I` catch this. Extreme aspect ratios are the place this would bite.
2. **`pointGroup_robust` disagrees with the true lattice point group** due to the tolerance choices we discussed in `plan.md`/`research.md`. If the user's lattice is near a Bravais boundary (tetragonal/cubic with ε distortion, etc.), `lattice_tol` controls how `pointGroup_robust` sees it, and `spacegroup` inherits that. Not a Phase 2 bug — it's a known robustness issue pushed to Phase 4 with a `verify_stable` analog.
3. **The `c_red` Crystal build fails the singular-A check** if numerical noise in `A_red` makes its det tiny. Unlikely — `minkReduce` output has the same det magnitude as input.
4. **Loose `pos_tol` over-promotes** analogous to the BaTiO₃-class case already discussed. Same mitigation: tight default + `verify_stable` in Phase 4.
5. **Off-by-one in `mod.(U_or * c.r, 1.0)`** at the `[0, 1)` boundary. Could show up as an atom at exactly `r = 1.0` not being folded to `0.0`. The `mod.()` in Julia handles this for positive values, but it's worth a defensive regression test — one that is not strictly needed if the code is obviously correct, but that pins this low-probability failure mode in case a future refactor touches the folding step. Implemented as test §4.2.10.

---

## 8. Commit plan

One commit:
- `src/Spacey.jl`: `SpacegroupOp` struct + methods, `toCartesian`, `spacegroup` implementation, exports updated.
- `test/runtests.jl`: two new testsets.
- `CLAUDE.md`: stale-line removals + new API blurb.
- `spacegroup_plan.md`: annotate Phase 2 as landed with the commit hash.

Commit message template:

```
Phase 2: spacegroup finder + SpacegroupOp type

Implements spacegroup(c; lattice_tol, pos_tol) returning the Vector of
SpacegroupOps that make up the crystal's space group. Algorithm:
Minkowski-reduce the lattice, find the point group via pointGroup_robust
on the reduced basis, enumerate candidate τ per R via probe-atom
differences, verify each candidate with isSpacegroupOp, and transform
surviving ops back to the user's basis.

Introduces SpacegroupOp with operator overloads for composition (*),
inverse, callable application to positions, mod-1 equality, identity
(one), pretty-printing, and the toCartesian helper.

Phase 2 exit criterion met: simple cubic with one atom returns 48 ops
all with τ=0. Also tested on simple tetragonal (16), simple orthorhombic
(8), triclinic (2), CsCl (48), closure, and isSpacegroupOp cross-check.
Non-symmorphic cases (diamond, HCP, quartz, NaCl conventional) deferred
to Phase 3.
```

---

## 9. Out of scope for Phase 2 (explicit)

- Classification into the 230 space groups (never in this plan's scope at all).
- Phase 3's known-crystal-structure comparison tests.
- Phase 4's noise-robustness + `verify_stable`-style detection.
- Performance optimization beyond the naive O(N³) matching.
- Hashing & Dict-key support — deferred pending §6.1 decision.
- Non-primitive cell auto-reduction (user must supply the cell they want analysed).
- Any kind of Wyckoff-position analysis or special-position detection.
