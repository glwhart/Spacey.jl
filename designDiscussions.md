# Design Discussions

A running record of design choices made during Spacey.jl development, capturing the reasoning behind each decision so future work (or future re-examination) has context beyond the terse annotations in the planning documents. Companion to `plan.md`, `research.md`, and `spacegroup_plan.md`; companion to `aiCodeExplanation.md` for syntax notes.

---

## Default `pos_tol` value for `spacegroup`

**Context (2026-04-23):** Choosing the default position-matching tolerance for the `spacegroup` function planned in `spacegroup_plan.md` §4.3. The form `pos_tol = α · (V/N)^(1/3)` was already agreed (scales correctly under any unit system the user provides); this discussion is about the value of α.

### The initial proposal

**α = 0.05** (5% of characteristic atom separation). Rationale: noisy X-ray data should still find a fairly high symmetry; loose is safer than tight.

At α = 0.05, `pos_tol ≈ 0.1 Å` for typical crystals.

### Counter-argument: α = 0.05 silently over-promotes real materials

**BaTiO₃ at room temperature** (the canonical ferroelectric):
- Lattice `a = 3.994 Å, c = 4.038 Å`
- Ti at `(½, ½, 0.512)` — off-centred by ~0.048 Å
- O₂ at `(½, ½, 0.023)` — off-centred by ~0.09 Å
- True space group: **P4mm (8 ops)**
- At α = 0.05: `pos_tol ≈ 0.12 Å` — all displacements below `pos_tol`, all atoms treated as on the cubic prototype positions → Spacey would report **Pm3̄m (48 ops)**. 6× over-promotion, silent.
- At α = 0.01: `pos_tol ≈ 0.023 Å` — Ti displacement (0.048 Å) above tol → correctly **P4mm**.

The same silent-wrong-answer failure mode applies to many important material classes:
- **Perovskite octahedral tilts** (SrTiO₃ below 105 K, CaTiO₃, most orthorhombic perovskites). O displacements 0.02–0.08 Å.
- **Jahn–Teller distortions** (LaMnO₃, KCuF₃, Cu²⁺ salts). 0.05–0.2 Å octahedral distortion.
- **Peierls / dimerisation** (1D CDW materials, doped TCNQ). 0.03–0.08 Å pairing.
- **Martensites / shape-memory alloys** (NiTi, β-brass). <0.1 Å tetragonal/monoclinic distortion from cubic austenite.
- **Antiferroelectric distortions** (PbZrO₃, NaNbO₃).

These are all small-displacement symmetry-breaking — exactly the class of problems for which a user runs a space-group finder in the first place.

### Counter-counter-argument: is α = 0.01 too loose too?

Yes, in different regimes. No single α is "correct" across all data types; the trade-off just shifts:

| α | pos_tol (typical) | Safe for | Silently wrong on |
|---|---|---|---|
| 0.05 | ~0.1 Å | rough/powder/H-atom data | ferroelectrics, perovskite tilts, JT, martensites |
| 0.01 | ~0.02 Å | well-refined X-ray, normal DFT | quantum paraelectrics, tight DFT, magnetostriction |
| 0.001 | ~0.002 Å | tight DFT, synchrotron single-crystal | routine X-ray noise → under-reports |

**Where α = 0.01 still over-promotes:**
- **Quantum paraelectrics** (SrTiO₃ below ~4 K; the question of whether it's truly cubic or weakly tetragonal hinges on sub-0.01-Å displacements).
- **Relaxor ferroelectrics** (Pb(Mg,Nb)O₃, Pb(Sc,Ta)O₃) — sub-0.01-Å off-centring from cubic prototype.
- **Well-converged DFT data** (modern tight-convergence settings reach 10⁻⁴ Å precision; meaningful displacements can be 0.005 Å).
- **Magnetostrictive distortions** (~10⁻³ Å) — marginal for Spacey since it doesn't do magnetic symmetry, but could matter in post-processed workflows.
- **Small-amplitude phases near a transition** — displacements go continuously to zero approaching Tc.
- **High-resolution synchrotron single-crystal data** (~0.001 Å ESDs for heavy atoms).

**Where α = 0.01 is too tight (would under-promote):**
- **Powder X-ray data** (typical ESDs 0.01–0.05 Å).
- **H-atom positions from X-ray** (uncertainty 0.05–0.1 Å).
- **High-temperature structures** with large thermal motion (`U_eq` ~ 0.05 Å at 1000 K can shift refined centroids).
- **Early-iteration DFT relaxations** before full convergence.
- **Disordered / split-site refinements.**

### The asymmetry that decides the default

Over-promotion is **silent** — user sees what they expected, never questions it. A wrong space group is quietly returned and propagates into whatever downstream analysis follows.

Under-promotion is **loud** — user sees fewer operations than expected, investigates, loosens `pos_tol`, and reaches the correct answer.

Defaults should err toward the obvious, recoverable failure mode rather than the silent one.

### Decision

**α = 0.01** (`pos_tol ≈ 0.02 Å` for typical crystals).

Rationale:
- The class of structures silently mis-identified at α = 0.01 (quantum paraelectrics, tight DFT, magnetostriction) is narrower than the class silently mis-identified at α = 0.05 (ferroelectrics, perovskite tilts, JT, martensites).
- Users operating in the "too tight" regime (powder, H-atom, high-T) get loud under-reporting, which self-corrects with one `pos_tol` override.
- Users operating in the sub-0.02-Å regime tend to be experts who know to pass an explicit tol anyway.

### Future polish (optional)

A `regime` convenience keyword could map to preset α values without forcing users to remember specific numbers:
```julia
spacegroup(c; regime=:xray)      # α = 0.01   (default)
spacegroup(c; regime=:loose)     # α = 0.05   (powder, H-atom, rough data)
spacegroup(c; regime=:dft)       # α = 0.001  (well-converged theory)
```

Not needed for Phase 1; can be added if users ask for it.

---

## Return type for `spacegroup` — tuple vs NamedTuple vs dedicated struct

**Context (2026-04-23):** `spacegroup_plan.md` §4.4. Three candidates for what `spacegroup(c::Crystal)` should return:

- **A. `Vector{Tuple{Matrix{Int}, Vector{Float64}}}`** — bare tuple per op.
- **B. `Vector{NamedTuple{(:R, :τ)}}`** — named fields, no type.
- **C. `Vector{SpacegroupOp}`** — dedicated struct type.

### What option C looks like

```julia
struct SpacegroupOp
    R::Matrix{Int}        # 3×3 integer rotation (lattice coords)
    τ::Vector{Float64}    # length-3 fractional translation
end
```

Data-wise identical to the tuple — Julia structs with concrete field types compile to the same memory layout. The only cost is the declaration.

### What the name unlocks

Method overloading. Because `SpacegroupOp` is a nominal type, Julia's multiple dispatch lets us attach operations cleanly — which you *can't* do with tuples or NamedTuples without type piracy (defining methods on `Base.Tuple` affects every tuple in the program).

The methods that come naturally on a space-group operation:

```julia
# Composition: (R1, τ1) ∘ (R2, τ2) = (R1·R2, R1·τ2 + τ1)
Base.:*(a::SpacegroupOp, b::SpacegroupOp) =
    SpacegroupOp(a.R * b.R, a.R * b.τ + a.τ)

# Inverse: (R, τ)⁻¹ = (R⁻¹, -R⁻¹·τ); |det R| = 1 so R⁻¹ is integer
Base.inv(op::SpacegroupOp) = begin
    Ri = round.(Int, inv(Float64.(op.R)))
    SpacegroupOp(Ri, -Ri * op.τ)
end

# Apply to a fractional position (callable struct)
(op::SpacegroupOp)(r) = mod.(op.R * r + op.τ, 1.0)

# Identity op
Base.one(::Type{SpacegroupOp}) = SpacegroupOp(I(3), zeros(3))

# Equality with mod-1 semantics on τ
function Base.:(==)(a::SpacegroupOp, b::SpacegroupOp; atol=1e-8)
    a.R == b.R || return false
    Δ = mod.(a.τ - b.τ .+ 0.5, 1.0) .- 0.5
    all(abs.(Δ) .< atol)
end
```

### What these enable cleanly

- **Group closure verification:** `all(op1 * op2 ∈ ops for op1 in ops, op2 in ops)`. Without `*`, you'd write `(op1[1] * op2[1], op1[1] * op2[2] + op1[2])` inline every time.
- **Applying the full group to a position:** `[op(r) for op in ops]`. The callable-struct syntax reads like crystallographic notation.
- **Membership / duplicate detection:** `op ∈ ops` requires `==`, which can't be overridden on tuples.
- **Composition tables, conjugacy classes, coset decompositions** — all standard group-theoretic constructions fall out of `*` and `inv`.

### Future extensibility

The struct is a natural home for future methods without cluttering module-level namespaces:
- `is_symmorphic(op) = all(op.τ .≈ 0)`
- `seitz(op)` — Seitz notation string like `{4+ 0 0 1 | 1/2 1/2 1/2}`
- `snap_rational(op; denoms=1:12)` — round τ to nearest small rational
- Hashing for `Set{SpacegroupOp}` / `Dict` keys

None Phase-1-required, but opting into the struct today lets them drop in cleanly later.

### Ecosystem precedent

Standard pattern in crystallographic libraries:
- **pymatgen**: `SymmOp`
- **ASE**: `Spacegroup` with `.operations` returning typed ops
- **Spglib** (C API): struct-based operation representation
- **Mathematica**: `SpaceGroupOperation`

### Gotchas (solvable, but worth flagging)

- **`inv` of an integer matrix.** Julia's generic `inv` returns `Float64`. For space-group R with `det = ±1` the true inverse is integer; we round and (should) verify. Minor precision care needed.
- **Hashing.** If users want `Set{SpacegroupOp}` or `Dict` keys, we need `Base.hash`. Should canonicalise τ at construction (e.g. `mod.(τ, 1.0)` with a tolerance-based snap to 0) so equal ops hash identically.
- **Equality on τ.** Periodic-boundary semantics (mod 1) require the custom `==` above; the default field-by-field comparison would miss physically-equivalent ops with `τ = 0.0` vs `τ = 1.0`.
- **Immutability.** Use `struct`, not `mutable struct`. Ops are values, not objects with identity. Enables safe use as Dict keys, safe parallelism.

### Comparison

| | A. Tuple | B. NamedTuple | C. `SpacegroupOp` |
|---|---|---|---|
| Field access | `op[1]`, `op[2]` | `op.R`, `op.τ` | `op.R`, `op.τ` |
| Destructuring | `R, τ = op` | `(; R, τ) = op` | `(; R, τ) = op` |
| Method overloading | ✗ (type piracy) | ✗ (type piracy) | ✓ clean |
| Composition syntax | manual | manual | `op1 * op2` |
| Application syntax | manual | manual | `op(r)` |
| Custom `==` | ✗ | ✗ | ✓ |
| Pretty printing | default | default (verbose) | customisable |
| Future extensibility | poor | poor | excellent |
| Lines of code | 0 | 0 | ~30 |

### Decision

**C.** The composition / inverse / application methods will be needed in the test suite almost immediately (verify closure, verify inverses are in the group, apply each op to probe positions). Having them as first-class overloads on a typed op pays off in every subsequent file that touches space-group ops. The ~30-line cost is negligible and invested once.

---

## Canonicalising `τ` in `SpacegroupOp` — round-to-digits vs snap-to-rational

**Context (2026-04-24):** Followed §4.1/§4.4 once we started running the full AFLOW corpus and ran into Test failures on trigonal (P3₁21) structures. The issue surfaced as closure tests (`a * b ∈ ops`) failing for the 3-fold-screw structures (γ-Se, α-Quartz, Cinnabar, CrCl₃) even though the operation set was clearly correct.

### What broke

Original implementation:
```julia
SpacegroupOp(R, τ) = new(R, mod.(round.(τ, digits=10), 1.0))
```

The intent was to kill floating-point drift accumulated through basis transforms and op composition: `mod` folds `τ` into `[0, 1)` and `round(…, digits=10)` clamps drift below the 10⁻¹⁰ scale. This worked for τ values like `0`, `½`, `¼` — round-trip-stable as 10-digit decimals.

It silently broke for thirds:
- `round(1/3, digits=10) = 0.3333333333` (10 threes)
- `round(2/3, digits=10) = 0.6666666667` (rounded *up*: 6,6,...,6,7)
- But `0.3333333333 + 0.3333333333 = 0.6666666666` (ten 6s, no carry)

Composing two ops with `τ_a = τ_b = 1/3` produced `τ_result = 0.6666666666` after rounding, which doesn't match the canonical 10-digit form of 2/3 (`0.6666666667`). So `(a * b)` wasn't strictly `==` to the op in the group with τ = 2/3, even though both represented the same physical translation. `∈` returned `false`. Closure failed.

This was completely silent on cubic / tetragonal / orthorhombic / hexagonal corpus prototypes — none of those have ⅓ τ-components.

### The fix

Replace the round-to-digits with a snap-to-rational with small denominator:

```julia
function _canonicalise_τ(τ; tol=1e-6)
    out = similar(τ, Float64)
    for i in eachindex(τ)
        x = mod(Float64(τ[i]), 1.0)
        snapped = x
        for q in 1:12
            p = round(Int, q * x)
            cand = (p == q) ? 0.0 : p / q
            if abs(x - cand) < tol
                snapped = cand
                break
            end
        end
        out[i] = snapped
    end
    return out
end
```

Every τ component in an ITA space group is an exact rational with denominator ≤ 12 (½, ⅓, ¼, ⅙, 1/12, …). Snapping to the nearest such value within 1e-6 catches every standard τ exactly while still tolerating accumulated floating-point drift up to 1e-6 — which is well above the worst case observed in practice.

### Why `tol = 1e-6` and `q ≤ 12`

The ITA standard τ denominators top out at 12 (P6₃ has c/2, P6₅ has c/6, etc.). Going higher (q ≤ 24, say) would catch nothing more in the AFLOW corpus and would risk falsely snapping a generic τ to a nearby rational accident. 1e-6 is generous enough for any drift accumulated through `O(20)` matrix multiplications at Float64 precision (machine eps × 10⁵ ≈ 2 × 10⁻¹¹) but tight enough to never confuse a "real" non-rational τ with a rational one.

### Decision and rollout

Replaced `round(…, digits=10)` with snap-to-rational in commit `c2ac9bf` (AFLOW 3-fold-screw cases). Every existing test continued to pass, all four trigonal-screw structures now pass closure tests, and the full AFLOW corpus run uses the same canonicalisation. No exposed knob — `tol` is hardcoded to 1e-6 inside `_canonicalise_τ`. Could be made a kwarg later if needed.

---

## Crystal-system identification — Layer 1 only (holohedry from point-group order)

**Context (2026-04-24):** User suggested that adding a Bravais-lattice identification would strengthen the AFLOW corpus tests by providing a second invariant alongside the operation count. Distinct lattices can collide on op-count alone (e.g. simple cubic / FCC / BCC primitive cells all give 48). A Bravais-or-system check would catch lattice mis-identification that the op-count test misses.

### Two layers considered

| Layer | Output | Cost | Coverage |
|---|---|---|---|
| 1 | Crystal **system** (7 classes: tri / mono / ortho / tet / trig / hex / cubic) | Trivial — wrapper over `pointGroup_robust` | Catches lattice over-promotion (`a/b ≈ 1`) and lattice mis-identification within a system |
| 2 | Full **Bravais lattice** (14 classes: P / C / I / F variants of each system) | Non-trivial — needs Niggli or metric-tensor classification with edge cases | Catches centering mis-identification on top of Layer 1 |

### Why Layer 1 alone is enough for now

Each crystal system's holohedry has a unique point-group order:

```
order  →  system          holohedry
  2       triclinic       C_i
  4       monoclinic      C_2h
  8       orthorhombic    D_2h
 12       trigonal        D_3d
 16       tetragonal      D_4h
 24       hexagonal       D_6h
 48       cubic           O_h
```

So `crystal_system(A)` is just: run `pointGroup_robust` on `minkReduce(A)` and look up the order. ~10 lines of code, trivially correct.

The AFLOW prototype label's Pearson symbol encodes the expected system directly:
```
a*  →  triclinic    m*  →  monoclinic   o*  →  orthorhombic
t*  →  tetragonal   c*  →  cubic
hP  →  hexagonal    hR  →  trigonal
```

Cross-checking these two on the full AFLOW corpus (1095 prototypes) yields a different deviation set than the op-count test — exactly the value the user predicted. Some structures pass op-count but fail crystal-system (lattice over-promotion in Spacey but the atom count happens to come out right), and some pass crystal-system but fail op-count (atomic over-promotion or accidental symmetry, but lattice geometry matches the label). Together they catch more.

### Why Layer 2 is deferred

Going from "crystal system" (7 classes) to "full Bravais lattice" (14 classes) requires distinguishing P / C / I / F primitive cells *within* a system. For primitive cells (which AFLOW POSCARs use), the centering is reflected in the geometric pattern of the basis vectors — e.g. BCC primitive has equal-length vectors with arccos(1/3) inter-vector angles. The decision tree is roughly 30 branches over the metric tensor `G = AᵀA`, with a few subtle cases at boundaries (a hexagonal lattice with c/a coincidentally matching FCC's primitive ratio, etc.).

That's an afternoon of careful work, useful but not yet motivated. The Layer-1 cross-check has empirically caught a useful fraction of deviations; if more granularity is needed, Layer 2 can be added without breaking the existing API (just a new function `bravais(A)`).

### Decision

Layer 1 only, landed in commit `d8d1fbc`. Layer 2 deferred indefinitely.
