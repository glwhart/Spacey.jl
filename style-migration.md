# Style migration: snake_case + Julia style guide

Working notes for the v0.8 → v0.9 style sweep. Reply inline (✅ / ❌ / notes
under each row or item) and I'll execute from this file.

For convenience, where I list multiple options I've marked my recommendation
with **(rec)**. If a row has no ambiguity, the "Notes" column is blank.

> Let's treat "pointgroup" and "spacegroup" throughout as one word.

>isagroup=> is_group

---

## 🔸 Revisions in this pass (2026-05-15)

Marker: 🔸 = new/changed since last review. Stuff I still need from you is in
the **"Open for you"** block immediately below.

Resolved by your replies:
- `pointgroup` / `spacegroup` treated as one word → drives 6 renames in §1a/§1b.
- `isagroup` → `is_group` confirmed; options removed.
- Section 1c's "Where" column now uses the new (post-rename) names.

Other changes:
- §2a S13, S14, S15: my answers added inline.
- §2a S9 and §2b C1: I added a recommendation since you left them blank.
- §3 reframed for "all-at-once, no deprecation, release as v0.9.0".
- New §5: "Risk / known footguns of doing it all at once" — short.

### 🔸 Open for you

These are the only things I still need from you to start:

1. **§2a S9** — keyword-ify `pointgroup_simple`'s `debug` flag? (recommendation given)
> your rec
2. **§2a S10** — you wrote "∈"; I'm reading that as ✅ apply, but flag if you meant something else.
> yes, apply
3. **§2a S12** — same: "∈" → reading as ✅ apply.
> yes, apply
4. **§2b C1** — implicit vs explicit `return`? (recommendation given)
> your rec
5. **§2b C3** — I'll mark candidate docstring blocks during the edit, you review the diff. ✅ that's the workflow you wanted?
> yes

> Dont delete these:  debugging.jl, cubic_example.jl, snapExample.jl. Move them to `examples`. Make sure they still work.

### 🔸 Follow-up questions on the "move to examples" directive

Two things I need before I touch those three files:

6. **The three `2D_snap_example_*.jl` files at repo root** (`2D_snap_example_hex.jl`, `2D_snap_example_rectangle.jl`, `2D_snap_example_square.jl`) — these also look like examples. Move them into `examples/` too, or leave them at the root?
> 

7. **What does "make sure they still work" mean for `debugging.jl`?** Heads-up on what I'd be walking into:
   - `cubic_example.jl` and `snapExample.jl` reference `snapToSymmetry(...)` plain — a function that doesn't exist (current name is `snapToSymmetry_SVD`, returns 5 values, not 3). And `snapExample.jl` unpacks `pointGroup_robust` as a 2-tuple, which is the **pre-v0.8** API. So these aren't just out-of-date on naming — they're out-of-date on signatures. Fixable, but it's a real rewrite, not a rename.
   - `debugging.jl` is a 271-line scratch pad that re-defines `avgVecOverOps` and `snapToSymmetry_avg` locally (now in `Spacey` itself), interleaved with experimental snippets. "Still works" is ambiguous here.

   Pick one:
   - **(A)** Each file runs top-to-bottom under `julia --project examples/<file>.jl` without error. I update API calls and rewrite any block that can't be made to run.
   - **(B)** Each `include("examples/<file>.jl")` from the Spacey project loads without syntax errors — but I don't promise every individual snippet executes meaningfully.
   - **(C)** Just rename references to the new API (`pointgroup_robust` etc.); accept that already-broken snippets stay broken (i.e. preserve current behaviour, just snake-cased).
   - **(D)** Other — specify.
> 

Everything else is locked.

---

## 1. Renames

### 1a. Exported (public API) — 🔸 no deprecation aliases (per M1); v0.9.0 ships the rename

| Current             | Proposed                | Notes                                                                                                            |
|---------------------|-------------------------|------------------------------------------------------------------------------------------------------------------|
| `pointGroup`        | 🔸 `pointgroup`         | one word per your decision                                                                                       |
| `snapToSymmetry_SVD`| 🔸 `snap_to_symmetry_svd` | (chose option A; remaining options stripped — flag if you want B/C instead)                                       |
| `isagroup`          | 🔸 `is_group`           | confirmed                                                                                                        |
| `Crystal`           | *(keep)*                | already `UpperCamelCase`                                                                                         |
| `isSpacegroupOp`    | 🔸 `is_spacegroup_op`   | locked by "spacegroup as one word"                                                                               |
| `fractional`        | *(keep)*                |                                                                                                                  |
| `cartesian`         | *(keep)*                |                                                                                                                  |
| `default_pos_tol`   | *(keep)*                |                                                                                                                  |
| `crystal_system`    | *(keep)*                |                                                                                                                  |
| `SpacegroupOp`      | *(keep)*                | 🔸 locked — one-word `Spacegroup` + `Op`                                                                          |
| `toCartesian`       | `to_cartesian`          |                                                                                                                  |
| `spacegroup`        | *(keep)*                | 🔸 locked — already one word, already snake-case-compatible                                                       |
| `is_equiv_lattice`  | *(keep)*                |                                                                                                                  |
| `is_derivative`     | *(keep)*                |                                                                                                                  |
| `is_primitive`      | *(keep)*                |                                                                                                                  |
| `make_primitive`    | *(keep)*                |                                                                                                                  |
| `read_poscar`       | *(keep)*                |                                                                                                                  |

### 1b. Internal (reached via `Spacey.<name>`) — no deprecation needed

| Current                            | Proposed                            | Notes                                                                                |
|------------------------------------|-------------------------------------|--------------------------------------------------------------------------------------|
| `pointGroup_robust`                | 🔸 `pointgroup_robust`              | one word                                                                             |
| `pointGroup_fast`                  | 🔸 `pointgroup_fast`                | one word                                                                             |
| `pointGroup_simple`                | 🔸 `pointgroup_simple`              | one word                                                                             |
| `snapToSymmetry_avg`               | `snap_to_symmetry_avg`              |                                                                                      |
| `avgVecOverOps`                    | `avg_vec_over_ops`                  |                                                                                      |
| `aspectRatio`                      | `aspect_ratio`                      |                                                                                      |
| `threeDrotation`                   | 🔸 `rotate_basis_3d`                | (chose option A; flag if you want B/C — but it's test-only so impact is local)        |
| `_canonicalize_τ`                  | *(keep)*                            | 🔸 keeping unicode to match the `τ` field                                             |
| `_probe_atoms`                     | *(keep)*                            |                                                                                      |
| `_find_translations_for_rotation`  | *(keep)*                            |                                                                                      |
| `_find_self_translations`          | *(keep)*                            |                                                                                      |

### 1c. Local variables inside function bodies (no API impact)

| Current               | Proposed             | Where                                                                                              |
|-----------------------|----------------------|----------------------------------------------------------------------------------------------------|
| `inputVol`            | `input_vol`          | 🔸 `pointgroup_robust`                                                                              |
| `offDiag`             | `off_diag`           | 🔸 `snap_to_symmetry_svd`                                                                           |
| `maxl`                | `best_order`         | 🔸 `pointgroup_robust` — current name reads as "max-L"; this is the largest group order found      |
| `AiAiT`               | *(keep)*             | 🔸 `pointgroup_fast` — matrix-math shorthand, acceptable                                            |
| `U_basis`, `Uinv_basis`, `U_ro`, `U_or` | *(keep)* | math-convention names                                                                              |

---

## 2. Other Julia style-guide adherence (independent of the rename)

Mark each row ✅ apply / ❌ skip / ✏ modify.

### 2a. Mechanical fixes

| #   | Issue                                                                                                                                | Examples                                                                                          | Apply? |
|-----|--------------------------------------------------------------------------------------------------------------------------------------|---------------------------------------------------------------------------------------------------|--------|
| S1  | **Indent function bodies 4 spaces.** Several functions have bodies at column 1 (zero indent).| `pointGroup_simple` (lines 887–904), `pointGroup_fast` (930–962), `snapToSymmetry_SVD` (1456–86), `threeDrotation` (851–58), `aspectRatio` (1565–78) |    ✅    |
| S2  | **Use 4-space, not 5-space, indent.** `avgVecOverOps` and `snapToSymmetry_avg` use 5-space indent.                                    | lines 28–31, 50–53                                                                                |✅|
| S3  | **Spaces around binary operators.**                                                                                                  | `i<j` → `i < j`, `a≈b` → `a ≈ b`, `R*A` → `R * A` (numerous places)                              |✅|
| S4  | **Spaces after commas in arg lists.**                                                                                                | `function f(a,b,c)` → `function f(a, b, c)` — applies to ~10 function headers                     |✅|
| S5  | **Drop redundant `==true` after `findall(broadcast)`.**                                                                              | line 900: `findall([t≈I(3) for t ∈ T].==true)` → `findall(t ≈ I(3) for t in T)`                  |✅|
| S6  | **Expand single-line `if … end`.**                                                                                                   | line 898: `if debug return T end`; line 1076: `if il > length(idx) continue end`                  |don't apply|
| S7  | **Element-iteration variables shouldn't be named `i`.** `i` reads as an index; for elements prefer `v`, `x`, or a domain noun.       | `for i ∈ c` where `c::Vector{Vector{Float64}}` (lines 945–947, 953, 1053–55)                      |✅|
| S8  | **Pick one of `∈` / `in` and use consistently.** Currently mixed within the same file (e.g. `for i ∈ c` vs `for i in 1:N`).         | (rec) use `in` everywhere — easier to type and read, no semantic difference                       |    use ∈ everywhere    |
| S9  | **Promote positional-debug-flag to a keyword.**                                                                                      | `pointGroup_simple(a1, a2, a3, debug=false)` → `pointGroup_simple(a1, a2, a3; debug::Bool=false)` | 🔸 **rec ✅** — debug flag, no callers in production, keyword is idiomatic       |
| S10 | **Promote `_canonicalize_τ`'s `tol` to keyword.** (Currently second positional arg.)                                                  | line 689                                                                                          |∈|
| S11 | **Wrap long lines** (target ≤92 chars).                                                                                              | lines 1042, 1097, 1415 — all warning messages                                                     |leave as is|
| S12 | **Use `eachindex` / `axes` instead of `1:length(x)`.**                                                                                | lines 897, 956 — `for i ∈ 1:length(R)` → `for i in eachindex(R)`                                  |∈|
| S13 | **Remove `src/debugging.jl` from `src/`** (CLAUDE.md says it's not part of the module). Either delete or move to `support/` or `examples/`. | 🔸 verified: `debugging.jl`, `cubic_example.jl`, `snapExample.jl` are all orphaned (no `include` from anywhere). The two `*_example.jl` files reference `snapToSymmetry(...)` which no longer exists — they're actively broken. |🔸 **rec: delete all three.** None are loaded; two are stale. Canonical examples now live in `docs/src/`.|
| S14 | **`r in members` relies on exact `==` for `Matrix{Int}` (line 103)** — works, but the in-line comment "relies on exact == underneath" hints that the original author was uneasy. Document or guard.    | 🔸 **Not a real problem.** `Matrix{Int} == Matrix{Int}` is exact integer-tuple equality — well-defined and fast. Comment is reassurance to readers who land on the float branch (line 120) and wonder why this one's simpler. (rec) reword the comment to "(integer matrices: `==` is exact)" and move on. |Is it really a problem?|
| S15 | **`isagroup` distinctness loop uses `@view members[(k+1):end]`** but the integer-matrix variant uses `unique(members)` instead. Make the two methods structurally parallel. | 🔸 **Don't unify — leave both.** Integer `unique` is O(n) via hashing; float pairwise is O(n²) because `isapprox` has no consistent hash. Different idioms because the underlying type contracts differ. (rec) add a one-line comment in each method noting why they diverge. |which is better?|

### 2b. Style choices worth a deliberate decision (not strictly wrong now)

| #   | Issue                                                                                                                                                                                | Options                                                                                                              | Pick |
|-----|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|----------------------------------------------------------------------------------------------------------------------|------|
| C1  | **Implicit vs explicit `return`.** Currently mixed (e.g. `fractional`, `cartesian` are one-liners with no `return`; `isagroup` ends with `return true`).                              | Julia style accepts either; pick one for new/touched code. (rec) keep explicit `return` only where it improves clarity (early return, last expression is non-obvious) | 🔸 **rec: use explicit `return` for early returns (predicate functions) and at the end of any function >5 lines; drop it from one-liners and short utilities.** That matches what the file already does in most places — formalising the existing pattern means a *very* small diff.    |
| C2  | **`Crystal` constructor's `coords=` keyword as `Symbol`**. Consider a `@enum` or string literal type union — but `Symbol` is the de-facto Julia convention. (rec) leave as-is.        |your rec|      |
| C3  | **Comments inside the docstring `# Examples` blocks**: some examples are long. (rec) move heavy commentary into `# Extended help` sections so REPL `?` output stays terse.            | 🔸 plan: I'll mark candidate blocks (`pointgroup`, `spacegroup`, `Crystal`, `read_poscar`, `make_primitive`) in the diff and you decide per-function. ✅ that workflow?  |Mark these, I'll double check after they are updated|
| C4  | **`error("…")` everywhere vs typed exceptions.** Spacey throws `ErrorException` via `error(...)` for all argument validation. Julia's recommendation is `ArgumentError`/`DomainError`. | (rec) bulk-replace `error(...)` → `throw(ArgumentError(...))` for argument-validation sites only                       |✅|
| C5  | **`@warn` keyword-arg syntax** — line 1097/1415 use the rich-keyword form (`tol group_at_tol=…`); good, just verify it renders fine in tests.                                          | (rec) keep                                                                                                            |✅|
| C6  | **`pointGroup` matrix overload uses `A[:,1]`, `A[:,2]`, `A[:,3]`** rather than `eachcol(A)`.                                                                                          | (rec) `eachcol(A)...`                                                                                                |✅|

---

## 3. 🔸 Migration mechanics (revised for "all at once, release as v0.9.0")

Per your replies: change everything in one PR, no `@deprecate` aliases, ship as
v0.9.0 (not v1.0 yet).

| #   | Step                                                                                                                            | Status |
|-----|---------------------------------------------------------------------------------------------------------------------------------|--------|
| M1  | ~~Add `@deprecate` aliases~~ — **skipped per your reply.** Hard cut.                                                            | 🔸 dropped |
| M2  | Single PR with the full rename + style sweep.                                                                                   | ✅ |
| M3  | CHANGELOG entry for v0.9.0 with the full old → new mapping plus a "fix your code by" instruction.                               | ✅ |
| M4  | Update doctests in `src/` + all of `docs/src/` in lockstep. Doc build (`checkdocs = :exports`) will catch stragglers.            | ✅ |
| M5  | Tag **v0.9.0** when merged (not v1.0). v1.0 left for a later "we're really sure" release.                                       | 🔸 v0.9.0, not v1.0 |
| M6  | Final `git grep` for old names to catch leftovers.                                                                              | ✅ |
| 🔸 M7 | **Run the full test suite (`Pkg.test()`)** + **build docs** locally before commit. The latter is the only check that catches stale jldoctests. | ✅ added |
| 🔸 M8 | **Touch `Project.toml` version** to `0.9.0` in the same commit. (`2fe3397 Trigger registration of v0.8.0` shows the registration trigger pattern.) | ✅ added |

---

## 4. Things I deliberately did **not** suggest renaming

- Mathematical single-letters: `A`, `R`, `τ`, `u`, `v`, `w`, `T` — these are the standard symbols in crystallography and linear algebra; renaming `A` to `lattice_matrix` would harm readability for the target audience.
- Unicode operators (`∈`, `×`, `⋅`, `∛`) and unicode field name `τ` — idiomatic in mathematical Julia; the existing usage is fine.
- `Crystal` / `SpacegroupOp` type names — already conform to UpperCamelCase.

---

## 🔸 5. Risk / footguns of doing it all at once

Short list — these aren't reasons not to do it, just things I'll watch for.

| #  | Risk                                                                                                                          | Mitigation                                                                                  |
|----|-------------------------------------------------------------------------------------------------------------------------------|---------------------------------------------------------------------------------------------|
| R1 | **Downstream breakage**: any package depending on Spacey 0.8 will fail to load after upgrade with no deprecation hint.        | Loud CHANGELOG (M3) + version bump signals the break. Julia's `Pkg` resolver respects SemVer; 0.8 → 0.9 is a breaking bump and downstreams must opt in. |
| R2 | **Doctest drift**: `docs/src/**` and inline jldoctests reference old names in prose *and* code blocks. Code blocks fail loudly; prose stays stale silently. | Grep prose for old names (M6) — same regex catches both.                                    |
| R3 | **`SpacegroupOp` field name `R`** is referenced positionally in tests and may also appear in user docstrings as `op.R`. The struct field doesn't get renamed — but I should double-check no one is renaming it accidentally. | Field stays. Confirmed by Read of `src/Spacey.jl:674–747`.                                  |
| R4 | **AFLOW test generator** (`tools/generate_aflow_tests.jl`) likely emits the old function names; auto-regen will overwrite the test file if the generator isn't updated first. | Update the generator *before* re-running it; if the existing test file is hand-edited (no regen needed), just rename in place. |
| R5 | **`benchmark/benchmarks.jl`** references old names → PkgBenchmark `judge` between `main` and `HEAD` won't run cleanly across the rename. | Either rename benchmarks in the same PR, or accept that the v0.8 ↔ v0.9 benchmark comparison won't work and skip it for this release. |
