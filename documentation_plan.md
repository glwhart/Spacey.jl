# Documentation Plan for Spacey.jl

A planning document for building user-facing documentation for Spacey.jl. Follows the directions in `spacegroup_plan.md` §7. Reads of philosophy resources, the existing source docstrings, and our own theory documents (`research.md`, `plan.md`, `designDiscussions.md`) inform the proposal below.

This is a planning document — not the docs themselves. Decisions and trade-offs are recorded so the actual writing has a clear target.

---

## 1. Documentation philosophy: what we're adopting

### 1.1 Diátaxis as the spine

The Diátaxis framework (Daniele Procida) splits documentation into four mutually exclusive types based on the user's situation along two axes:

|              | **Skill acquisition** (learning) | **Skill application** (working) |
|---           |---                               |---                              |
| **Action**   | **Tutorials**                    | **How-to guides**               |
| **Cognition**| **Explanation**                  | **Reference**                   |

The four types serve different needs and require different writing styles. Mixing them is the most common source of confusing documentation.

The "compass" decision rule: when adding a piece of documentation, ask *what is the reader doing?* (action vs cognition) and *what state are they in?* (learning vs applying). The answer points at exactly one quadrant.

Adopting Diátaxis means the Spacey docs are organised into four top-level sections that match those quadrants. We will not invent a fifth category, nor split content randomly into "Getting Started", "API", "Cookbook" etc. as is common but anti-pattern.

### 1.2 Supporting principles from Write the Docs

Diátaxis tells you the *structure*; Write the Docs adds *quality* practices that complement it:

- **Skimmability.** Descriptive headings, lead with the concept, link aggressively. Readers rarely read sequentially.
- **Concrete examples.** Every page that allows it should have a runnable code block.
- **Source proximity.** Docstrings live next to the code they describe; written prose lives in `docs/src/`. Single source of truth for each.
- **Currency.** Out-of-date documentation is worse than missing — it actively misleads. Doctests + CI catch the drift.
- **Strategic repetition.** A small amount of repetition (e.g. example signatures in both reference and how-to) is better than aggressive de-duplication that forces readers to chase links.

### 1.3 Other approaches considered, and why we're not using them

- **DITA** (Darwin Information Typing Architecture, OASIS standard) — XML-schema based, designed for large enterprise documentation pipelines. Far too heavy for a Julia research package; we'd spend more time fighting the schema than writing content.
- **Information Mapping** — a commercial methodology around chunking and labelling. Useful for procedural manuals; less applicable to a scientific library where explanation/theory is a first-class concern.
- **The Good Docs Project** templates — high-quality templates per content type, but conceptually nested under Diátaxis (their templates explicitly map to the four quadrants). Worth borrowing structure from for individual pages, but not a separate framework.
- **Read-the-Docs default ("Getting Started / Tutorial / API / FAQ")** — the conventional layout. Comfortable, but it conflates how-to with tutorial and buries explanation. Diátaxis is a strict superset and identifies exactly where this default fails.
- **Documenter.jl's auto-`@autodocs`** (used by MinkowskiReduction.jl's docs today, see §3.4 below) — fine for a single-page reference but provides nothing for tutorials, how-to, or explanation. We will use `@autodocs` for the reference quadrant only.

The decision: **Diátaxis structure + Write-the-Docs quality + Documenter.jl tooling**. None of the alternatives offer something Diátaxis doesn't already.

---

## 2. Sources we already have

Most of the Spacey documentation content already exists in scattered form. Mapping each source to its proper Diátaxis quadrant tells us what to lift verbatim, what to adapt, and what to write fresh.

| Source                       | Best target quadrant            | What to do                                                                                        |
|---                           |---                              |---                                                                                                |
| Source code docstrings       | Reference                       | Lift via `@autodocs` — already-live and tested via doctests where present                          |
| `research.md`                | Explanation                     | Selectively lift theory: Minkowski-reduction completeness theorem, tolerance failure modes, near-boundary discussion. The research-narrative tone needs trimming |
| `plan.md`                    | Mostly internal, some Explanation | The §1 robustness-strategy overview is good Explanation source. The §2 / §3 backlog items are internal |
| `designDiscussions.md`       | Explanation (in part)           | The pos_tol-α, SpacegroupOp-rationale, τ-canonicalisation, and crystal_system Layer-1-vs-2 discussions all illuminate "why" — lift trimmed versions |
| `spacegroup_plan.md`         | Internal                        | Stays internal — planning history, not user-facing                                                 |
| `phase2_plan.md`             | Internal                        | Frozen design-review record, not user-facing                                                       |
| `aiCodeExplanation.md`       | Internal                        | Julia-syntax notes for the project owner; not for end users                                        |
| `CLAUDE.md`                  | Internal                        | Explicitly for AI-assisted development; out of scope for user docs                                 |
| `test/runtests.jl`           | Tutorial / How-to source        | Existing tests are a cookbook of "how to use this on a real crystal" — mine for tutorial examples  |
| `test/aflow_structures_part*.jl` | Reference (linked, not lifted) | Mention as a validation corpus; don't reproduce 1095 prototypes in user docs                       |
| `test/nearMissBoundary*.jl`  | Explanation                     | The diagnostic heatmaps are excellent illustrations of the over-promotion failure mode             |

---

## 3. Proposed documentation structure

A `docs/src/` tree organised by the four Diátaxis quadrants, with an entry-point `index.md` that orients new readers.

```
docs/src/
├── index.md                              ← landing: 1-minute orientation + 4 paths
├── tutorials/
│   ├── index.md                          ← 1-paragraph intro, links to:
│   ├── 01-first-pointgroup.md            ← "Find the point group of a cubic lattice"
│   └── 02-first-spacegroup.md            ← "Find the space group of NaCl"
├── how-to/
│   ├── index.md                          ← intro + table of contents
│   ├── construct-a-crystal.md            ← `Crystal` constructor recipes (frac & cart input, 3-vector vs matrix)
│   ├── find-pointgroup.md                ← Choose between simple / fast / robust; tolerances
│   ├── find-spacegroup.md                ← spacegroup() invocation patterns
│   ├── classify-bravais.md               ← crystal_system + (future) full Bravais
│   ├── handle-noisy-data.md              ← Pos_tol + lattice_tol selection from real-data noise levels
│   ├── detect-tolerance-dependence.md    ← verify_stable usage on both pointGroup_robust and spacegroup
│   ├── compose-and-apply-ops.md          ← Working with `SpacegroupOp`: composition, inverse, apply, toCartesian
│   └── snap-to-symmetry.md               ← snapToSymmetry_SVD / _avg workflows
├── reference/
│   ├── index.md                          ← `@autodocs` of public API in canonical order
│   ├── crystals.md                       ← Crystal type, fractional/cartesian, default_pos_tol, crystal_system
│   ├── point-groups.md                   ← pointGroup_robust, _fast, _simple, pointGroup
│   ├── space-groups.md                   ← spacegroup, isSpacegroupOp, SpacegroupOp methods, toCartesian
│   ├── snap.md                           ← snapToSymmetry_SVD, snapToSymmetry_avg
│   └── helpers.md                        ← isagroup, threeDrotation, aspectRatio
└── explanation/
    ├── index.md                          ← intro: "why Spacey decides what it decides"
    ├── why-minkowski.md                  ← Theorem: 27-neighbour search is complete on Mink-reduced bases
    ├── tolerances.md                     ← The two-tolerance model (lattice vs position); near-boundary failure
    ├── over-promotion.md                 ← The classic failure mode; how `verify_stable` flags it; near-miss heatmaps
    ├── canonicalising-tau.md             ← Why τ snaps to small rationals; the round-vs-snap episode
    ├── crystal-system-vs-bravais.md      ← Why we provide system (Layer 1) but not centering (Layer 2 deferred)
    ├── algorithm-overview.md             ← End-to-end pipeline: minkReduce → 27-neighbour → filter → close → spacegroup
    └── validation-strategy.md            ← Hand-built tests + AFLOW corpus; what each catches
```

12 user-facing pages, plus four section indexes and the landing page = 17 markdown files. Each page targets a single Diátaxis quadrant; cross-links between quadrants are heavy but content is not duplicated.

### 3.1 Landing page — what `index.md` does

Three jobs:
1. **Tagline + one-paragraph what-is-Spacey:** a sentence the user can read in 10 seconds and decide whether they're in the right place.
2. **Four-path navigation** mirroring the Diátaxis quadrants, with one-sentence descriptions:
   - "I'm new — *follow a tutorial*"
   - "I have a specific task — *find a how-to guide*"
   - "I need to look something up — *consult the reference*"
   - "I want to understand how Spacey decides — *read an explanation*"
3. **One install-and-test snippet** so the reader can verify the package loads in 30 seconds.

### 3.2 Tutorials (2 of them, narrowly scoped)

A tutorial in Diátaxis is a *guided learning experience*, not a tour. It must produce a concrete, visible result and never explain. Two are enough for a small library:

- **`01-first-pointgroup.md`** — "Find the point group of a simple cubic lattice." Walks through importing Spacey, defining `A = I`, calling `pointGroup_robust`, getting 48 ops, trying it on a tetragonal lattice and getting 16. ~5 minutes.
- **`02-first-spacegroup.md`** — "Find the space group of NaCl." Constructs a `Crystal` with two atoms, runs `spacegroup`, gets 192 ops. Touches identity-at-index-1, shows applying an op to an atom position. ~10 minutes.

Each tutorial:
- Shows the input *and* the expected output verbatim, so readers know they're on track.
- Avoids vocabulary beyond what's needed for the next step.
- Does NOT explain *why* the count is 48 or 192 — that's an explanation page, linked from the wrap-up.

### 3.3 How-to guides (8 of them, problem-focused)

A how-to guide assumes the user knows what they want; it shows the steps to get there. The 8 above cover the canonical questions a working scientist would have:

1. **Construct a crystal** — `Crystal` constructor patterns; required `coords` kwarg; matrix vs three-vector input; integer / rational / Float input.
2. **Find a point group** — when to use `_simple` (validation), `_fast` (clean inputs), `_robust` (real data); tolerance choice.
3. **Find a space group** — `spacegroup(c)` workflow, choosing `pos_tol` and `lattice_tol`.
4. **Classify the Bravais system** — `crystal_system(A)` usage and what its 7 outputs mean.
5. **Handle noisy real-world data** — translate experimental noise levels into appropriate tolerances; the `default_pos_tol` formula; what to do when noise exceeds typical bond lengths.
6. **Detect tolerance-dependent answers** — `verify_stable=true` on both finders; when the warning fires; what to do about it.
7. **Compose and apply ops** — `SpacegroupOp` operator overloads; getting a Cartesian form via `toCartesian`.
8. **Snap to symmetry** — `snapToSymmetry_SVD` vs `snapToSymmetry_avg` trade-offs; volume preservation.

Each how-to is ≤ 1 page. Pure recipe form: numbered steps, code, expected output. No explanation; cross-link to the explanation pages where relevant.

### 3.4 Reference (auto + curated)

The reference quadrant is austere and complete. For Spacey:

- **Each reference page** uses `@docs` blocks to inject the docstrings of a coherent function group (Crystal-related on one page, point-group functions on another, etc.). This keeps source-of-truth in the source code.
- **Augmented with `@index` blocks** so each section has a typeable function index.
- **No tutorial or explanation content** in reference pages — strictly signatures, parameters, returns, doctests.

This is the section MinkowskiReduction.jl currently does as a single `@autodocs` page (see `~/.julia/dev/MinkowskiReduction/docs/src/index.md`). For Spacey we'll group by topic across multiple pages — the API is larger and topical grouping is more navigable.

**Pre-condition for the reference quadrant to land cleanly:** every public symbol needs a docstring in the source. Audit shows current coverage is high but inconsistent (some private helpers have docstrings, some public symbols have minimal ones). Pre-flight: docstring audit + cleanup before docs publication.

### 3.5 Explanation (7 pages, each illuminating a "why")

These are the most valuable pages for a research library — they're the difference between "a user can call the function" and "a user understands when the function will mislead them". Sources are pre-existing prose in our internal documents:

| Explanation page              | Lifts from                          | What it covers                                                                                                |
|---                            |---                                  |---                                                                                                            |
| Why Minkowski reduction       | `research.md` §2.1, Hart-Jorgensen-Morgan-Forcade 2019 §A.1 | The 27-neighbour completeness theorem; intuition for *why* a Minkowski-reduced basis is sufficient            |
| Tolerances                    | `research.md` §2.1, `designDiscussions.md` pos_tol section | The two-tolerance model; physical meaning; how to translate experimental noise to tolerance values |
| Over-promotion                | `research.md` near-boundary, `test/nearMissBoundary*.jl` heatmaps | The classic failure mode; how `verify_stable` flags it; the diagonal `tol ≈ ε` crossover |
| Canonicalising τ              | `designDiscussions.md` τ section    | Why τ snaps to rationals with q ≤ 12; the round-vs-snap episode and what it teaches about FP rationals       |
| Crystal system vs full Bravais| `designDiscussions.md` Layer-1 section | Why we provide the 7-class identification cheaply and defer the 14-class one                                 |
| Algorithm overview            | `plan.md` §1, `CLAUDE.md` Core Algorithm | The end-to-end pipeline at a single conceptual level                                                          |
| Validation strategy           | `plan.md` §3.13, this session's AFLOW work | Hand-built canonical tests, the AFLOW corpus, the op-count + crystal_system dual invariant, the deviation classes |

Each explanation page:
- Discusses, doesn't instruct.
- Acknowledges alternatives where they exist (e.g. "Niggli vs Minkowski reduction"; "Layer 1 vs Layer 2 Bravais").
- Includes references to literature where relevant (the AFLOW papers, Mehl, Hicks, Curtarolo et al.; Minkowski's original; Niggli; the FINDSYM / Spglib references).
- Links cross-quadrant: forward to relevant how-tos, back to relevant reference.

---

## 4. Theory pieces from `research.md` and `plan.md` worth surfacing

Going through both files line-by-line, the user-facing-relevant pieces:

### From `research.md`
- **§2.1 Minkowski reduction** — the 12-condition definition, the 27-neighbour theorem, the "epsilon-violation never breaks the search" result. **Surface in:** Explanation / Why Minkowski reduction.
- **§2.2 Niggli, §2.3 LLL** — alternatives compared. **Surface in:** Explanation / Why Minkowski reduction (closing paragraph).
- **§3 (FINDSYM, Spglib, AFLOW-SYM)** — comparison with related tools. **Surface in:** Explanation / Algorithm overview (closing paragraph: "Spacey vs adjacent tools").
- **The over-promotion section** — the canonical failure mode. **Surface in:** Explanation / Over-promotion.
- **The metric-tensor section** — useful background for users who'd want to extend toward Bravais Layer 2. **Surface in:** Explanation / Crystal system vs full Bravais (forward-looking note).

### From `plan.md`
- **§1 (Why the code is already strong)** — the high-level strategy. **Surface in:** Explanation / Algorithm overview.
- **§3.13 (AFLOW validation)** — already documented as landed. **Surface in:** Explanation / Validation strategy.

Most of `plan.md` §2 (the backlog of tolerance / consistency improvements) is internal-developer-facing and stays internal.

### From `designDiscussions.md`
- **pos_tol α discussion** — "why is it 1% of (V/N)^(1/3)?" is exactly the kind of question a serious user asks. **Surface in:** Explanation / Tolerances.
- **SpacegroupOp rationale** — surfaces only as much as needed in the SpacegroupOp how-to.
- **τ canonicalisation episode** — surfaces in Explanation / Canonicalising τ.
- **Crystal system Layer-1-vs-Layer-2** — surfaces in Explanation / Crystal system vs full Bravais.

Most of these go in trimmed form — the discussion-document narrative ("we considered X, then Y, then settled on Z") is internal. The user-facing version states the conclusion plus the reasoning.

---

## 5. Technical setup

### 5.1 Documenter.jl configuration

Replace the current minimal `docs/make.jl` (8 lines, no deployment) with a full setup:

```julia
using Documenter
using Spacey

DocMeta.setdocmeta!(Spacey, :DocTestSetup, :(using Spacey, LinearAlgebra); recursive=true)

makedocs(
    sitename = "Spacey.jl",
    authors  = "Gus Hart and contributors",
    repo     = "https://github.com/glwhart/Spacey.jl/blob/{commit}{path}#{line}",
    format   = Documenter.HTML(
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical  = "https://glwhart.github.io/Spacey.jl",
        assets     = String[],
    ),
    modules  = [Spacey],
    pages = [
        "Home" => "index.md",
        "Tutorials" => [
            "tutorials/index.md",
            "tutorials/01-first-pointgroup.md",
            "tutorials/02-first-spacegroup.md",
        ],
        "How-to guides" => [
            "how-to/index.md",
            "how-to/construct-a-crystal.md",
            "how-to/find-pointgroup.md",
            "how-to/find-spacegroup.md",
            "how-to/classify-bravais.md",
            "how-to/handle-noisy-data.md",
            "how-to/detect-tolerance-dependence.md",
            "how-to/compose-and-apply-ops.md",
            "how-to/snap-to-symmetry.md",
        ],
        "Reference" => [
            "reference/index.md",
            "reference/crystals.md",
            "reference/point-groups.md",
            "reference/space-groups.md",
            "reference/snap.md",
            "reference/helpers.md",
        ],
        "Explanation" => [
            "explanation/index.md",
            "explanation/why-minkowski.md",
            "explanation/tolerances.md",
            "explanation/over-promotion.md",
            "explanation/canonicalising-tau.md",
            "explanation/crystal-system-vs-bravais.md",
            "explanation/algorithm-overview.md",
            "explanation/validation-strategy.md",
        ],
    ],
)

deploydocs(
    repo      = "github.com/glwhart/Spacey.jl",
    devbranch = "main",
    devurl    = "dev",
    target    = "build",
    branch    = "gh-pages",
    versions  = ["stable" => "v^", "v#.#"],
)
```

### 5.2 Doctests

Every example in tutorials, how-to, and (especially) reference should be a `jldoctest` block. CI runs `Pkg.test()` followed by a `doctest = true` makedocs and fails on stale examples. This is the primary defence against documentation rot.

### 5.3 GitHub Actions

A new workflow `.github/workflows/Documentation.yml`:
- Triggers on push to `main` and on tags.
- Builds docs with `julia --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate(); include("docs/make.jl")'`.
- Deploys to gh-pages via Documenter's built-in deploy.

### 5.4 Pre-launch source-side prerequisites

- **Docstring audit.** Every public symbol must have a docstring with at least: signature line, one-sentence summary, parameter list, return description, and one example (preferably as a `jldoctest`). Currently incomplete on a few low-level helpers.
- **Doctest passing.** Existing doctests in the source need to be updated where they reflect old API (the old `pointGroup_robust` doctest example uses syntax that wouldn't run today).
- **Public/private boundary.** Some currently-`export`ed symbols are arguably internal (`isagroup`, `aspectRatio`?). Either document them or stop exporting. Decide before the reference quadrant is published.

---

## 6. Phasing

Don't write everything at once. Phased rollout:

### Phase D1 — Infrastructure (~½ day)
- New `docs/make.jl` with the full page tree (mostly empty pages).
- GitHub Actions workflow for docs build + deploy.
- Confirms the full skeleton compiles and deploys to gh-pages before content writing.

### Phase D2 — Reference (~1 day)
- Cleanest quadrant to do first: lifts directly from existing source docstrings.
- Audit + improve docstrings as part of this phase.
- Add `@docs` / `@autodocs` blocks per reference page.
- Doctest fixes as needed.

### Phase D3 — How-to guides (~1–2 days)
- Each how-to is a small, focused, runnable recipe.
- The 8 pages can be written largely independently; mining `test/runtests.jl` provides realistic examples for several.

### Phase D4 — Tutorials (~½ day)
- Two pages, narrowly scoped.
- Probably the highest-impact-per-line documentation in the project, since first-time-users often bounce if a tutorial is missing or unclear.

### Phase D5 — Explanation (~1–2 days)
- Most expensive in writing time but most valuable for serious users.
- Heavy lift from existing prose in `research.md` / `designDiscussions.md`, but careful trimming of internal-narrative tone.

### Phase D6 — Polish + cross-linking (~½ day)
- Add cross-references throughout (every reference page has "see also" pointing to the relevant explanation; every how-to has links to relevant reference).
- Landing-page polish.
- Search-engine sitemap, badge in the README pointing at the docs URL.

**Total: ~5 days of concentrated work.** Phase D1 should land first as a standalone commit so the deployment is verified before content investment.

---

## 7. Open questions

Things to resolve before Phase D1 begins:

1. **Public API boundary.** ~~Does `isagroup` belong in the public reference, or is it an internal helper? Same question for `aspectRatio`, `threeDrotation`.~~ **Resolved (2026-04-22):** `isagroup` stays exported (some users may find it useful as a sanity check on a custom op set). `aspectRatio` and `threeDrotation` are internal helpers — un-exported, docstrings re-headed as `Spacey.aspectRatio(...)` / `Spacey.threeDrotation(...)`, callers in `test/runtests.jl` updated, CLAUDE.md API table trimmed. Reference quadrant will document only the public set; the two internal helpers can still appear under `## Internals` on `reference/helpers.md` for those who want to reach for them via the qualified name.

2. **Doctest output style.** ~~Some existing doctest examples used floating-point output that's machine-dependent.~~ **Resolved (2026-04-22):** all 16 jldoctests now use exact, machine-independent output (integer matrices, structurally-zero comparisons, or values that are exactly representable in Float64). No `output =` blocks needed. Where exactness would have made an example "miss the point" (e.g. illustrating that a near-orthogonal basis has aspect ratio close to but not exactly 1), the example was rewritten with an obviously-exact input instead.

3. **Whether to publish the AFLOW corpus stats** — explicitly the 1095 / 56-deviation breakdown — in user-facing docs. Probably yes (in Explanation / Validation), with the disclaimer that the broken-pinned set is a known-quirk catalog, not a bug list.

4. **Who's the target reader?** ~~First-year-grad student in computational materials? Methodologist who'd be interested in the algorithmic pieces?~~ **Resolved (2026-04-22):** target reader is "scientist who knows enough crystallography to read Mehl et al." — i.e. someone comfortable with conventional vs primitive cells, point/space groups, the 14 Bravais types, and Wyckoff notation, but who doesn't necessarily think about finite-precision algorithm design daily. Tutorials should not teach crystallography; how-to guides assume working vocabulary; explanation pages can pull in algorithmic depth without apologising.

5. **Versioning policy.** ~~Decide whether we maintain `dev` + `stable` + per-tag versions or simpler `stable`-only.~~ **Resolved (2026-04-22):** `stable`-only. Documenter's `versions = ["stable" => "v^"]` setting publishes a single tracked version that follows the highest tag. Simpler to reason about and matches the size of the project.
---

## 8. Out of scope for this doc plan

- Interactive in-browser examples (e.g. via JuliaHub's Pluto or similar). Beyond the typical Documenter-on-GitHub-Pages baseline.
- Full localisation. English-only.
- Reference cross-linking with MinkowskiReduction.jl's docs. Could add `inventory` or external-references later; not in initial scope.
- Migration of MinkowskiReduction.jl's docs to a similar structure. See §8.2 below for what that would entail; ultimately it's a separate decision for that package.

### 8.1 Visual diagrams — candidates for inclusion (was: out of scope)

The user's note on §8 asked for a list of candidate diagrams. The recommendation now is: pick a small subset (3–5) for the initial publication, defer the rest. Each candidate is rated for **value** (how much it clarifies for a reader) and **cost** (effort to produce + maintain).

**Diagrams sourced from existing publications (low cost — just embed):**

1. **Hart-Jorgensen-Morgan-Forcade 2019, Figure A1** — 2D illustration that every point on the boundary of `U_R` (union of basis cells around origin) has an interior cousin closer to the origin when the basis is Minkowski-reduced. **Value: high** — directly motivates *why* a bounded search space suffices. **Cost: low** — paper is open access (CC-BY 3.0), so a clean reproduction or a re-drawn close analogue is straightforward. **Place in docs:** `explanation/why-minkowski.md`.

2. **Hart-Jorgensen-Morgan-Forcade 2019, Figure A2** — the 3D version of A1: planes through the origin, showing the proof structure. **Value: medium-high** — useful for readers who want to follow the proof step into 3D. **Cost: low** — same open-access source. **Place:** same page as A1, as a "see also in 3D" supplement.

3. **Hart-Jorgensen-Morgan-Forcade 2019, Figure 8** — the 2D "closest cousin" guarantee: Brillouin zone (blue) is contained in the union of 4 basis cells (red) when the basis is Minkowski-reduced; if not, regions can spill into a green cell outside the union. **Value: high** — the single most intuitive picture of what Minkowski reduction *buys* you. **Cost: low**. **Place:** `explanation/why-minkowski.md`, leading figure.

   **(Note for MinkowskiReduction.jl):** Figures 8, A1, and A2 would be at least as valuable in MinkowskiReduction.jl's docs as in Spacey's. Currently that library's docs are a single `@autodocs` page with no explanation of the *purpose* of Minkowski reduction. A "Why Minkowski reduction" page there, anchored on Figure 8 plus A1, would close that gap. This is a recommendation for a sibling-package PR, not part of the Spacey docs build.

**Diagrams to generate from Julia tooling (medium cost):**

4. **27-cell {-1, 0, 1}³ search space.** A 3D wireframe of 27 unit cells around a central origin cell, with the candidate vectors (the {-1,0,1} integer combinations of basis vectors) drawn as arrows from the origin. Highlight the at-most-26 non-zero candidates Spacey actually checks. **Value: high** — the central algorithmic claim made concrete. **Cost: medium** — ~30 lines of GLMakie. **Place:** `explanation/algorithm-overview.md`.

5. **Pipeline flowchart.** Boxes-and-arrows: `input A` → `minkReduce` → `27-candidate generation` → `norm filter` → `volume check` → `integer-matrix check` → `group-closure search` → `(LG, G)`. Annotated with the tolerance that gates each step. **Value: high** — gives a reader who's read the docstrings a single mental model. **Cost: low** — Mermaid (Documenter has built-in support). **Place:** `explanation/algorithm-overview.md`.

6. **Crystal-system decision tree.** Holohedry order → crystal system: 2 → triclinic, 4 → monoclinic, 8 → orthorhombic, 12 → trigonal, 16 → tetragonal, 24 → hexagonal, 48 → cubic. **Value: medium** — a small concept, but the table form may be enough; a Mermaid flowchart adds clarity for first-time readers. **Cost: low** — Mermaid. **Place:** `explanation/crystal-system-vs-bravais.md` and/or `how-to/classify-bravais.md`.

7. **Over-promotion heatmap.** The `(tol, ε)` 2D scan (already produced by `test/nearMissBoundary*.jl`): x-axis tolerance, y-axis perturbation magnitude, cell colour = group order found. Diagonal failure region is visually striking. **Value: high** — turns a subtle numerical claim into a picture. **Cost: medium** — figures already exist as test artefacts; need cleanup pass to publication quality. **Place:** `explanation/over-promotion.md`.

**Diagrams that would be nice but probably not worth the cost for initial publication:**

8. **All 14 Bravais lattices, 3D rendered.** Useful pedagogy, but redundant with Wikipedia and standard textbooks. Skip in favour of linking out (e.g. to Wikipedia or DoITPoMS) until/unless we ship Layer 2 Bravais classification.

9. **Brillouin-zone visualisations.** [Brillouin.jl](https://thchr.github.io/Brillouin.jl/stable/) (by Thomas Christensen) constructs and plots Wigner-Seitz cells for arbitrary lattices using PlotlyJS or GLMakie. Useful only if Spacey grows a reciprocal-space story (which it currently doesn't). Defer.

10. **SpacegroupOp action visualisation.** Animated rotation/translation of an atom under successive `(R, τ)` ops. **Value: medium** for the right reader, **cost: high** (requires animation, substantial Makie work). Defer.

11. **Crystal-with-atoms 3D rendering.** [LatPhysPlottingMakie.jl](https://github.com/ffreyer/LatPhysPlottingMakie.jl) provides Makie plot recipes for lattices with sites and bonds; useful for tutorial NaCl/diamond examples. **Cost: medium** (requires adding a docs-only dep). Reasonable to fold in during Phase D4 (Tutorials) if a reader-facing picture would clarify the `Crystal` constructor.

**Tooling notes:**

- **Mermaid** (flowcharts, decision trees) — supported natively by Documenter; lowest-friction option.
- **GLMakie/CairoMakie** (3D and publication-quality 2D plots) — Julia ecosystem standard. CairoMakie produces SVG/PDF suitable for static docs; GLMakie for interactive screenshots.
- **PGFPlotsX / TikzPictures** — high-quality, slow; overkill unless we want LaTeX-grade typography.
- **Existing figures (Figs 8, A1, A2 from kpointFolding paper)** — re-use under CC-BY 3.0 with attribution. Cleaner than re-drawing if we can get the originals from the corresponding author (the project owner).

**Recommended initial set for publication:** items 1, 3, 4, 5, 7. This covers the four highest-value pictures (the two boundary-cousin proofs, the 27-cell search space, the algorithm pipeline, the over-promotion heatmap) with a mix of low and medium cost. Items 2, 6 are low-cost additions if time allows.

### 8.2 What "migrating MinkowskiReduction.jl docs to a similar structure" would mean

Mentioned in passing earlier and worth being explicit about. MinkowskiReduction.jl currently has a single docs page (`docs/src/index.md`) that uses one `@autodocs` block to list all exported symbols. That gives reference coverage but provides:

- no tutorial (a new user has to read the source to learn how to call `minkReduce`),
- no how-to guide (no recipe for "I have a noisy basis and want to verify it's reduced"),
- no explanation (no answer to "why Minkowski reduction over Niggli or LLL?", "what does `orthogonalityDefect` mean physically?", or "when does the 15-iteration limit get hit?").

Migrating to a Diátaxis structure would mean splitting that one page into the same four-quadrant tree we're proposing here for Spacey, with content drawn from:
- the existing docstrings (already lift directly to reference),
- the 2019 Hart-Jorgensen-Morgan-Forcade paper §A.1 + Figures 8/A1/A2 (explanation page on "why Minkowski reduction"),
- a couple of small worked examples for the tutorial (one on `minkReduce`, one on verifying with `isMinkReduced`).

It's a smaller package than Spacey, so the Diátaxis tree would be lighter (probably ~6–8 pages vs 17). The argument for doing it is that Spacey's "why Minkowski reduction" explanation page logically belongs in MinkowskiReduction.jl's docs and Spacey would link to it. The argument against is that it's a separate package owned by the same author, so it's effectively a parallel project rather than a dependency of this one.

**Recommendation:** treat the MinkowskiReduction.jl migration as a separate, lower-priority follow-up. Spacey ships its own self-contained explanation pages first; if/when MinkowskiReduction.jl gets the same treatment, Spacey's docs swap their inline explanation for a cross-link.

---

## 9. References / further reading

- Daniele Procida, *[Diátaxis](https://diataxis.fr/)*. Primary framework, full doc-philosophy site.
- *[Diátaxis compass](https://diataxis.fr/compass/)*, action-vs-cognition × acquisition-vs-application decision rule.
- *[Write the Docs Documentation Principles](https://www.writethedocs.org/guide/writing/docs-principles/)*. Quality-of-writing guidance complementary to Diátaxis structure.
- *[Documenter.jl Guide](https://documenter.juliadocs.org/stable/man/guide/)*. Julia ecosystem standard for package docs.
- *[Julia Manual: Documentation](https://docs.julialang.org/en/v1/manual/documentation/)*. Docstring conventions and the `@doc` macro.
- *[MinkowskiReduction.jl docs](https://glwhart.github.io/MinkowskiReduction.jl/)*. Sibling library; current docs are minimal `@autodocs` only — a baseline to exceed rather than match.
- Mehl, Hicks, Toher, Levy, Hanson, Hart, Curtarolo (2017), *The AFLOW Library of Crystallographic Prototypes: Part 1.* Comput. Mater. Sci. **136**, S1–S828.
- Hicks et al. (2019), Part 2. Comput. Mater. Sci. **161**, S1–S1011.
- Hicks et al. (2021), Part 3. Comput. Mater. Sci. **199**, 110450.

---

## 10. Status

**Plan only.** No documentation files have been written or modified by this proposal. The current `docs/src/index.md` (two lines) and `docs/make.jl` (8 lines) are unchanged.

Next step pending user approval: begin Phase D1 (infrastructure).
