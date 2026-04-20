# Session Summary

This file summarizes all prompts and responses from the Claude Code session on this project.

---

## Prompt 1: Initial Analysis and Planning

**User:** Asked Claude to carefully read all code in Spacey.jl, deeply understand why it is robust, then write two files:
- `plan.md`: Suggestions for improving robustness and new unit tests.
- `research.md`: Survey of other algorithms and codes for point group / space group finding.

Context given: The code calculates the point group of a lattice using Minkowski reduction for numerical robustness. Future work will add space group functionality. The immediate goal is to make the point group finder "bulletproof" against finite precision errors, edge cases, and high aspect ratio unit cells.

**Response:** Read all source files (`src/Spacey.jl`, `test/runtests.jl`, examples, debugging scripts). Wrote both documents:
- `plan.md` with 11 code change suggestions and 12 new test suites.
- `research.md` surveying spglib, AFLOW-SYM, FINDSYM, SOFI, Symmetrizer, and various numerical techniques.

---

## Prompt 2: Metric Tensor Formulation

**User:** Asked for more detail about the metric tensor formulation. Is avoiding `inv(A)` ever important for a 3x3 matrix? Do I need to worry about conditioning?

**Response:** No. For 3x3 matrices, `inv(A)` conditioning is negligible—even at aspect ratio 500, the error is ~10⁻¹³, far below any useful tolerance. Furthermore, the metric tensor formulation (`G' ≈ G` check) is actually *worse* than Spacey.jl's current `T ≈ I` check for high aspect ratios, because `isapprox` with `rtol` gets dominated by the large b² entries of G, becoming blind to the short-axis differences. The current T ≈ I check is uniformly scaled (all entries ~1) and better suited. The metric tensor formulation is useful for lattice *identification* (Bravais type lookup) and rational arithmetic, but not for the symmetry *search* step.

---

## Prompt 3: Aspect Ratio Warning

**User:** Suggested adding a `@warn` for high aspect ratio inputs, since large aspect ratios are unlikely in practice and a warning is simpler than making the code robust against extreme cases.

**Response:** Agreed. Proposed a warning at aspect ratio > 100 (5× safety margin before the ~500 failure zone). Added it to `pointGroup_robust` after volume normalization.

---

## Prompt 4: Don't mention BigFloats in the warning

**User:** Rejected the initial edit that mentioned BigFloats as an option. Just wants the warning without suggestions.

**Response:** Simplified the warning to just report the ratio and note unreliability above ~500. Edit accepted and applied.

---

## Prompt 5: Session Recap

**User:** Asked for a recap of the whole session.

**Response:** Provided structured summary of all work done: initial analysis, plan.md, research.md, the metric tensor discussion conclusion, and the `@warn` code change.

---

## Prompt 6: Continuing on Another Computer

**User:** Asked how to continue the session on another computer.

**Response:** Explained that session files are stored in `~/.claude/projects/` and must be manually copied to the same path on the other machine, then `claude --continue` from the project directory.

---

## Prompt 7: Clarification (where to type claude --continue)

**User:** Asked whether to type `claude --continue` in the project folder or in `.claude/projects`.

**Response:** In the project folder. Claude associates sessions with the working directory.

---

## Prompt 8: "3"

**User:** Sent "3" (unclear reference).

**Response:** Asked for clarification.

---

## Prompt 9: Where to copy session files

**User:** Confused about where to copy session files on the other machine.

**Response:** Copy to the same `~/.claude/projects/` path (matching subfolder structure) on the other machine. Then cd into the project directory and run `claude --continue`.

---

## Prompt 10: Google Drive shared projects

**User:** Asked if running on Google Drive means the session file is already there, so just `claude --continue` is enough.

**Response:** Yes, if both machines share the same Google Drive path for `~/.claude/projects/`, just cd into the project and `claude --continue`.

---

## Prompt 11: "where am i"

**User:** Asked current directory.

**Response:** `/Users/hart/home/codes/Spacey.jl`

---

## Prompt 12: Open VS Code

**User:** Asked to open VS Code in the project folder.

**Response:** Ran `code /Users/hart/home/codes/Spacey.jl`.

---

## Prompt 13: Commit new files

**User:** Asked to commit the new files so they're available on other computers.

**Response:** Committed `src/Spacey.jl` (aspect ratio warning), `plan.md`, `research.md`, and `aiuse.md` with message "Add robustness analysis, research notes, and aspect ratio warning". Pushed to `origin/main`.

---

## Prompt 14: Comments on plan.md

**User:** Provided detailed comments on each suggestion in plan.md:
- 2.1, 2.2: Accepted as is.
- 2.3: Intentional footgun for testing, leave but don't act.
- 2.4: 10% is generous, leave alone.
- 2.5: Explain more verbosely.
- 2.6: Plan on adding tol parameter to snapToSymmetry_SVD.
- 2.7: Plan to implement, but asked whether high aspect ratios are the only failure mode.
- 2.8, 2.9: Leave in place, no action.
- 2.10: Good, but is it computationally expensive?
- 2.11: Good, employ later.
- 3.1: Good, will do.
- 3.2: Explain how it differs from the current test.
- 3.3, 3.4: Good as is.
- 3.5: Explain rationale.
- 3.6: Good as is.
- 3.7: Explain rationale more.
- 3.8: Explain utility.
- 3.9: Good, document rationale better.
- 3.10: Good, but variants will be removed eventually.
- 3.11, 3.12: Good as is.
- Section 4: Good, leave for later.

**Response:** Updated plan.md with all feedback:
- Added status markers (Accepted/Future/No action) to each suggestion.
- Expanded 2.5 with detailed failure mode explanation and concrete fix.
- Noted 2.7 also applies to near-degenerate bases and angle-boundary cases.
- Answered 2.10: O(n²), n≤48, negligible cost.
- Replaced 3.2 with "Structured noise patterns" (axis-aligned, angular, length-only).
- Added rationale for 3.5 (validates similarity transformation precision).
- Renamed 3.7 to "Non-primitive cell invariance" with expanded rationale.
- Expanded 3.8 with pseudosymmetry explanation and added orthorhombic-near-tetragonal test.
- Added rationale to 3.9 (documents working envelope, shows improvement from fixes).
- Noted 3.10 variants will be removed.

---

## Prompt 15: Comments on research.md

**User:** Provided detailed comments asking for expansions:
- Section 1: Explain more about the metric tensor.
- Section 2: Not every code does reduction (VASP uses hard-coded bounds). Discuss how noise affects Minkowski reduction.
- Section 3.1: More on spglib tolerance, known failures. Only cite verifiable references.
- Section 3.2: Can AFLOW-SYM over-identify? Boundary cases? Does it snap or pick nice orientations?
- Section 3.4: More on SOFI. Could ideas be adapted?
- Section 3.5: More on Symmetrizer. How could it help?
- Section 5.1: Pitfalls of adaptive tolerance in Spacey.jl. Over-identification risk.
- Section 5.3: More on metric tensor and rational arithmetic.
- Section 5.4: More details on interval arithmetic.
- Section 5.5: More details on eigenvalue-based frames.

**Response:** Researched all topics online for verifiable information, then rewrote research.md with substantial expansions to all requested sections. Key additions: noise effects on Minkowski reduction, spglib crash statistics, AFLOW-SYM tolerance landscape table, SOFI's non-principal-axis design, bottom-up symmetry detection workflow, over-identification safeguards, rational arithmetic two-step process, interval arithmetic step-by-step for Spacey.jl with costs/benefits, eigenvalue-based frames with degeneracy analysis.

---

## Prompt 16: Summarize session and commit

**User:** Asked to summarize all prompts and responses, store in `aisummary.md`, add to repo, and commit all changes.

**Response:** This file. Committed along with updated plan.md and research.md.
