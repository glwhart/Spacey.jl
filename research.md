# Literature Review: Point Group and Space Group Symmetry Detection

This document surveys algorithms, software packages, and mathematical results relevant to the crystallographic point-group-finding problem that Spacey.jl solves. The goal is to surface ideas that could make Spacey.jl more robust and to lay groundwork for future space-group functionality.

---

## 1. The Core Mathematical Problem

Given three basis vectors **a₁, a₂, a₃** ∈ ℝ³ defining a Bravais lattice Λ, find all linear maps R : ℝ³ → ℝ³ such that R(Λ) = Λ and R is an isometry (i.e., ‖Rv‖ = ‖v‖ for all v). The set of all such R forms the **point group** G(Λ), a finite subgroup of O(3).

### The Metric Tensor Formulation

The **metric tensor** (or Gram matrix) of a lattice basis A = [**a₁** | **a₂** | **a₃**] is:

> **G = AᵀA**

Its entries are the dot products between basis vectors: G_ij = **aᵢ** · **aⱼ**. The diagonal entries are squared lengths (G_ii = ‖**aᵢ**‖²) and the off-diagonal entries encode angles (G_ij = ‖**aᵢ**‖ ‖**aⱼ**‖ cos θ_ij).

The metric tensor encodes the full geometry of the lattice—lengths, angles, and volume—without reference to any external coordinate system. Two bases A and B describe the same lattice if and only if there exists an integer matrix M with det(M) = ±1 such that B = AM, in which case their metric tensors are related by G_B = MᵀG_A M.

The point group in lattice coordinates is exactly the **automorphism group of the metric tensor**:

> **Aut(G) = { M ∈ GL(3,ℤ) : MᵀGM = G }**

This is because MᵀGM = G expands to (AM)ᵀ(AM) = AᵀA, which says the candidate basis AM has the same metric as A—meaning the map A → AM preserves all distances and angles, i.e., it is a symmetry.

**Why the metric tensor matters algorithmically:**
- It reduces the problem from ℝ³ (continuous, orientation-dependent) to the space of 3×3 symmetric positive-definite matrices (6 independent parameters, orientation-independent).
- The condition MᵀGM = G is a finite system of integer equations once G is known—no transcendental functions, no trigonometry.
- The metric tensor is well-defined even for degenerate or nearly-degenerate bases where the basis matrix A is ill-conditioned.

**Relationship to Spacey.jl's current approach:** Spacey.jl works in Cartesian coordinates and checks R = A'A⁻¹ for orthogonality (RᵀR ≈ I). This is mathematically equivalent to checking A'ᵀA' ≈ AᵀA (i.e., same metric tensor). As discussed elsewhere, the T ≈ I check is actually better-conditioned for high aspect ratios because all entries of T are ~1, whereas the metric tensor has entries spanning a²...b² for aspect ratio b/a.

However, the metric tensor formulation has distinct advantages for **lattice classification** (identifying the Bravais type) and for **rational arithmetic approaches** (see Section 5.3), where working with G directly avoids the need to ever compute floating-point square roots or inversions.

---

## 2. Why Lattice Reduction Is Fundamental (But Not Universal)

Most modern algorithms for lattice symmetry start with lattice reduction, because the bounds that make the search finite and computationally tractable depend on it. However, **not all codes require reduction.** For example, the VASP code (Vienna Ab initio Simulation Package) uses hard-coded search bounds over a large range of integer coefficient combinations—this is inefficient but fast enough in practice for typical DFT cells, and avoids dependence on a reduction algorithm. The downside is that the hard-coded bounds are not provably complete for all possible inputs: pathological cases (extreme aspect ratios, very non-orthogonal cells) could in principle escape the search window.

Lattice reduction guarantees completeness: for a Minkowski-reduced basis, the search over ±1 integer coefficients is provably exhaustive, and the search space is minimal (27 candidate vectors). This is the approach Spacey.jl takes.

### 2.1 Minkowski Reduction

A basis {**b₁**, **b₂**, **b₃**} is Minkowski-reduced if each basis vector is as short as possible subject to the constraint that it must complete the previous vectors to a basis. More precisely, the 12 Minkowski conditions must hold:

- **b₁** is the shortest nonzero lattice vector: ‖**b₁**‖ ≤ ‖**b₂**‖ and ‖**b₂**‖ ≤ ‖**b₃**‖  
- **b₃** is at least as short as all 12 combinations **b₃ ± b₁**, **b₃ ± b₂**, **b₃ ± b₁ ± b₂**  
- **b₂** is at least as short as **b₂ ± b₁**

These are exactly the 12 conditions checked by `isMinkReduced` in `MinkowskiReduction.jl`. The key property they guarantee is that each basis vector is shorter than any other vector that occupies the same "position" in the lattice (the same coset of the sublattice spanned by the shorter vectors). This is the reduction used by Spacey.jl (via the `MinkowskiReduction.jl` package).

Key property for point-group algorithms: a Minkowski-reduced basis satisfies tight geometric bounds on the angles and length ratios. Specifically, for any two basis vectors: cos(θᵢⱼ) ≤ 1/2, meaning all angles between basis vectors lie in [60°, 120°]. This is exactly the condition needed to prove that any symmetry operation maps basis vectors to ±1 linear combinations of themselves—the 27-neighbor search is sufficient.

**Computational note:** Minkowski reduction is NP-hard in `n` dimensions but polynomial (O(n³)) in fixed dimension 3, and practical implementations (including `MinkowskiReduction.jl`) are fast. However, floating-point implementations lose precision when vector lengths span many orders of magnitude—this is the source of Spacey.jl's known failure at extreme aspect ratios (~2²⁵).

#### How numerical noise affects Minkowski reduction

Minkowski reduction is a discrete optimization problem (find the shortest lattice vectors), but its implementation involves floating-point comparisons of vector lengths and angles. Numerical noise in the input basis vectors can corrupt the reduction in several ways:

1. **Length-comparison instability:** The reduction algorithm must decide which of two vectors is shorter. When two candidate vectors have nearly equal lengths (differing by less than floating-point precision), the algorithm may choose arbitrarily, leading to different "reduced" bases for inputs that differ by noise below machine epsilon. This is most problematic for lattices where multiple vectors have the same length (cubic, rhombohedral).

2. **Angle-bound violations:** The Minkowski bounds require angles in [60°, 120°]. When the true angle is exactly 60° or 120° (which happens for hexagonal and rhombohedral lattices), noise can push the computed angle just outside the bound, causing the reduction to perform an unnecessary swap. The resulting basis is in a different orientation, but this alone does not cause symmetry-finding to fail — Spacey.jl's search over all 27 candidate vectors doesn't depend on any particular orientation. The real risk is that the swap produces a basis that no longer satisfies all 12 Minkowski conditions. If any condition is violated, the {-1,0,1}³ search space becomes incomplete (see point 4 below).

3. **Non-idempotence:** Ideally, reducing an already-reduced basis returns the same basis. But noise can make `minkReduce(minkReduce(u,v,w))` differ from `minkReduce(u,v,w)` because the first reduction changes the basis slightly (reordering), and the second reduction sees a slightly different numerical input. This is the issue Spacey.jl guards against with the `orthogonalityDefect` check.

4. **Propagation to symmetry finding:** If the reduction produces a basis that doesn't tightly satisfy the Minkowski bounds (due to any of the above), the 27-neighbor search may be incomplete. A symmetry operation that maps a basis vector to a ±2 combination (which shouldn't be possible for a truly reduced basis) would be missed.

   **Concrete example:** Consider a 2D square lattice with the non-reduced basis **a₁** = [1, 0], **a₂** = [3, 1]. The true reduced basis is **a₁** = [1, 0], **a₂** = [0, 1]. The 4-fold rotation (90°) maps [1, 0] → [0, 1] in Cartesian. But in the non-reduced basis, [0, 1] = -3**a₁** + **a₂**, so the operation requires a coefficient of -3. The {-1,0,1}² search would miss it entirely. In 3D, a hexagonal lattice with a non-reduced choice of second basis vector (e.g., **a₁** + **a₂** instead of **a₂**) produces the same situation: the 6-fold rotation maps the long vector to a combination with coefficient -2, which falls outside the search window.

   **This is the strongest argument for verifying that Minkowski reduction actually succeeded**, using `isMinkReduced` as a post-condition check. A unit test that deliberately supplies a nearly-reduced (but not quite reduced) basis and checks that the point-group finder either (a) gives the correct answer or (b) emits a warning, would expose this failure mode. No such test currently exists in `test/runtests.jl`.

**Approaches to making reduction more robust:**

- **Tolerance-aware reduction:** Treat vectors as "same length" when they differ by less than a threshold. Use a deterministic tie-breaking rule (e.g., lexicographic ordering of components) when lengths are indistinguishable. This prevents arbitrary choices driven by noise.
- **Verify post-conditions:** After reduction, explicitly check that all angles are in [60°-δ, 120°+δ] for some small δ. If violated, the basis is at a degenerate boundary and the symmetry finder should extend its search radius (see plan.md §2.7).
- **Canonical ordering:** Always sort the reduced basis vectors by length (then by a deterministic secondary criterion), so that the output is unique regardless of input noise.
- **Snap-then-reduce:** For noisy inputs, run a preliminary symmetry search, snap to symmetry, then re-reduce. The snapped basis will be exactly on the Minkowski bounds, giving a cleaner reduction. **However, this strategy is hazardous at symmetry boundaries.** A lattice near the cubic/tetragonal boundary (c/a ≈ 1) could be genuinely tetragonal due to physics (symmetry breaking), or it could be cubic with numerical noise making c/a depart from 1. Snapping to the tetragonal group and re-reducing would destroy the cubic symmetry in the first case — which is correct. But it would also destroy it in the second case — which is an error. The algorithm has no way to resolve this without physics input. The honest behavior is to report the ambiguity (which symmetry group is consistent, and at what tolerance) rather than silently committing to one side.

### 2.2 Niggli Reduction

The **Niggli reduced cell** (Niggli 1928; Křívý & Gruber 1976) provides a *canonical* representative for every lattice—the unique reduced form under a precise set of main conditions and special conditions. While Minkowski reduction gives the *shortest* basis, Niggli reduction gives the *unique* description.

The key advantage for symmetry work: once a Niggli cell is computed, the Bravais lattice type can be determined by a table lookup using the metric tensor components, and the point group follows from the lattice type. Modern implementations (e.g., Andrews & Bernstein 2022, *J. Appl. Cryst.*) have improved the numerical stability of Niggli reduction, particularly for handling near-degenerate cases where multiple Niggli forms are nearby in parameter space.

**Relevance to Spacey.jl:** Niggli reduction is an alternative to Minkowski reduction. It could be used as a preprocessing step or as a cross-check. The NIST table-lookup approach after Niggli reduction is orthogonal to Spacey.jl's exhaustive-search approach and could serve as an independent validation.

### 2.3 LLL Reduction

The Lenstra–Lenstra–Lovász (LLL) algorithm is a polynomial-time approximation to Minkowski reduction. It does not find the theoretically shortest basis but is significantly faster for high-dimensional lattices and has well-understood numerical behavior. For 3D crystallography, true Minkowski reduction is both tractable and preferable; LLL is only be relevant in more than four dimensions.

---

## 3. Major Software Packages

### 3.1 Spglib

**Language:** C with Python/Julia/Fortran wrappers  
**Reference:** Togo & Tanaka, *Sci. Technol. Adv. Mater.: Methods* 2024; arXiv:1808.01590

Spglib is the most widely used crystallographic symmetry library. Its algorithm:

1. Finds the primitive cell using Delaunay or Niggli reduction.
2. Searches for all translation symmetries (to find symmorphic operations).
3. Tests candidate point-group operations using a position-based tolerance: two atomic positions are equivalent if their Euclidean distance is less than `symprec` (default 1×10⁻⁵ Å).
4. For space groups: identifies Wyckoff positions and the Hermann–Mauguin symbol.

#### Tolerance strategy in detail

According to the arXiv paper and documentation, spglib's tolerance adjustment is not a simple binary search. The algorithm has multiple stages, and tolerance is adjusted differently at each:

- **Primitivization stage:** The algorithm tries to find pure translations that relate equivalent atoms. Sometimes, the initial tolerance fails because candidate translations don't produce a lattice volume consistent with the number of translations found. When this happens, the tolerance is *reduced* and the search is repeated.
- **Symmetry search stage:** Candidate rotations are tested against the atomic structure. If too few or too many operations are found, the tolerance can be adjusted.
- **Consistency check:** After finding a candidate space group, spglib verifies internal consistency (e.g., that the found operations form a group, that Wyckoff positions are consistent). If this fails, tolerance is reduced and the process restarts.
- **Dual-stage optimization:** Some pure translations that don't satisfy constraints are discarded by re-examining with tightened tolerance, which is faster than restarting from scratch.

The key design choice is that **tolerance is generally only decreased, not increased.** The philosophy is to start with the user-specified tolerance and tighten it until a consistent answer is found. This means spglib will never report *higher* symmetry than consistent with the given tolerance, but it may report *lower* symmetry if the iterative tightening overshoots.

#### Known failure modes

Based on the spglib GitHub issue tracker and published benchmarks:

- **Tolerance-dependent crashes:** Testing on 194,017 crystal structures showed crash rates increase sharply at higher tolerances (137 crashes at 0.1 Å vs. 0 at 0.00001 Å). The crash path typically involves the symmetry search finding too many candidates, leading to memory or index-bounds issues.
- **Threshold parameter bug (post-1.14.1):** `angle_tolerance` parameter was silently ignored by the symmetry operation search while still being used by group identification, leading to inconsistent results.
- **Unit-dependence:** Because `symprec` is an absolute tolerance in Ångströms, it doesn't scale with the size of the unit cell. A tolerance of 10⁻⁵ Å is tight for a protein crystal (100 Å cell) but generous for a hydrogen molecule (0.7 Å bond).
- **Heuristic decisions:** Several internal tolerance criteria are "determined heuristically" (their words) without rigorous mathematical justification. The paper acknowledges that applying tolerance to mathematically perfect structures is inherently ill-defined.
- **Near-degenerate lattices:** When the primitive cell search produces a cell at the boundary between two Bravais types, the iterative tolerance reduction can oscillate or converge to a subgroup.

**Key difference from Spacey.jl:** Spglib works on *atomic structures* (positions + species) rather than bare lattices. Its tolerance applies to atomic positions, not to the lattice itself. For purely lattice point groups, Spacey.jl's approach (tolerance applied to geometric properties of the lattice) is more direct and avoids the complexity of atom matching.

### 3.2 AFLOW-SYM

**Language:** C++ (part of the AFLOW framework)  
**Reference:** Hicks et al., *Acta Cryst. A* 74, 184–203 (2018); arXiv:1802.07977

AFLOW-SYM's primary contribution is **self-consistent tolerance determination**:

> Rather than asking the user to supply a tolerance, AFLOW-SYM searches for a tolerance ε such that the set of symmetry operations found at ε is consistent with a valid group (closed under composition, contains identity and inverses).

#### Algorithm and boundary behavior

AFLOW-SYM scans across a range of tolerance values and identifies which space groups are consistent at each tolerance. The paper presents detailed case studies showing the tolerance landscape. For example, in the AgBr structure:

| Tolerance range (Å) | Consistent space group |
|---------------------|----------------------|
| 1.0×10⁻⁶ to 4.09×10⁻² | Monoclinic (No. 11) |
| 4.09×10⁻² to 1.64×10⁻¹ | Orthorhombic (No. 59) |
| 1.64×10⁻¹ to 2.46×10⁻¹ | **Gap: no consistent group** |
| 2.46×10⁻¹ to 6.70×10⁻¹ | Rhombohedral (No. 166) |
| 6.70×10⁻¹ to 1.0 | FCC cubic (No. 225) |

The "gap" regions are where no self-consistent group exists—the candidates found at that tolerance don't close under multiplication. This is an important feature: it means the algorithm can explicitly *refuse* to assign a symmetry when the input is ambiguous.

#### Risk of over-identification

**Can AFLOW-SYM find a group that is too large?** Yes, in principle. The self-consistent tolerance is the *loosest* tolerance for which the candidates form a valid group. If the input structure has noise that happens to look like a higher-symmetry distortion, the algorithm will find the higher symmetry as "self-consistent" at a loose tolerance.

**Protection mechanisms:**
- The algorithm reports the *tightest* self-consistent tolerance (smallest ε where a group is found), not the loosest. This favors the lowest consistent symmetry.
- The "gap" regions act as natural barriers between symmetry assignments. At a boundary, neither the higher nor lower group is self-consistent, forcing the algorithm to pick one side.
- Validation against 54,000+ ICSD entries provides empirical evidence that over-identification is rare.

**However**, for near-boundary cases (structures genuinely between two symmetries), the algorithm's choice is inherently ambiguous. The paper documents this as a feature (the tolerance itself characterizes "how symmetric" the structure is), not as a failure mode. But for a user who needs a definitive answer, the boundary behavior is a limitation.

#### Snapping and standard settings

The AFLOW-SYM paper describes the algorithm as providing "various representations for the point, factor and space groups, site symmetries and Wyckoff positions." The algorithm includes transformation to standard settings (ITA conventions). It's unclear from the available references whether AFLOW-SYM performs explicit symmetrization (snapping the structure to perfect symmetry) as a separate step, or whether the self-consistent tolerance implicitly defines the symmetrized structure. The focus appears to be on *identifying* symmetry rather than *enforcing* it on the coordinates.

### 3.3 FINDSYM (ISOTROPY Software Suite)

**Reference:** Stokes & Hatch, *J. Appl. Cryst.* 38 (2005)

FINDSYM identifies space groups from crystal structures and transforms them to a standard setting. Unlike Spglib/AFLOW, FINDSYM provides symmetry in a standardized basis (ITA convention). Recently extended to handle **superspace groups** for modulated/quasicrystalline structures.

**Key algorithmic idea:** FINDSYM uses the Niggli cell for reduction, then a table lookup for the Bravais type, then an exhaustive search over the corresponding point-group operations. The table-lookup approach is fast and deterministic, but relies on first classifying the Bravais type correctly (which itself requires handling near-degenerate cases). Spacey.jl's approach differs: it searches exhaustively over all possible symmetry operations without first committing to a Bravais type, making it more robust for inputs that are near a boundary between two lattice types.

### 3.4 SOFI: Point Groups by Shape Matching

**Reference:** Gunde et al., *J. Chem. Phys.* 161, 062503 (2024); arXiv:2408.06131

SOFI (Symmetry Operations FInder) formulates point-group detection as a **degenerate shape-matching problem**. Given a set of points (atomic positions), it finds rotation/reflection operations that map the set to itself.

#### Algorithm details

The core insight is that a symmetry operation is a *degenerate solution* to the shape-matching problem: there are multiple rotations R that satisfy R(structure) = structure (up to permutation of equivalent atoms). The algorithm:

1. **Formulate shape matching:** Given two copies of the atomic cluster, find the optimal rotation R that maps one to the other (minimizing RMSD). This is a standard problem solvable by SVD (Kabsch algorithm).
2. **Exploit degeneracy:** For a symmetric structure, the shape-matching problem has multiple solutions (one for each symmetry operation). All solutions have RMSD = 0 for a perfect structure, or RMSD < threshold for a distorted one.
3. **Enumerate solutions:** Rather than finding just the optimal match, find ALL rotations within a threshold of zero RMSD. Each such rotation is a symmetry operation.
4. **Avoid principal axes:** Unlike earlier methods (SymMol, Symmetrizer), SOFI does *not* use the inertia tensor to define a reference frame. The paper explicitly notes that principal axes "can be ambiguous in some specific cases (spheroidal shapes)." Instead, SOFI searches along all possible axes.
5. **Non-modifying:** SOFI returns symmetry operations in the original reference frame without modifying or symmetrizing the input.

#### Adaptation potential for Spacey.jl

SOFI is designed for molecular/cluster point groups (finite number of atoms, no periodicity). Direct adaptation to periodic lattices is not discussed in the paper. However, several ideas could transfer:

- **Threshold-as-quality:** SOFI interprets the RMSD of each symmetry operation as a quality measure. Operations with lower RMSD are "better" symmetries. This is analogous to Spacey.jl's `norm(T - I(3))` sorting, and the idea of returning a quality metric (plan.md §2.11).
- **Avoid canonical frames:** SOFI's design choice to avoid inertia-tensor-based frames (because of degeneracy issues) is relevant. For cubic lattices, any eigenvector-based frame is degenerate. SOFI's approach of directly searching over all candidate operations, rather than first defining a frame, is closer to what Spacey.jl already does.
- **Exhaustive-then-filter:** SOFI generates all candidate rotations within a threshold, then identifies which ones form a group. This is the same philosophy as Spacey.jl's `pointGroup_robust`.

The main limitation for lattice adaptation is that SOFI works with discrete point sets, while a lattice is infinite. Spacey.jl's use of Minkowski reduction to bound the search space is the key piece that makes the problem finite for lattices—this is something SOFI doesn't need (for molecules).

### 3.5 Symmetrizer and SYMMOL: Bottom-Up Symmetry Element Detection

**References:**
- Largent et al., *J. Comput. Chem.* 33, 2249 (2012) — Symmetrizer
- Pilati & Forni, *J. Appl. Cryst.* 31, 557 (1998) — SYMMOL

These programs take a fundamentally different approach from Spacey.jl. Rather than generating all candidate operations and filtering, they **detect individual symmetry elements first** and then compose them to determine the full group.

#### How the bottom-up approach works

1. **Detect candidate rotation axes:** For each pair of equivalent atoms (or, for lattices, lattice vectors), the midpoint and connecting vector define a candidate mirror plane or rotation axis. For a lattice, the short vectors after Minkowski reduction define natural candidate axes.

2. **Test each candidate element:** For each candidate axis/plane, check whether the corresponding rotation/reflection maps the structure to itself (within tolerance). Score each element by its deviation from perfect symmetry.

3. **Identify the point group:** The detected elements constrain the possible point groups. For example, if a 4-fold axis and two perpendicular 2-fold axes are detected, the group must be D₄ₕ or a subgroup. The algorithm selects the largest group consistent with all detected elements.

4. **Quantify fit quality:** Because each element is tested independently, the deviation for each element is known. This provides per-element quality information that an exhaustive search doesn't give.

#### How this could make Spacey.jl more robust

The bottom-up approach has complementary strengths to Spacey.jl's current top-down approach:

- **Faster failure diagnosis:** If the algorithm fails to find the expected group, it can report *which specific elements* are missing and by how much they miss the tolerance. This is more informative than "found 12 operations instead of 48."
- **Partial symmetry:** For structures with broken symmetry (defects, distortions), the bottom-up approach naturally finds the remaining symmetry elements without needing to find a complete group.
- **Independent validation:** Detecting elements bottom-up and composing them to a group could serve as a cross-check against the exhaustive search. If both methods agree, confidence is high.
- **Reduced search for high symmetry:** For cubic symmetry (48 operations), detecting a 4-fold axis, a 3-fold axis, and inversion is sufficient to determine the full group. This is 3 tests instead of checking 48 candidate operations.

**Limitations for lattice application:** The bottom-up approach requires enumerating candidate symmetry axes/planes, which for a lattice comes from the lattice vectors themselves. After Minkowski reduction, there are at most 26 non-zero candidate directions (the 27 neighbors minus the zero vector). Testing each as a rotation axis of order 2, 3, 4, or 6 is cheap (at most 4×26 = 104 tests). This could be faster than the current approach for high-symmetry lattices.

---

## 4. Algorithms for Space Groups (Future Spacey.jl Work)

### 4.1 The Space Group Problem

The **space group** G of a crystal is an extension of the point group by lattice translations:

> G = { (R, t) : R ∈ G(Λ), t ∈ ℝ³, (R, t) maps the crystal to itself }

Finding the space group requires knowing the atomic basis (fractional coordinates and atomic species). There are 230 space groups in 3D.

The relationship to point groups: every space group has an associated point group obtained by ignoring translations. The point group is a factor group of the space group. Spacey.jl already finds the point group; the additional work for space groups is:

1. Find the set of fractional translations `t` such that `(R, t)` is a symmetry for each point-group operation R.
2. Classify the resulting set of (R, t) pairs into one of the 230 space-group types.

### 4.2 Wyckoff Positions and Standardization

To identify the space group uniquely (not just as an abstract group), the structure must be transformed to the standard setting (as in the International Tables for Crystallography, ITA). This requires:

1. Finding the **transformation matrix** from the current basis to the ITA standard basis.
2. Identifying the **Wyckoff positions** occupied by each atomic species.
3. Matching the resulting (lattice type, Wyckoff occupancies) to the ITA table.

FINDSYM and Spglib both implement this; it would be a significant but well-understood addition to Spacey.jl.

### 4.3 Magnetic Space Groups

For magnetic materials, the relevant symmetry group includes time-reversal operations. There are 1651 magnetic space groups (Shubnikov groups). FINDSYM has been extended to handle these. This is a direction worth keeping in mind for future Spacey.jl development.

---

## 5. Numerical Robustness Techniques from the Literature

### 5.1 Adaptive Tolerance (Spglib / AFLOW-SYM)

Both Spglib and AFLOW-SYM converge on the idea that a fixed tolerance is fragile. The right tolerance depends on the structure:
- For highly symmetric, perfect crystals: very tight (10⁻¹⁰)
- For DFT-relaxed structures with residual forces: moderate (10⁻⁵)
- For experimental structures with disorder: loose (10⁻²)

The self-consistent approach (find the largest ε for which the candidates form a group) is optimal in the sense that it finds the highest possible symmetry consistent with the input precision.

#### Implementation in Spacey.jl: pitfalls and challenges

An adaptive/self-consistent tolerance for Spacey.jl (plan.md §2.8) would work roughly as: sweep tolerance from large to small, at each value run the point-group finder, return the result from the finest tolerance that gives a valid group. However, there are specific dangers:

**Over-identification risk:** The self-consistent approach inherently favors the highest consistent symmetry. For a slightly distorted cubic cell (c/a = 1.0001), a loose tolerance will see cubic (48 ops) as self-consistent. The algorithm must be designed to prefer the *tightest* tolerance (lowest symmetry), not the loosest.

**Proposed safeguard:** Start from the tightest tolerance and increase, rather than starting loose and decreasing. Return the *first* tolerance at which a valid group is found. This way:
- For a perfect cubic cell: even the tightest tolerance gives 48 ops.
- For a slightly tetragonal cell: tight tolerance gives 16 ops, which is already a valid group, so the algorithm stops there and never "discovers" the spurious cubic at a looser tolerance.

**Boundary oscillation:** Near the boundary between two symmetry assignments (e.g., a structure that is barely tetragonal vs barely cubic), the algorithm could oscillate: at tolerance ε₁ it finds 16, at ε₁+δ it finds 48, at ε₁-δ it finds 16 again. The AFLOW-SYM paper documents "gap" regions where no consistent group exists. For Spacey.jl, this could manifest as the algorithm returning different answers for nominally identical inputs with different floating-point rounding.

**Proposed safeguard:** If the group size changes between consecutive tolerance values, flag the result as "borderline" and report both possibilities to the user. This is more honest than silently picking one.

**Computational cost:** Running the point-group finder multiple times (one per tolerance value) multiplies the cost. For Spacey.jl, where a single run takes ~100μs, this is acceptable (10 runs = 1ms). But it would need to be benchmarked.

### 5.2 Metric Tensor Approach (Expanded)

The metric tensor **G = AᵀA** is a positive-definite 3×3 symmetric matrix. Point-group operations are exactly the integer matrices M with MᵀGM = G.

**Working directly with G has several advantages for certain sub-problems:**

1. **Orientation independence:** G doesn't change if you rotate the entire lattice. Spacey.jl is also orientation-independent (after Minkowski reduction, the candidate search and group closure test don't depend on the initial orientation). The distinction is that G makes this orientation-independence *explicit*: because G has no rotational degrees of freedom at all, an algorithm based directly on G can work with a 6-parameter representation (the 6 independent entries of the symmetric matrix) rather than a 9-parameter matrix A whose 3 extra degrees of freedom correspond to rigid rotations. For the specific problem of Bravais type classification, this is a real advantage: checking which entries of G are equal is simpler than checking which angles between Cartesian basis vectors are equal.

2. **Natural norm for comparison:** The Frobenius distance between two metric tensors ‖G₁ - G₂‖_F is a natural measure of "how different" two lattices are. This is used by the minimum-strain symmetrization approach (Phys. Rev. Research 2, 013077, 2021) to find the closest lattice of a given Bravais type.

3. **Rational arithmetic compatibility:** If the input basis vectors have rational coordinates (or can be scaled to have integer coordinates), then G has rational entries, and the condition MᵀGM = G can be checked in exact rational arithmetic. See Section 5.3.

4. **Bravais type determination:** After Niggli reduction, the Bravais type can be determined by checking which constraints the metric tensor satisfies (e.g., G₁₁ = G₂₂ for tetragonal, G₁₁ = G₂₂ = G₃₃ for cubic). This is the "table lookup" approach used by FINDSYM.

**Why Spacey.jl's current approach (T ≈ I) is still better for the symmetry search itself:**

As discussed in Section 1, checking RᵀR ≈ I is better-conditioned for high aspect ratios because all entries of T = RᵀR are ~1 (for R approximately orthogonal), making the `isapprox` comparison uniform. The metric tensor G has entries spanning a²...b², causing `isapprox` with `rtol` to be dominated by the large entries. For the *search* step, T ≈ I wins. But for *validation*, *classification*, and *rational arithmetic*, G is the right object.

### 5.3 Exact Arithmetic with Integer Coefficients

After Minkowski reduction, all point-group operations have integer matrix representations. The group closure test `isagroup(LG)` in Spacey.jl already exploits this by using exact integer arithmetic. This is the strongest possible validation—no tolerance needed, no false positives.

**Extending exact arithmetic to the search step:**

The metric tensor G = AᵀA has entries that are dot products of basis vectors. If the basis vectors are given as floating-point numbers, G is also floating-point. However, for the purpose of checking MᵀGM = G for integer M, we can:

1. Compute G from the input (floating-point).
2. For each candidate integer matrix M, compute MᵀGM exactly (since M is integer, this is just integer linear combinations of the entries of G—the floating-point operations are additions and multiplications of G entries, with no division or square roots).
3. Check if MᵀGM ≈ G component-wise.

This is subtly different from Spacey.jl's current approach, which computes `inv(A)` (introducing floating-point error proportional to κ(A)) and then rounds to integer. The metric tensor approach never inverts A—it only multiplies G entries by integers.

**Rational arithmetic approach (used in mathematical software like Sage):**

If the metric tensor G has *rational* entries (which happens when the basis vectors have rational coordinates—common in crystallographic databases where coordinates are given as fractions), then the entire computation can be done in exact rational arithmetic:

1. Represent G as a matrix of `Rational{BigInt}` values.
2. For each candidate integer matrix M, compute MᵀGM in exact rational arithmetic.
3. Check MᵀGM == G exactly (no tolerance needed!).

This completely eliminates floating-point error. The cost is ~10× slower than Float64 operations (rational arithmetic involves GCD computations), but for a 3×3 matrix with only 48 candidate operations to test, the total cost would be negligible (microseconds).

**Practical applicability:** Most real-world inputs to Spacey.jl are floating-point (from DFT calculations, experiments), so exact rational arithmetic isn't directly applicable to the raw input. However, it could be used in a two-step process:
1. Use the current floating-point algorithm to find candidate integer operations.
2. Rationalize G (round entries to nearby rationals using `rationalize(g, tol=...)`).
3. Verify the candidates against the rationalized G using exact arithmetic.

This gives the speed of floating-point search with the certainty of exact validation. The only risk is that `rationalize` introduces error if the input G isn't actually close to a rational matrix—but for reasonable lattice parameters, this should always work.

### 5.4 Interval Arithmetic

**The idea:** Replace all floating-point numbers with *intervals* [a, b] that are guaranteed to contain the true value. Every arithmetic operation propagates the intervals correctly (e.g., [a,b] + [c,d] = [a+c, b+d], [a,b] × [c,d] = [min(ac,ad,bc,bd), max(ac,ad,bc,bd)]). The final result is an interval guaranteed to contain the correct answer. If the interval is narrow enough to make a decision (e.g., "this determinant is in [0.99, 1.01], which is ≈ 1"), the result is rigorous. If the interval is too wide, the algorithm says "I don't know" rather than giving a wrong answer.

**Relevance to lattice symmetry:**

For Spacey.jl's algorithm, interval arithmetic would work as follows:

1. Represent input basis vectors as intervals: `u = [u₁±ε, u₂±ε, u₃±ε]` where ε reflects measurement/computation uncertainty.
2. Compute the Minkowski reduction using interval arithmetic. The reduction algorithm involves comparisons (`is |v₁| < |v₂|?`); with intervals, a comparison may be *undecidable* (the intervals overlap). In that case, both branches must be explored.
3. Compute candidate operations: since the 27 neighbors are integer combinations, this is exact (integers have zero-width intervals). The matrix products and comparisons propagate the input uncertainty.
4. Check orthogonality: `T = RᵀR` is an interval matrix. Check whether `I(3)` is contained in the interval matrix T. If yes: the operation is definitely a symmetry. If no: definitely not. If the intervals overlap with I(3) boundary: undecidable—need tighter input or more precision.

**Benefits:**
- No tolerance parameter needed—the algorithm derives its own "tolerance" from the input uncertainty.
- Guaranteed correctness: if the algorithm returns an answer, it is provably correct for the given input intervals.
- Explicit "I don't know" answers for genuinely ambiguous cases, rather than a potentially wrong definitive answer.

**Challenges and costs:**
- **Interval explosion:** Each arithmetic operation can widen intervals. For a sequence of multiplications (like matrix products), intervals grow exponentially in the worst case. For a well-conditioned 3×3 matrix, this is manageable, but for ill-conditioned cases (high aspect ratio), the intervals may grow too wide to be useful.
- **Branch explosion in reduction:** If the Minkowski reduction can't decide which vector is shorter (overlapping intervals), it must explore both choices. For degenerate lattices with many nearly-equal vectors, this could be exponentially expensive.
- **Computational overhead:** Interval arithmetic costs 2–4× per operation (each number is now two floats, and each operation requires rounding-mode control). For Spacey.jl's already-fast algorithm (~100μs), this is acceptable.
- **Julia ecosystem support:** `IntervalArithmetic.jl` is a mature, IEEE 1788-2015 compliant package. It supports all basic operations and linear algebra. Using it with Spacey.jl would require making the algorithm generic over number types (using Julia's parametric type system), which is a moderate refactoring effort.

**Existing work:** Espitau & Joux (arXiv:1905.11743) developed a certified lattice reduction algorithm using interval arithmetic for the LLL algorithm. Their key result: the cost of interval arithmetic with floating-point bounds is at most 4× that of classical large-precision floating-point. They use adaptive precision—starting with narrow intervals and widening only when needed. This same adaptive approach could be applied to Spacey.jl's Minkowski reduction.

**Practical recommendation for Spacey.jl:** Rather than converting the entire algorithm to interval arithmetic (high effort, moderate benefit), a targeted application would be valuable:
- Use interval arithmetic *only* for the Minkowski reduction step, which is where precision loss causes the known failures at high aspect ratios.
- Or: use interval arithmetic as a *post-hoc validation*—after finding candidate operations with standard floats, re-check them with intervals to get rigorous bounds on "how close to orthogonal" each operation really is.

### 5.5 Eigenvalue-Based Canonical Frames

**The idea:** The metric tensor G = AᵀA is a symmetric positive-definite matrix. Its eigenvectors define a natural coordinate frame for the lattice. In this frame, G is diagonal (G = diag(λ₁, λ₂, λ₃) where λᵢ are eigenvalues, i.e., squared lengths along principal axes). Symmetry operations in this frame take a simpler form.

**How it works in detail:**

1. **Diagonalize G:** Compute G = VDVᵀ where D = diag(λ₁, λ₂, λ₃) and V is orthogonal (columns are eigenvectors).

2. **Transform to principal frame:** The basis A_principal = A·V has metric tensor G_principal = Vᵀ·G·V = D (diagonal). In this frame, the basis vectors are along the eigenvectors of G, and their lengths are √λᵢ.

3. **Simplify the symmetry search:** A point-group operation M (in lattice coordinates) transforms G as MᵀGM = G. In the eigenbasis, this becomes Mᵀ_p D M_p = D, where M_p = Vᵀ M V. If all eigenvalues are distinct (λ₁ ≠ λ₂ ≠ λ₃, the general case for orthorhombic and lower symmetry), then M_p must be a signed permutation matrix: each column has exactly one nonzero entry, which is ±1, and it can only swap axes with the same eigenvalue.

4. **Enumerate possible operations:** For distinct eigenvalues, the only operations are sign changes on each axis: (±1, ±1, ±1) → 8 operations (orthorhombic). For two equal eigenvalues (λ₁ = λ₂ ≠ λ₃, tetragonal), operations also include rotations in the degenerate subspace (the 1-2 plane): this adds 4-fold rotations → 16 operations. For all equal eigenvalues (λ₁ = λ₂ = λ₃, cubic), the full permutation group of three axes is allowed → 48 operations.

**The degenerate eigenvalue problem (cubic and higher symmetry):**

When eigenvalues are degenerate (equal), the eigenvectors are not unique—any rotation within the degenerate subspace gives equally valid eigenvectors. This means the "canonical frame" isn't canonical for high-symmetry lattices:
- **Cubic (3-fold degeneracy):** Any orthogonal frame is an eigenframe. The eigenvector decomposition provides no information about orientation.
- **Tetragonal (2-fold degeneracy):** The unique axis (c-axis, along the non-degenerate eigenvector) is well-defined, but the orientation within the basal plane is arbitrary.
- **Orthorhombic and lower (no degeneracy):** The eigenframe is unique (up to sign), giving a true canonical orientation.

**Near-degenerate case (pseudosymmetry):** For a slightly distorted cubic cell (λ₁ ≈ λ₂ ≈ λ₃ but not exactly equal), the eigenvalues are distinct so the eigenframe is technically unique—but the eigenvectors are numerically unstable (a tiny perturbation can rotate them by 45°). This is precisely the pseudosymmetry problem: the eigenframe "wants" to be degenerate but isn't quite.

**Potential value for Spacey.jl:**

- **Fast pre-screening:** The eigenvalues of G immediately tell you the *maximum possible* point group order: all distinct → ≤ 8, two equal → ≤ 16, all equal → ≤ 48. This can short-circuit the search early.
- **Degeneracy detection:** If two eigenvalues are within tolerance of each other, the lattice is near a higher-symmetry type. This could trigger a pseudosymmetry warning or an extended search.
- **Canonical orientation for non-degenerate cases:** For orthorhombic and lower symmetry (the majority of real crystals), the eigenframe provides a unique orientation without arbitrary choices. This could simplify the `snapToSymmetry_SVD` logic.
- **Minimum-strain symmetrization:** Recent work (Phys. Rev. Research 2, 013077, 2021) uses the eigenvalue structure of G to find the minimum-strain transformation to a higher-symmetry Bravais type. This is a more principled version of "snap to symmetry."

**Limitations:**
- Adds no value for cubic lattices (the most symmetric case) due to full degeneracy.
- The eigenvector instability near degeneracies makes it dangerous to rely on the eigenframe for lattices that are close to higher symmetry.
- The eigenvalue decomposition itself has O(n³) cost (negligible for 3×3) but introduces its own floating-point errors.

---

## 6. High Aspect Ratio: What Other Codes Do

High aspect ratio unit cells are a known challenge for all crystallographic codes. The literature distinguishes:

- **Intrinsic high aspect ratio:** the actual crystal has lattice vectors that differ by orders of magnitude (e.g., layered materials with van der Waals gaps). Spacey.jl's Minkowski-reduction approach handles this well up to aspect ratio ~500.
- **Artificial high aspect ratio:** the user provides a non-primitive supercell. Reducing to a primitive cell (another Minkowski reduction) eliminates this; the primitive cell always has an aspect ratio bounded by a dimension-dependent constant.

Spglib handles high aspect ratios by first finding the primitive cell. AFLOW-SYM does the same. VASP's hard-coded search bounds can fail for extreme cases. For intrinsically anisotropic structures, all codes eventually struggle and rely on user-supplied tolerance tuning.

**Layered materials (2D physics):** Codes like Phonopy that target 2D materials handle the "vacuum direction" (large lattice parameter out of the plane) by using a 2D symmetry analysis that ignores the vacuum direction. This is orthogonal to Spacey.jl's current scope but worth noting.

**Higher-precision arithmetic:** For extreme aspect ratios, `BigFloat` (arbitrary precision) or `Double64` (from `DoubleFloats.jl`, roughly 128-bit, ~2× slower) are viable options in Julia. The cost is moderate; acceptable if confined to a fallback path triggered only when `aspectRatio > 500`.

---

## 7. Open Problems and Algorithmic Frontiers

### 7.1 Pseudosymmetry

Many real crystal structures are "pseudosymmetric": they appear highly symmetric but have subtle distortions that break the ideal symmetry. For example, a slightly compressed cubic cell looks cubic but is actually tetragonal. Managing the tolerance so that:
- Genuine structural distortions are not misidentified as symmetry
- Measurement noise is not misidentified as broken symmetry

is an open challenge. AFLOW-SYM's self-consistent tolerance addresses this. Research groups using machine learning to distinguish pseudosymmetry from true distortions (e.g., Bayesian symmetry detection) represent a frontier approach.

### 7.2 Generalized Symmetry (Layer Groups, Rod Groups)

For 2D materials (graphene, MoS₂, etc.), the relevant symmetry is a **layer group** (80 groups, compared to 230 space groups). For 1D materials (nanotubes, nanowires), it is a **rod group**. Both are subgroups of the full 3D space group obtained by restricting certain translation directions. FINDSYM handles layer groups; no Julia package currently does.

### 7.3 Quasicrystal Symmetry

Quasicrystals have long-range order without periodicity, and their symmetry is described by **superspace groups** (higher-dimensional periodic structures projected to 3D). FINDSYM has been extended to handle these. This is a significant extension and likely out of scope for Spacey.jl in the near term.

### 7.4 Point Groups of Molecules vs. Lattices

Molecular point groups (C₂ᵥ, Td, Oₕ, etc.) are found by different algorithms than lattice point groups, because molecules have a finite number of atoms rather than an infinite periodic lattice. SOFI, Symmetrizer, and MSGSym are designed for this setting. Spacey.jl's algorithms are specific to the lattice case and exploit the periodicity structure.

---

## 8. Key References

| Reference | Relevance |
|-----------|-----------|
| Togo & Tanaka, *Sci. Technol. Adv. Mater.: Methods* 2024 | Spglib algorithm and tolerance strategy |
| Hicks et al., *Acta Cryst. A* 74 (2018) | AFLOW-SYM self-consistent tolerance |
| Stokes & Hatch, *J. Appl. Cryst.* 38 (2005) | FINDSYM space group identification |
| Křívý & Gruber, *Acta Cryst. A* 32 (1976) | Niggli reduction algorithm |
| Andrews & Bernstein, *J. Appl. Cryst.* 55 (2022) | Improved Niggli reduction numerical stability |
| Gunde et al., *J. Chem. Phys.* 161 (2024) | SOFI: shape-matching point groups |
| Largent et al., *J. Comput. Chem.* 33 (2012) | Symmetrizer: bottom-up symmetry elements |
| Pilati & Forni, *J. Appl. Cryst.* 31 (1998) | SYMMOL: maximum symmetry within tolerance |
| Grosse-Kunstleve & Adams, *Acta Cryst. A* 58 (2002) | Lattice symmetry from metric tensor |
| Espitau & Joux, arXiv:1905.11743 (2019) | Certified lattice reduction with interval arithmetic |
| Hart & Forcade, *Phys. Rev. B* 77 (2008) | Minkowski reduction and derivative structures |
| Lenstra, Lenstra & Lovász, *Math. Ann.* 261 (1982) | LLL lattice basis reduction |
| Phys. Rev. Research 2, 013077 (2021) | Minimum-strain symmetrization via metric tensor |

---

## 9. Summary of Actionable Ideas for Spacey.jl

| Idea | Source | Difficulty | Impact |
|------|--------|------------|--------|
| Self-consistent tolerance (`pointGroup_auto`) | AFLOW-SYM | Medium | High |
| Start-tight-then-loosen (avoid over-identification) | Spacey analysis | Low | High |
| Metric tensor for Bravais type classification | Literature | Low | Medium |
| Niggli reduction as cross-check | FINDSYM | Low | Medium |
| Extend neighbor search for degenerate cases | Spacey analysis | Low | Medium |
| `DoubleFloats.jl` fallback for aspect ratio > 500 | Julia ecosystem | Low | Medium |
| Rational arithmetic validation of integer ops | Sage/literature | Low | Medium |
| Interval arithmetic for Minkowski reduction | Espitau & Joux | Medium | Medium |
| Bottom-up element detection as cross-check | Symmetrizer/SYMMOL | Medium | Medium |
| Eigenvalue pre-screening of max group order | Standard linear algebra | Low | Low |
| Explicit identity/inverse checks in `isagroup` | Standard algebra | Low | Low |
| Space group identification (future) | Spglib/FINDSYM | High | High |
| Layer/rod group support (future) | FINDSYM | High | Medium |
