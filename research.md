# Literature Review: Point Group and Space Group Symmetry Detection

This document surveys algorithms, software packages, and mathematical results relevant to the crystallographic point-group-finding problem that Spacey.jl solves. The goal is to surface ideas that could make Spacey.jl more robust and to lay groundwork for future space-group functionality.

---

## 1. The Core Mathematical Problem

Given three basis vectors **a₁, a₂, a₃** ∈ ℝ³ defining a Bravais lattice Λ, find all linear maps R : ℝ³ → ℝ³ such that R(Λ) = Λ and R is an isometry (i.e., ‖Rv‖ = ‖v‖ for all v). The set of all such R forms the **point group** G(Λ), a finite subgroup of O(3).

Equivalently, in the basis matrix A = [**a₁** | **a₂** | **a₃**], find all integer matrices M ∈ GL(3, ℤ) with det(M) = ±1 such that AM Aᵀ = Aᵀ, i.e., AᵀA commutes with M in a specific way. Writing this out: AMᵀ(AᵀA)⁻¹Mᵀ⁻¹ = I.

The metric tensor **G = AᵀA** encodes everything; two bases give the same Bravais lattice if and only if they have the same metric tensor (up to a GL(3,ℤ) transformation). The point group in lattice coordinates is exactly the automorphism group of the metric tensor:

> **Aut(G) = { M ∈ GL(3,ℤ) : MᵀGM = G }**

---

## 2. Why Lattice Reduction Is Fundamental

Every major modern algorithm for lattice symmetry starts with lattice reduction, because the bounds that make the search finite and computationally tractable depend on it.

### 2.1 Minkowski Reduction

A basis {**b₁**, **b₂**, **b₃**} is Minkowski-reduced if each **bᵢ** is the shortest vector in the lattice that can extend {**b₁**, ..., **bᵢ₋₁**} to a basis. This is the reduction used by Spacey.jl (via the `MinkowskiReduction.jl` package).

Key property for point-group algorithms: a Minkowski-reduced basis satisfies tight geometric bounds on the angles and length ratios. Specifically, for any two basis vectors: cos(θᵢⱼ) ≤ 1/2, meaning all angles between basis vectors lie in [60°, 120°]. This is exactly the condition needed to prove that any symmetry operation maps basis vectors to ±1 linear combinations of themselves—the 27-neighbor search is sufficient.

**Computational note:** Minkowski reduction is NP-hard in general dimension but polynomial (O(n³)) in fixed dimension 3, and practical implementations (including `MinkowskiReduction.jl`) are fast. However, floating-point implementations lose precision when vector lengths span many orders of magnitude—this is the source of Spacey.jl's known failure at extreme aspect ratios (~2²⁵).

### 2.2 Niggli Reduction

The **Niggli reduced cell** (Niggli 1928; Křívý & Gruber 1976) provides a *canonical* representative for every lattice—the unique reduced form under a precise set of main conditions and special conditions. While Minkowski reduction gives the *shortest* basis, Niggli reduction gives the *unique* description.

The key advantage for symmetry work: once a Niggli cell is computed, the Bravais lattice type can be determined by a table lookup using the metric tensor components, and the point group follows from the lattice type. Modern implementations (e.g., Andrews & Bernstein 2022, *J. Appl. Cryst.*) have improved the numerical stability of Niggli reduction, particularly for handling near-degenerate cases where multiple Niggli forms are nearby in parameter space.

**Relevance to Spacey.jl:** Niggli reduction is an alternative to Minkowski reduction. It could be used as a preprocessing step or as a cross-check. The NIST table-lookup approach after Niggli reduction is orthogonal to Spacey.jl's exhaustive-search approach and could serve as an independent validation.

### 2.3 LLL Reduction

The Lenstra–Lenstra–Lovász (LLL) algorithm is a polynomial-time approximation to Minkowski reduction. It does not find the theoretically shortest basis but is significantly faster for high-dimensional lattices and has well-understood numerical behavior. For 3D crystallography, true Minkowski reduction is both tractable and preferable; LLL is primarily relevant if Spacey.jl is ever extended to higher dimensions.

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

**Tolerance strategy:** Spglib uses an *adaptive* tolerance. When a symmetry check fails, the tolerance is iteratively adjusted (usually decreased, sometimes both increased and decreased) until the found symmetries satisfy all crystallographic constraints, or until a minimum tolerance floor is reached. This is more sophisticated than Spacey.jl's single-tolerance approach.

**Key difference from Spacey.jl:** Spglib works on *atomic structures* (positions + species) rather than bare lattices. It finds the space group, of which the point group is derived. For purely lattice point groups, Spglib can still be used but requires the atomic basis. Spacey.jl's advantage is operating purely on the lattice (no atomic basis required) and its explicit use of the Minkowski theorem to limit the search.

**Numerical robustness:** Spglib's adaptive tolerance is its main robustness mechanism. The iterative adjustment can take many iterations for near-symmetric structures and may not converge uniquely. The 2024 paper documents several historical bugs related to tolerance handling.

### 3.2 AFLOW-SYM

**Language:** C++ (part of the AFLOW framework)  
**Reference:** Hicks et al., *Acta Cryst. A* 2018; arXiv:1802.07977

AFLOW-SYM's primary contribution is **self-consistent tolerance determination**:

> Rather than asking the user to supply a tolerance, AFLOW-SYM searches for a tolerance ε such that the set of symmetry operations found at ε is consistent with a valid group (closed under composition, contains identity and inverses).

The algorithm:
1. Begins with a coarse tolerance and finds candidate operations.
2. Refines the tolerance until the candidates form a valid group.
3. Reports the "mapping tolerance" as a measure of how symmetric the structure is.

This is exactly the concept suggested in `plan.md §2.8` as `pointGroup_auto`. Implementing this for Spacey.jl would bring it to the state of the art.

**Additional AFLOW-SYM innovations:**
- **Symmetrization of input:** after finding approximate operations, the input structure is snapped to perfect symmetry (analogous to `snapToSymmetry_SVD`), and the tolerance is recomputed—an iterative refinement loop.
- **Validation at scale:** tested against 54,000+ ICSD entries and 1.7M AFLOW structures, finding that many structures in databases had incorrect symmetry assignments.

### 3.3 FINDSYM (ISOTROPY Software Suite)

**Reference:** Stokes & Hatch, *J. Appl. Cryst.* 2005

FINDSYM identifies space groups from crystal structures and transforms them to a standard setting. Unlike Spglib/AFLOW, FINDSYM provides symmetry in a standardized basis (ITA convention). Recently extended to handle **superspace groups** for modulated/quasicrystalline structures.

**Key algorithmic idea:** FINDSYM uses the Niggli cell for reduction, then a table lookup for the Bravais type, then an exhaustive search over the corresponding point-group operations. The table-lookup approach is fast and deterministic, though less flexible than Spacey.jl's metric-tensor-based approach.

### 3.4 SOFI: Point Groups by Shape Matching

**Reference:** Pangerl et al., *J. Chem. Phys.* 2024

SOFI (Symmetry Operations FInder) formulates point-group detection as a **degenerate shape-matching problem**. Given a set of points (atomic positions), it finds rotation/reflection operations that map the set to itself.

The key innovation: SOFI uses the **inertia tensor** of the point set to define a coordinate frame, then searches for symmetry operations within that frame. This is robust to numerical noise because the inertia tensor is a smooth function of the positions.

**Relevance:** SOFI is designed for molecular point groups (no periodicity), but the shape-matching formulation is interesting. For lattices, the "points" are the lattice vectors and their images; the inertia-tensor approach could provide an alternative method for defining a canonical frame before the search.

### 3.5 Symmetrizer

**Reference:** Pilati & Forni, *J. Comput. Chem.* 2012

Symmetrizer finds symmetry elements (rotation axes, mirror planes, inversion centers) by a **bottom-up** approach: detect candidate symmetry axes from the geometry, then compose them to determine the full group. This is the approach used in molecular symmetry programs like SYMMOL.

For lattices, the three shortest lattice vectors define natural candidate axes/planes. This bottom-up approach could complement Spacey.jl's top-down (exhaustive candidate + filter) approach.

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

### 5.2 Metric Tensor Approach

The metric tensor **G = AᵀA** is a positive-definite 3×3 symmetric matrix. Point-group operations are exactly the integer matrices M with MᵀGM = G. Comparing metric tensors is more numerically stable than comparing basis matrices because:
- The metric tensor is independent of the orientation of the basis (any rotation of A gives the same G up to a congruence transformation).
- The comparison MᵀGM ≈ G avoids the need to compute `inv(A)`, which amplifies errors for ill-conditioned bases.

Spacey.jl currently computes `AᵢAᵢᵀ` (related to the inverse metric) rather than `AᵀA` directly. Working with the metric tensor in forward form (not inverse) could improve conditioning for nearly-singular bases.

### 5.3 Exact Arithmetic with Integer Coefficients

After Minkowski reduction, all point-group operations have integer matrix representations. The group closure test `isagroup(LG)` in Spacey.jl already exploits this by using exact integer arithmetic. This is the strongest possible validation.

For the *finding* step (not just validation), exact arithmetic is possible if all computations are done in terms of the metric tensor **G** (rational entries) rather than the Cartesian basis (floating-point). Specifically:
- `G = AᵀA` has entries that are dot products of basis vectors.
- The condition MᵀGM = G can be checked with rational arithmetic if G is rational.
- For practical inputs (floating-point basis vectors), G is floating-point, but the comparison MᵀGM ≈ G can be made tolerant at a single step, after which everything is integer.

This is the **rational lattice** approach and is used in some mathematical software (e.g., Sage's `Lattice.automorphisms()`). It is more complex to implement but eliminates accumulated rounding error.

### 5.4 Interval Arithmetic

In principle, replacing all floating-point operations with **interval arithmetic** (e.g., using `IntervalArithmetic.jl` in Julia) would give guaranteed results: the true answer is guaranteed to lie within the computed interval. For crystallographic computations, the intervals would need to be very tight (the basis vectors are measured or computed to many significant figures), making this approach practical only for very well-conditioned problems. No major crystallographic package uses this approach.

### 5.5 Eigenvalue-Based Canonical Frames

The **eigenvectors of the metric tensor G** define a canonical coordinate frame for the lattice (independent of the input basis orientation). Expressing the symmetry search in this frame could improve robustness:
1. Diagonalize G = VDVᵀ.
2. Work in the eigenbasis, where G is diagonal (a "spherical" metric).
3. All symmetry operations in this frame are permutations of eigenvectors (with possible sign changes).

This approach is related to how SOFI uses the inertia tensor. It is exact for lattices with distinct eigenvalues; for cubic lattices (degenerate eigenvalues), an additional step is needed.

---

## 6. High Aspect Ratio: What Other Codes Do

High aspect ratio unit cells are a known challenge for all crystallographic codes. The literature distinguishes:

- **Intrinsic high aspect ratio:** the actual crystal has lattice vectors that differ by orders of magnitude (e.g., layered materials with van der Waals gaps). Spacey.jl's Minkowski-reduction approach handles this well up to aspect ratio ~500.
- **Artificial high aspect ratio:** the user provides a non-primitive supercell. Reducing to a primitive cell (another Minkowski reduction) eliminates this; the primitive cell always has an aspect ratio bounded by a dimension-dependent constant.

Spglib handles high aspect ratios by first finding the primitive cell. AFLOW-SYM does the same. For intrinsically anisotropic structures, both codes can struggle at extreme aspect ratios and rely on user-supplied tolerance tuning.

**Layered materials (2D physics):** Codes like Phonopy that target 2D materials handle the "vacuum direction" (large lattice parameter out of the plane) by using a 2D symmetry analysis that ignores the vacuum direction. This is orthogonal to Spacey.jl's current scope but worth noting.

**Higher-precision arithmetic:** For extreme aspect ratios, `BigFloat` (arbitrary precision) or `DoubleFloats.jl` (roughly 128-bit, hardware-accelerated on modern CPUs) are viable options in Julia. The cost is 2–10× slower; acceptable if confined to a fallback path triggered only when `aspectRatio > 500`.

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
| Pangerl et al., *J. Chem. Phys.* 2024 | SOFI: shape-matching point groups |
| Pilati & Forni, *J. Comput. Chem.* 33 (2012) | Symmetrizer: bottom-up symmetry elements |
| Grosse-Kunstleve & Adams, *Acta Cryst. A* 58 (2002) | Lattice symmetry from metric tensor |
| Lenstra, Lenstra & Lovász, *Math. Ann.* 261 (1982) | LLL lattice basis reduction |
| Hart & Forcade, *Phys. Rev. B* 77 (2008) | Minkowski reduction and derivative structures |

---

## 9. Summary of Actionable Ideas for Spacey.jl

| Idea | Source | Difficulty | Impact |
|------|--------|------------|--------|
| Self-consistent tolerance (`pointGroup_auto`) | AFLOW-SYM | Medium | High |
| Adaptive tolerance adjustment (iterative) | Spglib | Medium | High |
| Metric tensor formulation (avoid `inv(A)`) | Literature | Low | Medium |
| Niggli reduction as cross-check | FINDSYM | Low | Medium |
| Extend neighbor search for high aspect ratio | Spacey analysis | Low | Medium |
| `BigFloat` fallback for aspect ratio > 500 | Julia ecosystem | Low | Medium |
| Explicit identity/inverse checks in `isagroup` | Standard algebra | Low | Low |
| Eigenvalue-based canonical frame | SOFI/literature | High | Medium |
| Space group identification (future) | Spglib/FINDSYM | High | High |
| Layer/rod group support (future) | FINDSYM | High | Medium |
