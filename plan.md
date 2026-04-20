# Robustness Improvement Plan for Spacey.jl

This document analyzes the current implementation of `pointGroup_robust` (and related functions) and suggests targeted improvements to increase numerical robustness, correctness guarantees, and test coverage. No changes are proposed here; this is a planning document only.

---

## 1. Deep Understanding of the Current Robustness Strategy

### Why the code is already strong

The core insight is that **Minkowski reduction makes the problem well-conditioned**. A Minkowski-reduced basis has the shortest possible vectors spanning the lattice, which minimises the ratio of the largest to smallest scale in the computation. Any symmetry operator of the lattice must permute the shortest lattice vectors among themselves, so searching only among the 27 combinations `A * [i, j, k]` (i,j,k ∈ {-1,0,1}) is provably sufficient for a reduced basis—this is the key theorem that makes the search finite and correct.

The layered filtering in `pointGroup_robust` then narrows a worst-case 27³ = 19,683 candidate bases down to a small set using three cheap tests (norm match → volume conservation → orthogonality), each removing the bulk of false candidates before the more expensive next step.

Volume normalization (`u, v, w ./= ∛|det A|`) before all comparisons brings the computation to a canonical scale, eliminating one source of floating-point magnitude variation.

Finally, choosing the **largest subset of candidates that closes under multiplication** (checking `[48, 24, 16, 12, 8, 4, 2]` in decreasing order) enforces the group axioms algebraically—a spurious near-miss that doesn't compose with anything else is automatically discarded.

### Known weak points

The following sections identify specific fragilities in the existing code.

---

## 2. Suggested Code Changes

### 2.1 Inconsistent tolerance in the candidate-distinctness check

**Status:** Accepted, implement.

**Location:** `pointGroup_robust`, line 235

```julia
A′ = [[i j k] for i ∈ c1 for j ∈ c2 if !(i≈j) for k ∈ c3 if !(i≈k) && !(j≈k)]
```

`!(i≈j)` uses Julia's default `isapprox` tolerance (1e-8 relative), completely ignoring the user-supplied `tol`. When `tol` is large (noisy inputs) this check is far stricter than everything else in the pipeline, so nearly-duplicate candidate vectors pass the norm-match filter but are then dropped as "same"—reducing the candidate pool too aggressively. Conversely, when `tol` is very small the default 1e-8 could be too loose.

**Fix:** Replace `!(i≈j)` with `!isapprox(i, j; rtol=tol)` throughout.

---

### 2.2 Suspicious extra factor in the volume tolerance

**Status:** Accepted, implement.

**Location:** `pointGroup_robust`, line 236

```julia
A′ = A′[findall([isapprox(abs(det(i)),vol,rtol=tol*min(norms...)) for i in A′])]
```

After volume normalization `vol ≈ 1` and `min(norms...) ≈ 1`, so the factor `tol*min(norms...)` is numerically `≈ tol`. However, for non-unit-volume inputs (which should not occur after the normalization, but might in edge cases) this formula mixes dimensionless `tol` with a dimensional length, producing an inconsistent relative tolerance. The intent appears to be just `rtol=tol`, and the extra factor should either be documented as deliberate (with a mathematical justification) or removed.

**Fix:** Replace with `isapprox(abs(det(i)), vol; rtol=tol)` and add a comment explaining why.

---

### 2.3 `pointGroup_robust` silently requires a pre-reduced basis but errors obscurely

**Status:** Noted for future consideration. No action now.

**Location:** `pointGroup_robust`, lines 216–218

```julia
if !(orthogonalityDefect(u,v,w)≈orthogonalityDefect(minkReduce(u,v,w)[1:3]...))
    error("Input basis for 'pointGroup' is not reduced. ...")
end
```

The error message says to "use `minkReduce`" but the public-facing wrapper `pointGroup(A)` already reduces; the danger is for direct callers of `pointGroup_robust`. This creates a footgun: a user who calls `pointGroup_robust` directly with an unreduced basis gets a hard error. The comparison `orthogonalityDefect(u,v,w)≈orthogonalityDefect(minkReduce(u,v,w)[1:3]...)` uses default tolerance, which can fail for degenerate cases where reduction is idempotent but floating-point jitter changes the defect slightly.

**Note:** `pointGroup_robust` is intended primarily for internal/testing use; typical users call `pointGroup(A)` which handles reduction automatically. This footgun is acceptable for now.

**Fix option A (if revisited):** Auto-reduce inside `pointGroup_robust` and emit a `@warn` rather than an error.

**Fix option B (if revisited):** Use a relative tolerance for the defect comparison, e.g. `rtol=1e-6`, to tolerate floating-point noise in the defect calculation itself.

---

### 2.4 Hardcoded 10% tolerance in `avgVecOverOps`

**Status:** No action. The 10% cutoff is very generous and is unlikely to be the source of any real failure.

**Location:** `avgVecOverOps`, line 15

```julia
cands = filter(i->norm(i-vec)<.1*norm(vec), cands)
```

The 10% cutoff is unrelated to any user-specified tolerance, but in practice it is so generous that it should never be too tight. This function is part of the `snapToSymmetry_avg` path, which is less critical than the main `pointGroup_robust` path. We can revisit this if `snapToSymmetry` is reworked later.

---

### 2.5 `pointGroup_robust` group-selection exits at first valid group size

**Status:** Accepted, implement.

**Location:** `pointGroup_robust`, lines 252–259

The code iterates `[48, 24, 16, 12, 8, 4, 2]` and returns the first (largest) size for which the sorted top-k candidates form a group. This is correct in spirit but has a subtle failure mode.

**The problem in detail:** Candidates are sorted by their deviation from orthogonality (`norm(T - I(3))`). If there are more candidates than the group size being tested (e.g., 50 candidates passing the orthogonality filter, testing for a 48-element group), only the top 48 are checked. A spurious candidate—one that passed the orthogonality check but isn't a real symmetry—could wedge itself into the top-48 by having a slightly smaller `norm(T - I(3))` than a legitimate operation. This bumps a legitimate operation out, so the top 48 don't form a group, and the algorithm falls through to the next smaller group size (24), returning an incorrect subgroup.

**When could this happen?** When `tol` is generous (the intended use case for noisy inputs), many near-miss candidates pass the orthogonality filter. If one of those near-misses happens to have a slightly better `T ≈ I` deviation than a true operation that was distorted by noise, the top-k selection is corrupted.

**Fix:** For each group size `il`, instead of testing only the top `il` candidates, test a slightly wider window. Specifically, try the top `il` first (fast path). If that fails, try replacing each of the bottom few candidates in the top-`il` set with candidates ranked just outside it (positions `il+1`, `il+2`, etc., up to `il+4`). This catches the case where one or two spurious candidates displace legitimate ones without introducing combinatorial explosion. The cost is at most ~4 extra `isagroup` calls per group size tested, which is negligible since `isagroup` on 48 integer matrices is fast (exact arithmetic, no floating point).

---

### 2.6 `isagroup` float tolerance is independent of `pointGroup_robust` tolerance

**Status:** Accepted, plan to implement.

**Location:** `isagroup` (float overload), line 73

```julia
function isagroup(...; atol=1e-8, rtol=1e-8)
```

Inside `pointGroup_robust`, `isagroup(ops[tp[1:il]])` is called with integer matrices so the exact-integer overload is used—this is fine. But `isagroup(G)` (the Cartesian/float form called in tests) uses fixed 1e-8 tolerances. If the Cartesian rotations were computed from a noisy basis, their products may differ by more than 1e-8 from exact group elements, making `isagroup(G)` return false even for a correct result.

**Fix:** Add a `tol` parameter to `snapToSymmetry_SVD` and pass a consistent tolerance to `isagroup` when validating the float operations.

---

### 2.7 Neighbor search radius may be insufficient in certain cases

**Status:** Accepted, plan to implement.

**Location:** `pointGroup_fast` and `pointGroup_robust`, lines 174/228

```julia
c = [A*[i,j,k] for i ∈ (-1,0,1) for j ∈ (-1,0,1) for k ∈ (-1,0,1)]
```

The proof that 27 neighbors suffice relies on the basis being *exactly* Minkowski-reduced. The search can fail when the Minkowski reduction itself loses precision.

**Failure modes beyond just high aspect ratio:**
- **High aspect ratio (known):** When vector lengths span many orders of magnitude, the floating-point Minkowski reduction may not satisfy the theoretical bounds tightly.
- **Near-degenerate bases:** When two or more basis vectors have nearly identical lengths, the reduction algorithm may not consistently pick the "right" shortest vector. Small perturbations could flip which vector is chosen, and the resulting basis might not satisfy the tight angular bounds needed for the ±1 coefficient proof.
- **Angles near the 60°/120° Minkowski bounds:** When interaxial angles are close to the limits (π/3 or 2π/3), floating-point jitter in the reduction could produce a basis that technically violates the bounds by a tiny amount, allowing a symmetry operation to map a vector to a ±2 coefficient combination.

**Exploratory fix:** For inputs where `aspectRatio > threshold` (e.g., 100) OR where angles between reduced basis vectors are within 1° of 60°/120°, extend the search to `[i,j,k]` for i,j,k ∈ {-2,-1,0,1,2} (125 candidates). The extra cost (125 vs 27 candidate vectors, leading to more candidate bases) is manageable and correctness is still verified by group-closure testing.

---

### 2.8 Auto-try tolerance sweep (self-consistent tolerance)

**Status:** Future consideration. No action now.

The current API requires the user to supply `tol`. For automated pipelines, a wrong `tol` can silently return a subgroup or a too-large group. Following the AFLOW-SYM philosophy (see `research.md`), the function should be able to determine a self-consistent tolerance.

**Proposed algorithm:**

1. Run with a generous starting tolerance (e.g., `tol_start = 0.15`).
2. Record the group size found: `n_found`.
3. Halve the tolerance and rerun. If the group size changes, the original was suspect.
4. Continue until the group size is stable over two consecutive tolerance halvings, or tolerance is below a minimum floor (e.g., `1e-10`).
5. Return the result from the finest tolerance that still gives a valid group.

This can be exposed as a new function `pointGroup_auto` rather than changing the existing API.

---

### 2.9 Use higher-precision arithmetic for extreme aspect ratios

**Status:** Future consideration. No action now.

For aspect ratios > 500, Float64 mantissa bits are exhausted in distinguishing the short and long axes. One option: detect this case and promote the computation to `BigFloat` or `Double64` (from the `DoubleFloats.jl` package). The snap-to-symmetry workflow followed by recomputation is already the recommended workaround, but it could fail when the initial point group is wrong due to precision loss.

---

### 2.10 Validate identity and inverses explicitly in `isagroup`

**Status:** Accepted, implement.

The `isagroup` function only checks closure (and distinctness). For lattice point groups all of {identity, inverses, closure} follow from being a finite closed set of linear maps with determinant ±1, but being explicit makes the validation stronger and catches bugs where a candidate passes closure but is not really orthogonal.

**Computational cost:** Negligible. Identity check is O(n)—one comparison per element. Inverse check for 3×3 integer matrices with det = ±1 is cheap: the inverse is the adjugate (cofactor matrix), which requires only a few multiplications per matrix. Membership check after computing each inverse is O(n). Total cost is O(n²)—the same order as the closure check already performed. Since n ≤ 48 for 3D lattice point groups, this adds at most 48² = 2304 cheap integer comparisons.

**Fix:** Add identity check (`I(3) ∈ members`) and inverse check (`∀g, g⁻¹ ∈ members`) to `isagroup`. For integer matrices, also check `det(g) == ±1` as a fast pre-filter.

---

### 2.11 Return a quality metric alongside the operators

**Status:** Accepted for later implementation.

Currently `pointGroup_robust` returns `(ops, Rc)` with no indication of how reliable the result is. Adding a confidence measure (e.g., the maximum deviation `max(norm(t - I(3)) for t ∈ T)`) would help callers decide whether to trust the result or increase the tolerance.

---

## 3. Suggested New Unit Tests

### 3.1 All 14 Bravais lattices × random orientations

**Status:** Implement.

The current "LG/G tests" check correctness of each lattice but only in canonical orientation. Add a test that applies 20+ random Euler angle triples to each of the 14 lattices and verifies the point group size is unchanged.

```julia
@testset "All Bravais lattices, random orientations" begin
    for (name, (A, nops)) in BravaisLatticeList
        for _ in 1:20
            α, β, γ = 2π*rand(), π*rand(), 2π*rand()
            u, v, w = threeDrotation(A[:,1], A[:,2], A[:,3], α, β, γ)
            u, v, w = minkReduce(u, v, w)[1:3]
            LG, G = pointGroup_robust(u, v, w)
            @test length(LG) == nops
            @test isagroup(LG)
            @test LG_G_test(LG, G, [u v w])
        end
    end
end
```

### 3.2 Structured noise patterns

**Status:** Implement (replaces original "noise scaling" suggestion).

The current random-noise test (in `runtests.jl` lines 164–184) already explores the full (tol, ε) parameter space thoroughly. Rather than duplicating that with a slightly different parameterization, this test adds value by using **structured (non-random) perturbations** that target specific failure modes:

- **Axis-aligned perturbations:** Perturb only one basis vector at a time, which can push the lattice toward a different Bravais type (e.g., perturbing one axis of cubic toward tetragonal).
- **Angular perturbations:** Rotate one basis vector slightly, breaking orthogonality without changing lengths.
- **Length-only perturbations:** Scale one vector without rotating it, testing the norm-matching filter specifically.

These structured perturbations can trigger failures that random noise (which averages out directionally) might not.

```julia
@testset "Structured noise patterns" begin
    for (name, (A, nops)) in BravaisLatticeList
        u, v, w = minkReduce(eachcol(A)...)[1:3]
        ε = 1e-4
        tol = 5e-3

        # Axis-aligned: perturb one vector only
        for δ in ([ε,0,0], [0,ε,0], [0,0,ε])
            LG, _ = pointGroup_robust((u .+ δ), v, w; tol=tol)
            @test length(LG) == nops
        end

        # Length-only: scale one vector by (1+ε)
        LG, _ = pointGroup_robust(u * (1+ε), v, w; tol=tol)
        @test length(LG) == nops

        # Angular: rotate one vector slightly around a random axis
        θ = ε
        R_small = I(3) + θ * [0 -1 0; 1 0 0; 0 0 0]  # small rotation about z
        LG, _ = pointGroup_robust(R_small * u, v, w; tol=tol)
        @test length(LG) == nops
    end
end
```

### 3.3 Group axioms: identity and inverses

**Status:** Implement.

```julia
@testset "Group axioms: identity and inverses" begin
    for (name, (A, nops)) in BravaisLatticeList
        A_red = minkReduce(A)
        LG, G = pointGroup(A_red)
        @test I(3) in LG                     # identity present
        @test all(inv(g) in G for g in G)    # closed under inverse (float)
        @test all(round.(Int, g) in LG for g in [round.(Int, inv(g)*1.0) for g in LG])
    end
end
```

### 3.4 Integer matrix properties

**Status:** Implement.

All integer point-group operations must have determinant ±1 and entries in {-2, -1, 0, 1, 2} (for 3D Minkowski-reduced bases).

```julia
@testset "Integer operation properties" begin
    for (name, (A, nops)) in BravaisLatticeList
        LG, _ = pointGroup(minkReduce(A))
        @test all(abs(det(g)) == 1 for g in LG)
        @test all(all(-2 .≤ g .≤ 2) for g in LG)  # entries bounded for reduced basis
    end
end
```

### 3.5 Cartesian operations are orthogonal

**Status:** Implement.

**Rationale:** The integer operations (lattice coordinates) are validated exactly by `isagroup`. But the Cartesian rotation matrices `G` are computed via a floating-point similarity transformation `R = A * M * inv(A)`. This conversion can accumulate error, especially for ill-conditioned `A`. Testing `RᵀR ≈ I` to tight tolerance (1e-12) validates that the representation conversion hasn't introduced significant error. Since point-group operations are by definition isometries, this property must hold exactly — any deviation indicates a numerical problem in the conversion, not an approximate symmetry.

```julia
@testset "Cartesian operations are orthogonal" begin
    for (name, (A, nops)) in BravaisLatticeList
        _, G = pointGroup(minkReduce(A))
        @test all(isapprox(g'*g, I(3); atol=1e-12) for g in G)
        @test all(isapprox(abs(det(g)), 1.0; atol=1e-12) for g in G)
    end
end
```

### 3.6 Snap-then-recompute gives exact group

**Status:** Implement.

After `snapToSymmetry_SVD`, the snapped basis should satisfy the point group with zero noise. Test that the resulting point group is *exactly* the expected size and forms an exact group (using tight tolerances).

```julia
@testset "Snap-then-recompute: exact result" begin
    A0 = [1.0 0 0; 0 1 0; 0 0 1]  # cubic
    for _ in 1:10
        ε = 5e-3
        Anoisy = A0 + (2*rand(3,3).-1)*ε
        u, v, w = minkReduce(eachcol(Anoisy)...)[1:3]
        ops, _ = pointGroup_robust(u, v, w; tol=5e-2)
        a, b, c, iops, rops = snapToSymmetry_SVD(u, v, w, ops)
        @test length(iops) == 48
        @test isagroup(iops)
        @test all(isapprox(g'*g, I(3); atol=1e-12) for g in rops)
    end
end
```

### 3.7 Non-primitive cell invariance

**Status:** Implement.

**Rationale:** The point group is an intrinsic property of the *lattice*, not of the particular basis used to describe it. A non-primitive (conventional) cell describes the same lattice as the primitive cell — Minkowski reduction should recover equivalent shortest vectors, and the point group should be identical. This tests that the full pipeline (Minkowski reduction → point group finding) is invariant to the choice of generating vectors, which is a non-trivial integration test: the Minkowski reduction must correctly "undo" the supercell to find the short vectors, and the point group finder must then work correctly on whatever basis the reduction produces.

Specific cases where this is non-trivial:
- A 2×1×1 supercell of a tetragonal cell has one vector twice the length of the other — Minkowski reduction must find the original short vector.
- A face-centered conventional cell has off-diagonal vectors — reduction must find the primitive vectors.

```julia
@testset "Non-primitive cell invariance" begin
    # Tetragonal: primitive vs 2×1×1 supercell
    A_tet = [1.0 0 0; 0 1 0; 0 0 1.5]
    for (sx, sy, sz) in [(2,1,1), (1,2,1), (2,2,1), (1,1,2)]
        A_sc = diagm([sx, sy, sz]) * A_tet
        LG, _ = pointGroup(minkReduce(A_sc))
        @test length(LG) == 16
    end
    # Cubic: primitive vs supercells
    A_cub = [1.0 0 0; 0 1 0; 0 0 1]
    for (sx, sy, sz) in [(2,1,1), (2,2,1), (2,2,2), (3,3,3)]
        A_sc = diagm([sx*1.0, sy*1.0, sz*1.0]) * A_cub
        LG, _ = pointGroup(minkReduce(A_sc))
        @test length(LG) == 48
    end
end
```

### 3.8 Pseudosymmetry: over-identification of symmetry

**Status:** Implement.

**Rationale:** In real applications (DFT relaxations, experimental refinements), structures often emerge that are *almost* but not quite a higher-symmetry Bravais type. For example, a tetragonal structure with c/a = 1.001 is very close to cubic. The most common and dangerous failure mode for symmetry finders is **over-identification**: reporting higher symmetry than actually present. This matters because:
- Downstream calculations (phonons, electronic bands) use the symmetry to reduce computational effort. Wrong symmetry → wrong physics.
- The boundary between "noise on a cubic" and "genuinely tetragonal" is exactly the kind of ambiguity where tolerance settings cause silent failures.

This test verifies that the algorithm does **not** over-identify: a small but definite tetragonal distortion (c/a = 1.001 up to 1.1) must always be reported as tetragonal (16 ops), not cubic (48 ops), when the tolerance is set small enough to resolve the distortion.

```julia
@testset "Pseudosymmetry: cubic vs tetragonal" begin
    for c_ratio in [1.001, 1.01, 1.05, 1.1]
        A = [1.0 0 0; 0 1 0; 0 0 c_ratio]
        u, v, w = minkReduce(eachcol(A)...)[1:3]
        # tol must be smaller than the distortion to resolve it
        LG, _ = pointGroup_robust(u, v, w; tol=1e-3)
        @test length(LG) == 16  # tetragonal, not cubic
    end
    # Analogous: orthorhombic near tetragonal
    for b_ratio in [1.001, 1.01, 1.05]
        A = [1.0 0 0; 0 b_ratio 0; 0 0 1.5]
        u, v, w = minkReduce(eachcol(A)...)[1:3]
        LG, _ = pointGroup_robust(u, v, w; tol=1e-3)
        @test length(LG) == 8  # orthorhombic, not tetragonal
    end
end
```

### 3.9 Aspect ratio stress test (parametric)

**Status:** Implement.

**Rationale:** The current tests check a few specific high-aspect-ratio cases (256, 500, 512, 1024) with hard-coded noise. This test systematically maps the working envelope of the algorithm as a function of aspect ratio alone (no noise), documenting exactly where the algorithm transitions from reliable to unreliable. This serves two purposes: (1) it defines the guarantee we can make to users (currently: reliable up to ~500), and (2) if we improve robustness (e.g., via suggestion 2.7), this test will immediately show whether the working envelope has expanded.

```julia
@testset "Aspect ratio stress test" begin
    results = Dict()
    for log2_ar in 1:12  # aspect ratio 2 to 4096
        ar = 2^log2_ar
        a1 = [1.0, 0, 0]; a2 = [0, ar*1.0, 0]; a3 = [0, 0, 1.0]
        u, v, w = minkReduce(a1, a2, a3)[1:3]
        LG, _ = pointGroup_robust(u, v, w; tol=0.1)
        results[ar] = length(LG)
    end
    # Document the results
    println("Aspect ratio working envelope:")
    for (ar, n) in sort(collect(results))
        println("  AR=$ar: $n ops (expected 16)")
    end
    # These should work:
    @test results[2]   == 16
    @test results[4]   == 16
    @test results[16]  == 16
    @test results[256] == 16
end
```

### 3.10 Consistency across function variants

**Status:** Implement.

All three functions (`pointGroup_simple`, `pointGroup_fast`, `pointGroup_robust`) should give the same group size for exact (no-noise) inputs. Note: `pointGroup_simple` and `pointGroup_fast` will eventually be removed once `pointGroup_robust` is fully validated, but until then this cross-check is a useful safety net.

```julia
@testset "Consistency across variants" begin
    for (name, (A, nops)) in BravaisLatticeList
        u, v, w = minkReduce(eachcol(A)...)[1:3]
        @test length(pointGroup_simple(u, v, w)) == nops
        @test length(pointGroup_fast(u, v, w))   == nops
        @test length(pointGroup_robust(u, v, w)[1]) == nops
    end
end
```

### 3.11 Tolerance edge cases: straddling the boundary

**Status:** Implement.

Verify that at the exact threshold (where tolerance just barely admits correct operations), the algorithm returns the right group and not a subgroup or supergroup.

```julia
@testset "Tolerance boundary behavior" begin
    # Place noise at exactly 10% of tolerance; should succeed
    ε = 1e-4
    tol = 10 * ε
    A = [1.0 0 0; 0 1 0; 0 0 1]
    for _ in 1:50
        Atest = hcat(minkReduce(eachcol(A + (2*rand(3,3).-1)*ε)...)[1:3]...)
        LG, _ = pointGroup(Atest; tol=tol)
        @test length(LG) == 48
    end
    # Place noise at twice the tolerance; should fail (group size < 48)
    ε_big = 0.5 * tol
    failed = 0
    for _ in 1:50
        Atest = hcat(minkReduce(eachcol(A + (2*rand(3,3).-1)*ε_big)...)[1:3]...)
        LG, _ = pointGroup(Atest; tol=tol)
        if length(LG) != 48; failed += 1; end
    end
    @test failed > 0  # At least some trials should fail with noise > tol
end
```

### 3.12 `snapToSymmetry_SVD` volume conservation (all Bravais types)

**Status:** Implement.

Currently only tested for cubic/tetragonal/orthorhombic. Extend to all 14 lattices.

```julia
@testset "Snap volume conservation, all Bravais types" begin
    for (name, (A, nops)) in BravaisLatticeList
        ε = 2e-3
        Anoisy = A + (2*rand(3,3).-1)*ε
        u, v, w = minkReduce(eachcol(Anoisy)...)[1:3]
        ops, _ = pointGroup_robust(u, v, w; tol=0.05)
        if length(ops) == nops  # Only test if we got the right group
            a, b, c, iops, rops = snapToSymmetry_SVD(u, v, w, ops)
            @test isapprox(det([a b c]), det([u v w]); rtol=1e-10)
            @test length(iops) == nops
        end
    end
end
```

---

## 4. Architectural Observations

These are noted for future consideration but no action is planned now.

- The `spacegroup` function stub currently returns `true` unconditionally. This should either be removed or clearly marked as `TODO` to avoid confusion.
- `Crystal` struct is defined but never used by any exported function. If space group work is planned, this struct will be the basis—keeping it but documenting it as forward-looking would be helpful.
- The `2D_snap_example_*.jl` files in `src/` are not part of the module. Consider moving them to `examples/` to avoid confusion with source code.
- The timing tests in `test/runTimingTests.jl` are excluded from CI. A lightweight version that just checks the _ratio_ (not absolute time) could be made CI-compatible.
