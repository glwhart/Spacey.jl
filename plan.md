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

**Location:** `pointGroup_robust`, line 235

```julia
A′ = [[i j k] for i ∈ c1 for j ∈ c2 if !(i≈j) for k ∈ c3 if !(i≈k) && !(j≈k)]
```

`!(i≈j)` uses Julia's default `isapprox` tolerance (1e-8 relative), completely ignoring the user-supplied `tol`. When `tol` is large (noisy inputs) this check is far stricter than everything else in the pipeline, so nearly-duplicate candidate vectors pass the norm-match filter but are then dropped as "same"—reducing the candidate pool too aggressively. Conversely, when `tol` is very small the default 1e-8 could be too loose.

**Fix:** Replace `!(i≈j)` with `!isapprox(i, j; rtol=tol)` throughout.

---

### 2.2 Suspicious extra factor in the volume tolerance

**Location:** `pointGroup_robust`, line 236

```julia
A′ = A′[findall([isapprox(abs(det(i)),vol,rtol=tol*min(norms...)) for i in A′])]
```

After volume normalization `vol ≈ 1` and `min(norms...) ≈ 1`, so the factor `tol*min(norms...)` is numerically `≈ tol`. However, for non-unit-volume inputs (which should not occur after the normalization, but might in edge cases) this formula mixes dimensionless `tol` with a dimensional length, producing an inconsistent relative tolerance. The intent appears to be just `rtol=tol`, and the extra factor should either be documented as deliberate (with a mathematical justification) or removed.

**Fix:** Replace with `isapprox(abs(det(i)), vol; rtol=tol)` and add a comment explaining why.

---

### 2.3 `pointGroup_robust` silently requires a pre-reduced basis but errors obscurely

**Location:** `pointGroup_robust`, lines 216–218

```julia
if !(orthogonalityDefect(u,v,w)≈orthogonalityDefect(minkReduce(u,v,w)[1:3]...))
    error("Input basis for 'pointGroup' is not reduced. ...")
end
```

The error message says to "use `minkReduce`" but the public-facing wrapper `pointGroup(A)` already reduces; the danger is for direct callers of `pointGroup_robust`. This creates a footgun: a user who calls `pointGroup_robust` directly with an unreduced basis gets a hard error. The comparison `orthogonalityDefect(u,v,w)≈orthogonalityDefect(minkReduce(u,v,w)[1:3]...)` uses default tolerance, which can fail for degenerate cases where reduction is idempotent but floating-point jitter changes the defect slightly.

**Fix option A:** Auto-reduce inside `pointGroup_robust` and emit a `@warn` rather than an error.

**Fix option B:** Use a relative tolerance for the defect comparison, e.g. `rtol=1e-6`, to tolerate floating-point noise in the defect calculation itself.

---

### 2.4 Hardcoded 10% tolerance in `avgVecOverOps`

**Location:** `avgVecOverOps`, line 15

```julia
cands = filter(i->norm(i-vec)<.1*norm(vec), cands)
```

The 10% cutoff is unrelated to any user-specified tolerance, so `snapToSymmetry_avg` can fail silently for inputs with small perturbations that exceed 10% (e.g., lattice vectors rotated slightly under a lower-symmetry candidate operation). 

**Fix:** Pass `tol` as an argument and use it for this filter, defaulting to the current 0.1.

---

### 2.5 `pointGroup_robust` group-selection exits at first valid group size

**Location:** `pointGroup_robust`, lines 252–259

The code iterates `[48, 24, 16, 12, 8, 4, 2]` and returns the first (largest) size for which the sorted top-k candidates form a group. This is correct in spirit but has a subtle failure mode: if there are 50 candidates passing the orthogonality test (e.g., because `tol` is large), only the top 48 are tested for 48-element groups. A legitimate 48-element group might be the "second-best 48" if one spurious candidate wedges itself into rank 1–48.

**Fix:** For each group size `il`, test all `length(idx) choose il` subsets of the top `min(length(idx), 2*il)` candidates, not just the first `il`. Practically, since `il` is always a divisor of 48 and the spurious candidates rank poorly in the orthogonality sort, testing the top `il + k` candidates for small `k` (say `k = 4`) would catch most cases without combinatorial explosion.

---

### 2.6 `isagroup` float tolerance is independent of `pointGroup_robust` tolerance

**Location:** `isagroup` (float overload), line 73

```julia
function isagroup(...; atol=1e-8, rtol=1e-8)
```

Inside `pointGroup_robust`, `isagroup(ops[tp[1:il]])` is called with integer matrices so the exact-integer overload is used—this is fine. But `isagroup(G)` (the Cartesian/float form called in tests) uses fixed 1e-8 tolerances. If the Cartesian rotations were computed from a noisy basis, their products may differ by more than 1e-8 from exact group elements, making `isagroup(G)` return false even for a correct result.

**Fix:** Add a `tol` parameter to `snapToSymmetry_SVD` and pass a consistent tolerance to `isagroup` when validating the float operations.

---

### 2.7 Neighbor search radius may be insufficient for extreme aspect ratios

**Location:** `pointGroup_fast` and `pointGroup_robust`, lines 174/228

```julia
c = [A*[i,j,k] for i ∈ (-1,0,1) for j ∈ (-1,0,1) for k ∈ (-1,0,1)]
```

The proof that 27 neighbors suffice relies on the basis being Minkowski-reduced. However, when the Minkowski reduction itself loses precision (which happens around aspect ratio 2^25 for Float64), the "reduced" basis may not satisfy the theoretical bounds tightly. In those cases a symmetry operation could map a basis vector to a linear combination outside the ±1 range.

**Exploratory fix:** For inputs where `aspectRatio > threshold` (e.g., 100), also search `[i,j,k]` for i,j,k ∈ {-2,-1,0,1,2} (125 candidates), then verify correctness by group-closure testing. The extra cost is manageable because high-aspect-ratio cases are already slow.

---

### 2.8 Auto-try tolerance sweep (self-consistent tolerance)

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

For aspect ratios > 500, Float64 mantissa bits are exhausted in distinguishing the short and long axes. One option: detect this case and promote the computation to `BigFloat` or `Double64` (from the `DoubleFloats.jl` package). The snap-to-symmetry workflow followed by recomputation is already the recommended workaround, but it could fail when the initial point group is wrong due to precision loss.

---

### 2.10 Validate identity and inverses explicitly

The `isagroup` function only checks closure (and distinctness). For lattice point groups all of {identity, inverses, closure} follow from being a finite closed set of linear maps with determinant ±1, but being explicit makes the validation stronger and catches bugs where a candidate passes closure but is not really orthogonal.

**Fix:** Add identity check (`I(3) ∈ members`) and inverse check (`∀g, g⁻¹ ∈ members`) to `isagroup`. For integer matrices, `det(g) == ±1` is also cheap to check.

---

### 2.11 Return a quality metric alongside the operators

Currently `pointGroup_robust` returns `(ops, Rc)` with no indication of how reliable the result is. Adding a confidence measure (e.g., the maximum deviation `max(norm(t - I(3)) for t ∈ T)`) would help callers decide whether to trust the result or increase the tolerance.

---

## 3. Suggested New Unit Tests

### 3.1 All 14 Bravais lattices × random orientations

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

### 3.2 Noise scaling: verify tol ≫ noise implies correct result

For each Bravais lattice, add noise of magnitude `ε` and verify correctness at `tol = 10ε`. This is more systematic than the current random-noise test because it explicitly validates the `tol ≥ 12ε` heuristic documented in the code.

```julia
@testset "Noise scaling correctness" begin
    for (name, (A, nops)) in BravaisLatticeList
        for log_eps in -10:-1
            ε = 10.0^log_eps
            tol = 20 * ε
            tol > 0.3 && continue  # skip unreasonably large tolerance
            for _ in 1:30
                Atest = hcat(minkReduce(eachcol(A + (2*rand(3,3).-1)*ε)...)[1:3]...)
                LG, G = pointGroup(Atest; tol=tol)
                @test length(LG) == nops
                @test isagroup(LG)
            end
        end
    end
end
```

### 3.3 Group axioms: identity and inverses

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

### 3.7 Supercell invariance

A 2×2×2 supercell of a simple cubic lattice should have the same point group as the original.

```julia
@testset "Supercell invariance" begin
    for scale in [2, 3]
        A_sc = [scale 0 0; 0 scale 0; 0 0 scale] * 1.0
        LG, _ = pointGroup(minkReduce(A_sc))
        @test length(LG) == 48
    end
    # Non-cubic supercell
    A_tet = [1.0 0 0; 0 1 0; 0 0 1.5]
    for (sx, sy, sz) in [(2,2,1), (3,3,1), (2,2,2)]
        A_sc = [sx 0 0; 0 sy 0; 0 0 sz] .* A_tet
        LG, _ = pointGroup(minkReduce(A_sc))
        @test length(LG) == 16
    end
end
```

### 3.8 Near-degenerate (pseudosymmetry) cases

Test that a slightly distorted cubic lattice (ε-away from cubic) is correctly identified as tetragonal, not cubic, and that the threshold between these is sharp.

```julia
@testset "Pseudosymmetry: cubic vs tetragonal" begin
    for c_ratio in [1.001, 1.01, 1.05, 1.1]
        A = [1.0 0 0; 0 1 0; 0 0 c_ratio]
        u, v, w = minkReduce(eachcol(A)...)[1:3]
        LG, _ = pointGroup_robust(u, v, w; tol=1e-3)
        # c_ratio != 1 → tetragonal, not cubic
        @test length(LG) == 16
    end
end
```

### 3.9 Aspect ratio stress test (parametric)

Systematically test aspect ratios from 2 to the known failure threshold (~512), verifying the correct group is returned at each ratio. This documents the actual working envelope better than the current binary pass/fail tests.

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
    # Document the results (not a hard @test but printed)
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

All three functions (`pointGroup_simple`, `pointGroup_fast`, `pointGroup_robust`) should give the same group size for exact (no-noise) inputs.

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

- The `spacegroup` function stub currently returns `true` unconditionally. This should either be removed or clearly marked as `TODO` to avoid confusion.
- `Crystal` struct is defined but never used by any exported function. If space group work is planned, this struct will be the basis—keeping it but documenting it as forward-looking would be helpful.
- The `2D_snap_example_*.jl` files in `src/` are not part of the module. Consider moving them to `examples/` to avoid confusion with source code.
- The timing tests in `test/runTimingTests.jl` are excluded from CI. A lightweight version that just checks the _ratio_ (not absolute time) could be made CI-compatible.
