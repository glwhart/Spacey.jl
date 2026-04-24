using Spacey
using Test
using LinearAlgebra
using MinkowskiReduction
    
# Timing tests are not consistent on github virtual machines.
# Timing tests are only checked locally, see runTimingTests.jl

# This list is used in multiple tests. Additional lattices could be added. (At least) One instance of each of the 14 Bravais lattices is included.
BravaisLatticeList=Dict([("Centered monoclinic 1",     ([1.0  1.0 0.5; 1.1 -1.1 0.0; 0.0 0.0 0.7],4)),
("Simple monoclinic 1",      ( [1.0  0.0 0.1; 0.0  1.1 0.0; 0.0 0.0 0.7],4)),
("Base-centered orthorhombic 1",([1.0  0.0 0.5; 0.0  1.1 0.0; 0.0 0.0 0.7] ,8)),
("Triclinic 1",([1.0  0.1 0.2; 0.2  1.1 0.0; 0.3 0.0 0.7] ,2)),
("FCC unreduced 1",([0.0 0.5 1.0; 0.5 0.0 1.0; 0.5 0.5 1.0],48)),
("Rhombohedral 1",([1.0  1.0 .5; 1.0  .5 1.0; .5 1.0 1.0],12)),
("hexagonal 1",([1.0 0.5 0.0; 0.0  √(.75) 0.0; 0.0 0.0 1.6],24)),
("Base-Centered orthorhombic 2",([1.1 1.9 0.0; -1.1 1.9 0.0; 0.0 0.0 1.3],8)),
("Body-centered orthorhombic 1",([1.1 0.0 0.55; 0.0 1.9 0.95; 0.0 0.0 0.7],8)),
("Face-centered orthorhombic 1",([0.55 0.0 0.55; 0.95 0.95 0.0; 0.0 0.35 0.35],8)),
("BCTet",([0.0 0.5 0.5; 0.5 0.0 0.5; 0.54 0.54 0.0],16)),
("FCC",([0.0 0.5 0.5; 0.5 0.0 0.5; 0.5 0.5 0.0],48)),
("BCC",([-1.0 1.0 1.0; 1.0 -1.0 1.0; 1.0 1.0 -1.0],48)),
("Simple cubic",([1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0],48)),
("Simple orthorhombic",([1.1 0.0 0.0; 0.0 1.2 0.0; 0.0 0.0 1.3],8)),
("Simple tetragonal",([1.1 0.0 0.0; 0.0 1.1 0.0; 0.0 0.0 1.3],16))
])

""" LG_G_test(LG,G,A)
    Test if the pointgroup in lattice coordinates, LG, and the pointgroup Cartesian coordinates, G, are consistent with the lattice A. That is, A is a similarity transform between the two representations.
"""
function LG_G_test(LG,G,A)
    return all([norm(A*lg*inv(A)-g) for (lg,g) ∈ zip(LG,G)] .< 5e-15)
end

@testset "LG/G tests" begin
    for (name, (A,nops)) ∈ BravaisLatticeList
        println(name, ",  nOps: ", nops)
        A = minkReduce(A)
        LG,G = pointGroup(A)
        @test isagroup(LG)
        @test isagroup(G)
        @test all([norm(A*lg*inv(A)-g) for (lg,g) ∈ zip(LG,G)] .< 5e-15)
    end
end

@testset "pointGroup_simple, random rotations" begin
    u = [1, 0, 0]
    v = [0.5, √3 / 2, 0]
    w = [0, 0, √(8 / 3)]
    @test length(pointGroup_simple(u, v, w)) == 24
    #@test all(abs.(det.(pointGroup_robust(u,v,w)[2])).==1.0)
    # Add a bit of noise
    @test length(pointGroup_simple(u, v.+[0,1e-9,0], w)) == 24
    # A bit more noise
    @test length(pointGroup_simple(u, v.+[0,1e-8,0], w)) ≠ 24
    #ideal hex lattice, 120° between basal plane vectors
    a = -u + v
    b = v
    c = w
    @test length(pointGroup_simple(a, b, c)) == 24
    # simple cubic lattice in unique orientation
    d, e, f = threeDrotation([1, 0, 0], [0, 1, 0], [0, 0, 1], π / 3, π / 5, π / 7)
    @test length(pointGroup_simple(d, e, f)) == 48
    # hexagonal lattice in unique orientation
    g, h, i = threeDrotation([1, 0, 0], [0.5, √3 / 2, 0], [0, 0, √(8 / 3)], π / 7, π / 11, π / 3)
    @test length(pointGroup_simple(g, h, i)) == 24
end

@testset "pointGroup_fast, small noise, large aspect ratios" begin
    # # Rhombohedral case
    u = [1, 1, 2]
    v = [1, 2, 1]
    w = [2, 1, 1]
    u, v, w = threeDrotation(u, v, w, π / 3, π / 5, π / 7)
    @test length(pointGroup_fast(u, v, w)) == 12

    # Cases of small tiny noise in input
    ϵ = 1.0e-9; u = [1, 0, ϵ]; v = [ϵ,1,0]; w = [0,0,1]
    @test length(pointGroup_fast(u,v,w))==48
    ϵ = 1.0e-7; u = [1, 0, ϵ]; v = [ϵ,1,0]; w = [0,0,1]
    @test length(pointGroup_fast(u,v,w))≠48
    # High aspect ratio cases
    a = 2^25; u = [1, 0, 0]; v = [0,a,0];  w = [0,0,1]
    @test length(pointGroup_fast(u,v,w))==16
    a = 2^26; u = [1, 0, 0]; v = [0,a,0];  w = [0,0,1]
    @test length(pointGroup_fast(u,v,w))≠16

    # Simple cubic example of snap function, ~1% noise
    a1 = [1+.01,0,0]; a2 = [0.,1-.01,0]; a3 = [0,0,1-.001];
    u,v,w = minkReduce(a1,a2,a3)
    ops,_ = pointGroup_robust(u,v,w;tol=5e-2)
    @test isagroup(ops)
    a,b,c,iops,rops = snapToSymmetry_SVD(u,v,w,ops)
    @test det([a b c])≈det([a1 a2 a3])
    @test norm(a)≈norm(b)≈norm(c)
    @test length(iops)==48
    @test isagroup(rops)
    @test LG_G_test(iops,rops,[a b c])

    # Simple tetragonal example of snap
    a1 = [1+.01,0,0]; a2 = [0.,1-.01,0]; a3 = [0,0,1.5];
    u,v,w = minkReduce(a1,a2,a3)
    ops,R = pointGroup_robust(u,v,w;tol=5e-2)
    @test isagroup(ops)
    @test isagroup(R)
    a,b,c,iops,rops = snapToSymmetry_SVD(u,v,w,ops)
    @test length(iops)==16
    @test norm(a)≈norm(b)
    @test det([a b c])≈det([a1 a2 a3])
    @test isagroup(iops)
    @test isagroup(rops)
    @test LG_G_test(iops,rops,[a b c])

    # Simple orthorhombic example of snap
    println("Orthorhombic example of snap")
    a1 = [1,0.01,-.005]; a2 = [0.001,2,-0.01]; a3 = [0.02,-0.003,1.5];
    u,v,w = minkReduce(a1,a2,a3)
    ops, R = pointGroup_robust(u,v,w;tol=5e-2)
    a,b,c,iops,rops = snapToSymmetry_SVD(u,v,w,ops)
    @test length(iops)==8
    @test det([a b c])≈det([a1 a2 a3])
    @test isagroup(iops)
    @test isagroup(rops)
    @test LG_G_test(iops,rops,[a b c])


    println("Example of large aspect ratio: 256 (should still succeed)")
    a1 = [1/16 - .0001, .001, .0001]
    a2 = [-0.001, 16-.0001, -.0001]
    a3 = [-.0001, .0001, 1.51]
    ops,_ = pointGroup(minkReduce(hcat(a1,a2,a3)))
    @test length(ops)==8

    println("Example of large aspect ratio: 500 (should still succeed)")
    a1 = [1/10 - .0001, .001, .0001]
    a2 = [-0.001, 50-.0001, -.0001]
    a3 = [-.0001, .0001, 1.51]
    ops,_ = pointGroup(minkReduce(hcat(a1,a2,a3)))
    @test length(ops)==8

    println("Example of large aspect ratio: 512 (fails, even with large tolerance)")
    a1 = [1/16 - .0001, .001, .0001]
    a2 = [-0.001, 32-.0001, -.0001]
    a3 = [-.0001, .0001, 1.51]
    ops,_ = pointGroup_robust(minkReduce(a1,a2,a3)[1:3]...;tol=9e-1)
    @test length(ops)==4

    println("Example of large aspect ratio: 1024 (fails, even with large tolerance)")
    a1 = [1/32 - .0001, .001, .0001]
    a2 = [-0.001, 32-.0001, -.0001]
    a3 = [-.0001, .0001, 1.51]
    ops,_ = pointGroup_robust(minkReduce(a1,a2,a3)[1:3]...;tol=9e-1)
    @test length(ops)==4


    println("Example of large aspect ratio: 1024 (succeeds, smaller noise,)")
    a1 = [1/32 - .0001, .00001, .00001]
    a2 = [-0.0001, 32-.0001, -.00001]
    a3 = [-.0001, .00001, 1.51]
    ops,_ = pointGroup_robust(minkReduce(a1,a2,a3)[1:3]...;tol=1e-1)
    @test length(ops)==8
end

@testset "Random noise, all lattices" begin
    for (name, (A,nops)) ∈ BravaisLatticeList
        println(name, ",  nOps: ", nops)
        a = 1e-0; Navg =100; Nsteps = 40; maxtol = 1e-1
        tols = logrange(5e-12,maxtol,20)
        plim = logrange(5e-6,1e0,Nsteps) # Size of noise to test. noise bigger that 1/12 tol will get skipped
        data = Vector{Float64}(undef,Nsteps)
        for (i,tol) ∈ enumerate(tols)
            for ε ∈ plim
                if tol < 12*ε/a; break; end # Anything less than 5*ε/a will fail for many cases. Slight tetragonal distortions from bcc/fcc will be found as cubic if tol is too big or distortion is too small.
                Atest =  hcat(minkReduce(eachcol(A*a + (2*rand(3,3).-1)*ε*a)...)[1:3]...)
                nSuccess = count([length(pointGroup(Atest;tol=tol)[1])==nops for _ ∈ 1:Navg])
                if nSuccess != Navg
                    @show Atest
                    error("Pointgroup size is not $nops.  tol: ", tol, "  ε: ", ε,"  nSuccess: ", nSuccess)
                end
                @test nSuccess == Navg
            end
        end
    end
end

@testset "Aspect ratio tests" begin
    println("Rhombohedral case")

end

@testset "Near-boundary tetragonal" begin
    # A = diag(1, 1, 1+ε) is genuinely tetragonal for ε > 0 (16 ops). At tol ≫ ε,
    # pointGroup_robust over-promotes to cubic (48) because the cubic-only ops
    # (e.g. 3-folds about body diagonals) have integer-matrix residuals of O(ε)
    # and pass the tol filter. See research.md §2.1 and test/nearMissBoundary.jl.
    for ε ∈ [1e-3, 1e-5, 1e-8]
        A = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0+ε]
        u, v, w = minkReduce(eachcol(A)...)[1:3]
        tight_tol = ε / 100
        loose_tol = 100 * ε
        @test length(pointGroup_simple(u, v, w)) == 16
        @test length(pointGroup_robust(u, v, w; tol=tight_tol)[1]) == 16
        @test length(pointGroup_robust(u, v, w; tol=loose_tol)[1]) == 48
        @test_logs (:warn, r"near a symmetry boundary") match_mode=:any pointGroup_robust(u, v, w; tol=loose_tol, verify_stable=true)
        @test length(pointGroup_robust(u, v, w; tol=tight_tol, verify_stable=true)[1]) == 16
    end
end

@testset "Crystal construction" begin
    A = Matrix{Float64}(I, 3, 3)
    r = reshape([0.0, 0.0, 0.0], 3, 1)

    # Basic fractional construction
    c = Crystal(A, r, [1]; coords=:fractional)
    @test c.A == A
    @test c.r == r
    @test c.types == [1]
    @test c isa Crystal{Int}

    # Symbol types
    c_sym = Crystal(A, r, [:Fe]; coords=:fractional)
    @test c_sym isa Crystal{Symbol}
    @test c_sym.types == [:Fe]

    # String types
    c_str = Crystal(A, r, ["Fe"]; coords=:fractional)
    @test c_str isa Crystal{String}

    # Cartesian input converts to fractional
    A2 = 2.0 * Matrix{Float64}(I, 3, 3)
    r_cart = reshape([1.0, 1.0, 1.0], 3, 1)   # Cartesian (1,1,1) in 2× cubic
    c_cart = Crystal(A2, r_cart, [1]; coords=:cartesian)
    @test c_cart.r ≈ reshape([0.5, 0.5, 0.5], 3, 1)

    # Integer input converts to Float64
    A_int = [1 0 0; 0 1 0; 0 0 1]
    r_int = reshape([0, 0, 0], 3, 1)
    c_int = Crystal(A_int, r_int, [1]; coords=:fractional)
    @test c_int.A isa Matrix{Float64}
    @test c_int.r isa Matrix{Float64}

    # Three-vector constructor
    c_3v = Crystal([1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0], r, [1]; coords=:fractional)
    @test c_3v.A == A

    # Validation errors
    @test_throws Exception Crystal(A, r, [1])                             # missing coords
    @test_throws ErrorException Crystal(A, r, [1]; coords=:oops)         # bad coords
    @test_throws ErrorException Crystal([1.0 0.0; 0.0 1.0], r, [1]; coords=:fractional)  # A not 3×3
    @test_throws ErrorException Crystal(A, [0.0 0.0; 0.0 0.0], [1, 2]; coords=:fractional)  # r rows
    @test_throws ErrorException Crystal(A, r, [1, 2]; coords=:fractional) # types length mismatch

    # Empty crystal rejected (no atoms)
    empty_r = Matrix{Float64}(undef, 3, 0)
    @test_throws ErrorException Crystal(A, empty_r, Int[]; coords=:fractional)

    # Singular A rejected (coplanar basis vectors: rank 2, det = 0)
    A_sing = [1.0 0.0 1.0; 0.0 1.0 1.0; 0.0 0.0 0.0]
    @test_throws ErrorException Crystal(A_sing, r, [1]; coords=:fractional)
    # Same test at a different scale — the check is scale-invariant
    @test_throws ErrorException Crystal(1e6 .* A_sing, r, [1]; coords=:fractional)
    @test_throws ErrorException Crystal(1e-6 .* A_sing, r, [1]; coords=:fractional)
end

@testset "fractional / cartesian / default_pos_tol" begin
    A = [2.0 0.0 0.0; 0.0 3.0 0.0; 0.0 0.0 4.0]
    r_frac = reshape([0.5, 0.5, 0.5], 3, 1)
    c = Crystal(A, r_frac, [1]; coords=:fractional)

    @test fractional(c) == r_frac
    @test cartesian(c) ≈ reshape([1.0, 1.5, 2.0], 3, 1)

    # default_pos_tol: V = 24, N = 1 → 0.01 · 24^(1/3)
    @test default_pos_tol(c) ≈ 0.01 * 24^(1/3)

    # Two atoms, unit cell: V = 1, N = 2 → 0.01 · (1/2)^(1/3)
    A_unit = Matrix{Float64}(I, 3, 3)
    r2 = [0.0 0.5; 0.0 0.5; 0.0 0.5]
    c2 = Crystal(A_unit, r2, [1, 2]; coords=:fractional)
    @test default_pos_tol(c2) ≈ 0.01 * (1/2)^(1/3)
end

@testset "isSpacegroupOp: Phase 1 trivial cases" begin
    A = Matrix{Float64}(I, 3, 3)
    r = reshape([0.0, 0.0, 0.0], 3, 1)
    c = Crystal(A, r, [1]; coords=:fractional)
    tol = 1e-8
    I3 = Matrix{Int}(I, 3, 3)

    # Identity op
    @test isSpacegroupOp(I3, [0.0, 0.0, 0.0], c; tol=tol)

    # Pure lattice translation (should fold to identity mod 1)
    @test isSpacegroupOp(I3, [1.0, 0.0, 0.0], c; tol=tol)
    @test isSpacegroupOp(I3, [-1.0, 2.0, -3.0], c; tol=tol)

    # 4-fold rotation about z — preserves (0,0,0)
    R4 = [0 -1 0; 1 0 0; 0 0 1]
    @test isSpacegroupOp(R4, [0.0, 0.0, 0.0], c; tol=tol)

    # Non-trivial translation that leaves the atom's image unmatched
    @test !isSpacegroupOp(I3, [0.5, 0.0, 0.0], c; tol=tol)

    # CsCl: Cs at (0,0,0), Cl at (½,½,½), full Pm3̄m
    c_cscl = Crystal(A, [0.0 0.5; 0.0 0.5; 0.0 0.5], [:Cs, :Cl]; coords=:fractional)
    @test isSpacegroupOp(I3, [0.0, 0.0, 0.0], c_cscl; tol=tol)
    @test isSpacegroupOp(R4, [0.0, 0.0, 0.0], c_cscl; tol=tol)
    @test isSpacegroupOp(-I3, [0.0, 0.0, 0.0], c_cscl; tol=tol)  # inversion

    # Type-preservation: (I, (½,½,½)) maps Cs→Cl and Cl→Cs.
    # Should be false because types don't match after the map.
    @test !isSpacegroupOp(I3, [0.5, 0.5, 0.5], c_cscl; tol=tol)

    # Asymmetric two-atom crystal: only identity is a symmetry
    c_asym = Crystal(A, [0.0 0.3; 0.0 0.7; 0.0 0.2], [:A, :B]; coords=:fractional)
    @test isSpacegroupOp(I3, [0.0, 0.0, 0.0], c_asym; tol=tol)
    @test !isSpacegroupOp(R4, [0.0, 0.0, 0.0], c_asym; tol=tol)
    @test !isSpacegroupOp(-I3, [0.0, 0.0, 0.0], c_asym; tol=tol)

    # Default tol kwarg path (no explicit tol)
    @test isSpacegroupOp(I3, [0.0, 0.0, 0.0], c_cscl)
    @test !isSpacegroupOp(I3, [0.5, 0.5, 0.5], c_cscl)

    # R / τ shape validation
    @test_throws ErrorException isSpacegroupOp([1 0; 0 1], [0.0, 0.0, 0.0], c; tol=tol)
    @test_throws ErrorException isSpacegroupOp(I3, [0.0, 0.0], c; tol=tol)

    # Pure-translation symmetry (R = I but τ ≠ 0, and the result is TRUE).
    # Two same-type atoms at (¼,0,0) and (¾,0,0) are exchanged by τ = (½,0,0).
    # Exercises the (R=I, τ≠0, should-succeed) path and the injectivity
    # guard in a case where claimed-skip actually matters: image of atom 1
    # must claim atom 2, then image of atom 2 must claim atom 1 — skipping
    # the already-claimed atom 2 en route.
    c_halftrans = Crystal(A, [0.25 0.75; 0.0 0.0; 0.0 0.0], [:A, :A]; coords=:fractional)
    @test isSpacegroupOp(I3, [0.5, 0.0, 0.0], c_halftrans; tol=tol)
    @test isSpacegroupOp(I3, [-0.5, 0.0, 0.0], c_halftrans; tol=tol)
    @test !isSpacegroupOp(I3, [0.0, 0.5, 0.0], c_halftrans; tol=tol)

    # Mod-1 boundary: atom at 0, translation of (1 − 1e-12) — image sits at
    # ≈ 0.999999999999 after the wrap. The signed-diff formula must bring
    # |Δ| back near 0. Catches sign typos in `mod.(Δ .+ 0.5, 1.0) .- 0.5`.
    @test isSpacegroupOp(I3, [1.0 - 1e-12, 0.0, 0.0], c; tol=1e-8)
    @test isSpacegroupOp(I3, [-1.0 + 1e-12, 0.0, 0.0], c; tol=1e-8)
end

@testset "SpacegroupOp methods" begin
    I3 = Matrix{Int}(I, 3, 3)

    # Identity
    e = one(SpacegroupOp)
    @test e.R == I3
    @test e.τ == zeros(3)

    # Canonicalization of τ at construction
    op_a = SpacegroupOp(I3, [0.0, 0.0, 0.0])
    op_b = SpacegroupOp(I3, [1.0, 0.0, 0.0])         # folds to [0,0,0]
    op_c = SpacegroupOp(I3, [2.5, -0.5, 1.25])       # folds to [0.5, 0.5, 0.25]
    @test op_a.τ == zeros(3)
    @test op_b.τ == zeros(3)
    @test op_c.τ ≈ [0.5, 0.5, 0.25]
    @test op_a == op_b                                # default == works
    @test hash(op_a) == hash(op_b)                    # hash consistent with ==

    # Composition: op1 * op2 means "apply op2 first, then op1"
    R1 = [0 -1 0; 1 0 0; 0 0 1]    # 4-fold about z
    op1 = SpacegroupOp(R1, [0.5, 0.0, 0.0])
    op2 = SpacegroupOp(R1, [0.0, 0.5, 0.0])
    op12 = op1 * op2
    @test op12.R == R1 * R1
    @test op12.τ ≈ mod.(R1 * [0.0, 0.5, 0.0] + [0.5, 0.0, 0.0], 1.0)

    # Composition is associative
    op3 = SpacegroupOp(R1, [0.25, 0.0, 0.0])
    @test ((op1 * op2) * op3) == (op1 * (op2 * op3))

    # Identity composition
    @test (e * op1) == op1
    @test (op1 * e) == op1

    # Inverse: op * inv(op) == e, mod 1
    op = SpacegroupOp(R1, [0.25, 0.0, 0.0])
    @test (op * inv(op)) == e
    @test (inv(op) * op) == e
    @test inv(inv(op)) == op

    # Callable: op(r) applies the op to position r
    r = [0.1, 0.2, 0.3]
    @test op(r) ≈ mod.(R1 * r + op.τ, 1.0)

    # Equality mod-1
    @test SpacegroupOp(I3, [0.0, 0.0, 0.0]) == SpacegroupOp(I3, [1.0, 0.0, 0.0])
    @test SpacegroupOp(I3, [0.0, 0.0, 0.0]) != SpacegroupOp(I3, [0.5, 0.0, 0.0])
    @test SpacegroupOp(R1, [0.0, 0.0, 0.0]) != SpacegroupOp(I3, [0.0, 0.0, 0.0])

    # toCartesian on identity through a non-diagonal A
    A = [1.0 0.3 0.0; 0.0 1.1 0.0; 0.0 0.0 1.5]
    Rc, τc = toCartesian(one(SpacegroupOp), A)
    @test Rc ≈ Matrix{Float64}(I, 3, 3)
    @test τc ≈ zeros(3)

    # toCartesian on a non-trivial op: the Cartesian R must equal the
    # rotation "A applied to lattice coords then A⁻¹ back".
    Rc2, τc2 = toCartesian(SpacegroupOp(R1, [0.5, 0.0, 0.0]), A)
    @test Rc2 ≈ A * R1 * inv(A)
    @test τc2 ≈ A * [0.5, 0.0, 0.0]

    # inv validation: a non-unimodular R should error
    @test_throws ErrorException inv(SpacegroupOp([2 0 0; 0 1 0; 0 0 1], zeros(3)))
end

@testset "spacegroup: Phase 2 core cases" begin
    I3 = Matrix{Int}(I, 3, 3)

    # 4.2.1 Exit criterion: simple cubic, 1 atom at origin → 48 ops, τ=0 all
    A_cubic = Matrix{Float64}(I, 3, 3)
    c_sc = Crystal(A_cubic, reshape([0.0, 0.0, 0.0], 3, 1), [:X]; coords=:fractional)
    ops = spacegroup(c_sc)
    @test length(ops) == 48
    @test all(op.τ ≈ zeros(3) for op in ops)
    @test all(abs(det(op.R)) == 1 for op in ops)

    # 4.2.2 Simple cubic shifted to body-centre: still 48 ops
    c_sc_shift = Crystal(A_cubic, reshape([0.5, 0.5, 0.5], 3, 1), [:X]; coords=:fractional)
    @test length(spacegroup(c_sc_shift)) == 48

    # 4.2.3 Simple tetragonal (c ≠ a), 1 atom at origin → 16 ops
    A_tet = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.5]
    c_tet = Crystal(A_tet, reshape([0.0, 0.0, 0.0], 3, 1), [:X]; coords=:fractional)
    @test length(spacegroup(c_tet)) == 16

    # 4.2.4 Simple orthorhombic, 1 atom → 8 ops
    A_ortho = [1.0 0.0 0.0; 0.0 1.2 0.0; 0.0 0.0 1.5]
    c_ortho = Crystal(A_ortho, reshape([0.0, 0.0, 0.0], 3, 1), [:X]; coords=:fractional)
    @test length(spacegroup(c_ortho)) == 8

    # 4.2.5 Triclinic, 1 atom → 2 ops (E and inversion-with-compensating-τ)
    A_tri = [1.0 0.3 0.2; 0.0 1.1 0.4; 0.0 0.0 0.9]
    c_tri = Crystal(A_tri, reshape([0.0, 0.0, 0.0], 3, 1), [:X]; coords=:fractional)
    @test length(spacegroup(c_tri)) == 2

    # 4.2.6 CsCl: Cs at (0,0,0), Cl at (½,½,½), Pm3̄m → 48 ops, τ=0 all
    c_cscl = Crystal(A_cubic, [0.0 0.5; 0.0 0.5; 0.0 0.5], [:Cs, :Cl];
                     coords=:fractional)
    ops_cscl = spacegroup(c_cscl)
    @test length(ops_cscl) == 48
    @test all(op.τ ≈ zeros(3) for op in ops_cscl)

    # 4.2.7 Every returned op passes isSpacegroupOp (cross-check)
    for op in ops
        @test isSpacegroupOp(op.R, op.τ, c_sc)
    end
    for op in ops_cscl
        @test isSpacegroupOp(op.R, op.τ, c_cscl)
    end

    # 4.2.8 Group closure: op1 * op2 ∈ ops (spot check)
    for _ in 1:10
        a, b = rand(ops), rand(ops)
        @test (a * b) ∈ ops
    end
    for _ in 1:10
        a, b = rand(ops_cscl), rand(ops_cscl)
        @test (a * b) ∈ ops_cscl
    end

    # 4.2.9 Identity is at index 1
    @test ops[1] == one(SpacegroupOp)
    @test ops_cscl[1] == one(SpacegroupOp)
    @test spacegroup(c_tri)[1] == one(SpacegroupOp)
    @test spacegroup(c_tet)[1] == one(SpacegroupOp)
    @test spacegroup(c_ortho)[1] == one(SpacegroupOp)

    # 4.2.10 Inverse of every op is also in the group
    for op in ops
        @test inv(op) ∈ ops
    end
end

@testset "spacegroup: Phase 3 known crystals" begin
    I3 = Matrix{Int}(I, 3, 3)

    # 1) NaCl, conventional cubic cell (Fm3̄m = #225, order 192).
    # Lattice is simple-cubic (A = I), but the atomic basis has FCC
    # centering: 4 Na at (0,0,0), (½,½,0), (½,0,½), (0,½,½) and 4 Cl at
    # those + (½,½,½). The space group picks up 4× centering translations
    # on top of the 48 cubic rotations: 48 · 4 = 192.
    A_cubic = Matrix{Float64}(I, 3, 3)
    Na_pos = [0.0 0.5 0.5 0.0;
              0.0 0.5 0.0 0.5;
              0.0 0.0 0.5 0.5]
    Cl_pos = Na_pos .+ [0.5, 0.5, 0.5]
    c_nacl = Crystal(A_cubic, hcat(Na_pos, Cl_pos),
                     [fill(:Na, 4); fill(:Cl, 4)]; coords=:fractional)
    ops_nacl = spacegroup(c_nacl)
    @test length(ops_nacl) == 192
    @test ops_nacl[1] == one(SpacegroupOp)
    # Closure (spot check)
    for _ in 1:10
        a, b = rand(ops_nacl), rand(ops_nacl)
        @test (a * b) ∈ ops_nacl
    end
    # Every op must be a valid symmetry
    for op in ops_nacl
        @test isSpacegroupOp(op.R, op.τ, c_nacl)
    end
    # NaCl's 4 pure-translation ops (R = I, τ = FCC centering) must be
    # present with R = I.
    identity_R_ops = filter(op -> op.R == I3, ops_nacl)
    @test length(identity_R_ops) == 4
    # Among them, τ ∈ {(0,0,0), (½,½,0), (½,0,½), (0,½,½)} — verify
    # each of the three non-zero centering translations appears.
    τs_sorted = sort([op.τ for op in identity_R_ops]; by=v -> (v[1], v[2], v[3]))
    @test τs_sorted[1] ≈ [0.0, 0.0, 0.0]
    @test any(τ -> τ ≈ [0.5, 0.5, 0.0], τs_sorted)
    @test any(τ -> τ ≈ [0.5, 0.0, 0.5], τs_sorted)
    @test any(τ -> τ ≈ [0.0, 0.5, 0.5], τs_sorted)

    # 2) Diamond, conventional cubic cell (Fd3̄m = #227, order 192).
    # Same cubic lattice and FCC-centred arrangement, but the 2-atom
    # primitive basis is (0,0,0) and (¼,¼,¼). Non-symmorphic: glide
    # planes appear as ops with non-zero τ at non-lattice-centring values.
    fcc_pos = [0.0 0.5 0.5 0.0;
               0.0 0.5 0.0 0.5;
               0.0 0.0 0.5 0.5]
    fcc_shifted = fcc_pos .+ [0.25, 0.25, 0.25]
    c_diamond = Crystal(A_cubic, hcat(fcc_pos, fcc_shifted), fill(:C, 8);
                        coords=:fractional)
    ops_diamond = spacegroup(c_diamond)
    @test length(ops_diamond) == 192
    @test ops_diamond[1] == one(SpacegroupOp)
    for op in ops_diamond
        @test isSpacegroupOp(op.R, op.τ, c_diamond)
    end
    # Non-symmorphic signature: some ops have τ ≉ any FCC centering.
    # FCC centerings are (0,0,0), (½,½,0), (½,0,½), (0,½,½).
    function is_fcc_centering(τ, atol=1e-8)
        centerings = ([0.0,0.0,0.0], [0.5,0.5,0.0], [0.5,0.0,0.5], [0.0,0.5,0.5])
        any(c -> all(abs.(mod.(τ .- c .+ 0.5, 1.0) .- 0.5) .< atol), centerings)
    end
    # At least one op must have a glide-style τ (not a pure FCC centering)
    @test any(!is_fcc_centering(op.τ) for op in ops_diamond)

    # 3) Hexagonal close-packed (P6₃/mmc = #194, order 24). Non-symmorphic
    # — the 6₃ screw axis.
    # Ideal HCP: c/a = √(8/3). Lattice a1=(a,0,0), a2=(-a/2, a√3/2, 0),
    # a3=(0,0,c). Two atoms at Wyckoff 2c: (⅓, ⅔, ¼) and (⅔, ⅓, ¾).
    a_hcp = 1.0
    c_hcp = sqrt(8/3)
    A_hcp = [a_hcp  -a_hcp/2       0.0;
             0.0    a_hcp*√3/2     0.0;
             0.0    0.0            c_hcp]
    r_hcp = [1/3 2/3;
             2/3 1/3;
             1/4 3/4]
    c_hcp_crystal = Crystal(A_hcp, r_hcp, [:X, :X]; coords=:fractional)
    ops_hcp = spacegroup(c_hcp_crystal)
    @test length(ops_hcp) == 24
    @test ops_hcp[1] == one(SpacegroupOp)
    for op in ops_hcp
        @test isSpacegroupOp(op.R, op.τ, c_hcp_crystal)
    end
    # Non-symmorphic signature: at least one op has τ with a ½ in the c
    # direction (the 6₃ screw translates by c/2 along c).
    @test any(abs(op.τ[3] - 0.5) < 1e-8 for op in ops_hcp)
    # Group closure (spot check)
    for _ in 1:10
        a, b = rand(ops_hcp), rand(ops_hcp)
        @test (a * b) ∈ ops_hcp
    end
end

@testset "spacegroup: AFLOW Part 1 (§3.5 validation, seed cases)" begin
    # Structures transcribed from AFLOW Part 1 POSCAR appendix
    # (Mehl, Hicks, Toher, Levy, Hanson, Hart, Curtarolo, Comput. Mater.
    # Sci. 136 (2017) S1–S828). Starting minimal — one example to confirm
    # the transcribe-and-test pipeline; more to follow per §3.13 of
    # plan.md.

    # ── FeS₂ pyrite ─────────────────────────────────────────────────
    # Example POSCAR from paper p. S22.
    # Prototype AB2_cP12_205_a_c, space group Pa-3 (#205), order 24.
    # Simple cubic lattice, a = 5.417 Å. 12 atoms in the primitive cell
    # (4 Fe at Wyckoff 4a, 8 S at Wyckoff 8c with internal parameter
    # x₂ = 0.38510 — so S coords are (0.11490, 0.38510, 0.61490, 0.88510)
    # permutations generated by the Pa-3 space-group ops).
    # Non-symmorphic: Pa-3 contains glide / screw operations, so some
    # of the 24 ops have non-zero τ.
    let
        a = 5.417
        A = a * Matrix{Float64}(I, 3, 3)
        Fe_pos = [0.0 0.0 0.5 0.5;
                  0.0 0.5 0.0 0.5;
                  0.0 0.5 0.5 0.0]
        S_pos = [0.11490 0.11490 0.38510 0.38510 0.61490 0.61490 0.88510 0.88510;
                 0.61490 0.88510 0.11490 0.38510 0.61490 0.88510 0.11490 0.38510;
                 0.88510 0.38510 0.88510 0.38510 0.61490 0.11490 0.61490 0.11490]
        c_pyrite = Crystal(A, hcat(Fe_pos, S_pos),
                           [fill(:Fe, 4); fill(:S, 8)]; coords=:fractional)
        ops = spacegroup(c_pyrite)
        @test length(ops) == 24
        @test ops[1] == one(SpacegroupOp)
        for op in ops
            @test isSpacegroupOp(op.R, op.τ, c_pyrite)
        end
        # Non-symmorphic signature: at least one op with a non-zero τ
        # component at 1/2 (Pa-3 has (½,0,½)-style glide translations).
        @test any(op -> any(t -> abs(t - 0.5) < 1e-8, op.τ), ops)
        # Closure (spot check)
        for _ in 1:10
            a, b = rand(ops), rand(ops)
            @test (a * b) ∈ ops
        end
    end
end

