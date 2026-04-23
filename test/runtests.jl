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

@testset "spacegroup stub raises" begin
    A = Matrix{Float64}(I, 3, 3)
    c = Crystal(A, reshape([0.0, 0.0, 0.0], 3, 1), [1]; coords=:fractional)
    @test_throws ErrorException spacegroup(c)
end

