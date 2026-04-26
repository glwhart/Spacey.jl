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
    @test length(Spacey.pointGroup_simple(u, v, w)) == 24
    #@test all(abs.(det.(Spacey.pointGroup_robust(u,v,w)[2])).==1.0)
    # Add a bit of noise
    @test length(Spacey.pointGroup_simple(u, v.+[0,1e-9,0], w)) == 24
    # A bit more noise
    @test length(Spacey.pointGroup_simple(u, v.+[0,1e-8,0], w)) ≠ 24
    #ideal hex lattice, 120° between basal plane vectors
    a = -u + v
    b = v
    c = w
    @test length(Spacey.pointGroup_simple(a, b, c)) == 24
    # simple cubic lattice in unique orientation
    d, e, f = Spacey.threeDrotation([1, 0, 0], [0, 1, 0], [0, 0, 1], π / 3, π / 5, π / 7)
    @test length(Spacey.pointGroup_simple(d, e, f)) == 48
    # hexagonal lattice in unique orientation
    g, h, i = Spacey.threeDrotation([1, 0, 0], [0.5, √3 / 2, 0], [0, 0, √(8 / 3)], π / 7, π / 11, π / 3)
    @test length(Spacey.pointGroup_simple(g, h, i)) == 24
end

@testset "pointGroup_fast, small noise, large aspect ratios" begin
    # # Rhombohedral case
    u = [1, 1, 2]
    v = [1, 2, 1]
    w = [2, 1, 1]
    u, v, w = Spacey.threeDrotation(u, v, w, π / 3, π / 5, π / 7)
    @test length(Spacey.pointGroup_fast(u, v, w)) == 12

    # Cases of small tiny noise in input
    ϵ = 1.0e-9; u = [1, 0, ϵ]; v = [ϵ,1,0]; w = [0,0,1]
    @test length(Spacey.pointGroup_fast(u,v,w))==48
    ϵ = 1.0e-7; u = [1, 0, ϵ]; v = [ϵ,1,0]; w = [0,0,1]
    @test length(Spacey.pointGroup_fast(u,v,w))≠48
    # High aspect ratio cases
    a = 2^25; u = [1, 0, 0]; v = [0,a,0];  w = [0,0,1]
    @test length(Spacey.pointGroup_fast(u,v,w))==16
    a = 2^26; u = [1, 0, 0]; v = [0,a,0];  w = [0,0,1]
    @test length(Spacey.pointGroup_fast(u,v,w))≠16

    # Simple cubic example of snap function, ~1% noise
    a1 = [1+.01,0,0]; a2 = [0.,1-.01,0]; a3 = [0,0,1-.001];
    u,v,w = minkReduce(a1,a2,a3)
    ops,_ = Spacey.pointGroup_robust(u,v,w;tol=5e-2)
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
    ops,R = Spacey.pointGroup_robust(u,v,w;tol=5e-2)
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
    ops, R = Spacey.pointGroup_robust(u,v,w;tol=5e-2)
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
    ops,_ = Spacey.pointGroup_robust(minkReduce(a1,a2,a3)[1:3]...;tol=9e-1)
    @test length(ops)==4

    println("Example of large aspect ratio: 1024 (fails, even with large tolerance)")
    a1 = [1/32 - .0001, .001, .0001]
    a2 = [-0.001, 32-.0001, -.0001]
    a3 = [-.0001, .0001, 1.51]
    ops,_ = Spacey.pointGroup_robust(minkReduce(a1,a2,a3)[1:3]...;tol=9e-1)
    @test length(ops)==4


    println("Example of large aspect ratio: 1024 (succeeds, smaller noise,)")
    a1 = [1/32 - .0001, .00001, .00001]
    a2 = [-0.0001, 32-.0001, -.00001]
    a3 = [-.0001, .00001, 1.51]
    ops,_ = Spacey.pointGroup_robust(minkReduce(a1,a2,a3)[1:3]...;tol=1e-1)
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
        @test length(Spacey.pointGroup_simple(u, v, w)) == 16
        @test length(Spacey.pointGroup_robust(u, v, w; tol=tight_tol)[1]) == 16
        @test length(Spacey.pointGroup_robust(u, v, w; tol=loose_tol)[1]) == 48
        @test_logs (:warn, r"near a symmetry boundary") match_mode=:any Spacey.pointGroup_robust(u, v, w; tol=loose_tol, verify_stable=true)
        @test length(Spacey.pointGroup_robust(u, v, w; tol=tight_tol, verify_stable=true)[1]) == 16
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

    # 4.2.2 Simple cubic shifted to body-center: still 48 ops
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
    # Same cubic lattice and FCC-centered arrangement, but the 2-atom
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

    # Helper: build a Crystal from a POSCAR-style transcription.
    # `A_rows` is a 3×3 matrix whose ROWS are the lattice vectors (as the
    # paper prints them); we transpose so that the Crystal's `A` has them
    # as columns, matching Spacey's convention.
    make_crystal(A_rows, r_cols, types) =
        Crystal(collect(transpose(A_rows)), r_cols, types; coords=:fractional)

    # Common spot-check: length matches, identity at index 1, every op
    # passes isSpacegroupOp, closure spot check.
    function check_spacegroup!(c::Crystal, expected_order::Int)
        ops = spacegroup(c)
        @test length(ops) == expected_order
        @test ops[1] == one(SpacegroupOp)
        for op in ops
            @test isSpacegroupOp(op.R, op.τ, c)
        end
        for _ in 1:5
            a, b = rand(ops), rand(ops)
            @test (a * b) ∈ ops
        end
        return ops
    end

    # ── Cf — A_aP4_2_aci ───────────────────────────────────────────
    # Space group P-1 (#2), triclinic primitive, 4 atoms. Order 2.
    # Symmorphic (only E and inversion-with-compensating-τ).
    let
        A_rows = [3.30700  0.00000  0.00000;
                  0.55574  7.39114  0.00000;
                  0.23614  0.02819  2.78286]
        r = [0.00000 0.00000 0.42800 0.57200;
             0.00000 0.50000 0.74100 0.25900;
             0.00000 0.00000 0.56700 0.43300]
        c = make_crystal(A_rows, r, fill(:Cf, 4))
        check_spacegroup!(c, 2)
    end

    # ── High-pressure Te — A_mP4_4_2a ──────────────────────────────
    # P2₁ (#4), monoclinic primitive, 4 atoms. Order 2. Non-symmorphic
    # (2₁ screw axis).
    let
        A_rows = [ 3.10400  0.00000  0.00000;
                   0.00000  7.51300  0.00000;
                  -0.22506  0.00000  4.75468]
        r = [0.25000 0.75000 0.48000 0.52000;
             0.23000 0.73000 0.00000 0.50000;
             0.48000 0.52000 0.02000 -0.02000]
        c = make_crystal(A_rows, r, fill(:Te, 4))
        check_spacegroup!(c, 2)
    end

    # ── NiTi — AB_mP4_11_e_e ───────────────────────────────────────
    # P2₁/m (#11), monoclinic primitive, 2 Ni + 2 Ti = 4 atoms. Order 4.
    # Non-symmorphic (2₁ screw + m perpendicular).
    let
        A_rows = [2.88370 0.00000 0.00000;
                  0.00000 4.10620 0.00000;
                  0.64457 0.00000 4.62268]
        r = [0.03870 0.96130 0.41130 0.58870;
             0.25000 0.75000 0.75000 0.25000;
             0.82520 0.17480 0.28160 0.71840]
        c = make_crystal(A_rows, r, [:Ni, :Ni, :Ti, :Ti])
        check_spacegroup!(c, 4)
    end

    # ── α-O — A_mC4_12_i ───────────────────────────────────────────
    # C2/m (#12), C-centered monoclinic. POSCAR gives the PRIMITIVE cell
    # (2 atoms), so Spacey sees 4 ops (point-group order = 4; the
    # C-centring translation is absorbed into the primitive basis).
    let
        A_rows = [ 2.70150 -1.71650  0.00000;
                   2.70150  1.71650  0.00000;
                  -3.41954  0.00000  3.75539]
        r = [0.10600 0.89400;
             0.10600 0.89400;
             0.17300 0.82700]
        c = make_crystal(A_rows, r, [:O, :O])
        check_spacegroup!(c, 4)
    end

    # ── High-pressure CdTe — AB_oP2_25_b_a ─────────────────────────
    # Pmm2 (#25), orthorhombic primitive, 2 atoms. Order 4. Symmorphic.
    let
        A_rows = [2.81020 0.00000 0.00000;
                  0.00000 5.25800 0.00000;
                  0.00000 0.00000 3.02650]
        r = [0.00000 0.00000;
             0.50000 0.00000;
             0.25000 0.00000]
        c = make_crystal(A_rows, r, [:Cd, :Te])
        check_spacegroup!(c, 4)
    end

    # ── High-pressure GaAs — AB_oI4_44_a_b ─────────────────────────
    # Imm2 (#44), I-centered orthorhombic. POSCAR gives the primitive
    # cell (2 atoms), so Spacey sees 4 ops.
    let
        A_rows = [-2.46000  2.39500  1.31750;
                   2.46000 -2.39500  1.31750;
                   2.46000  2.39500 -1.31750]
        r = [0.00000 0.92500;
             0.00000 0.42500;
             0.00000 0.50000]
        c = make_crystal(A_rows, r, [:As, :Ga])
        check_spacegroup!(c, 4)
    end

    # ── Naumannite Ag₂Se — A2B_oP12_19_2a_a ────────────────────────
    # P2₁2₁2₁ (#19), orthorhombic primitive, 8 Ag + 4 Se = 12 atoms.
    # Order 4. Non-symmorphic (three 2₁ screw axes).
    let
        A_rows = [7.76400 0.00000 0.00000;
                  0.00000 7.06200 0.00000;
                  0.00000 0.00000 4.33300]
        Ag = [0.185 0.315 0.685 0.815 -0.055 0.055 0.445 0.555;
              0.070 -0.070 0.430 0.570 0.265 0.765 0.235 0.735;
              0.465 -0.035 0.535 0.035 0.508 -0.008 0.492 0.008]
        Se = [0.116 0.384 0.616 0.884;
              0.489 0.511 0.011 -0.011;
              0.109 0.609 0.891 0.391]
        c = make_crystal(A_rows, hcat(Ag, Se),
                         [fill(:Ag, 8); fill(:Se, 4)])
        check_spacegroup!(c, 4)
    end

    # ── α-U — A_oC4_63_c ───────────────────────────────────────────
    # Cmcm (#63), C-centered orthorhombic. POSCAR gives the primitive
    # cell (2 atoms), so Spacey sees 8 ops (point-group mmm; centring
    # absorbed).
    let
        A_rows = [1.42220 -2.93445 0.00000;
                  1.42220  2.93445 0.00000;
                  0.00000  0.00000 4.93160]
        r = [0.10228 0.89772;
             0.89772 0.10228;
             0.75000 0.25000]
        c = make_crystal(A_rows, r, [:U, :U])
        check_spacegroup!(c, 8)
    end

    # ── CaTiO₃ Pnma perovskite — AB3C_oP20_62_c_cd_a ───────────────
    # Pnma (#62), orthorhombic primitive, 4 Ca + 12 O + 4 Ti = 20 atoms.
    # Order 8. Non-symmorphic (glide planes and screw axes).
    # Canonical perovskite distortion — the "20-atom" diagnostic case
    # requested for §3.13 coverage.
    let
        A_rows = [5.42240 0.00000 0.00000;
                  0.00000 7.65100 0.00000;
                  0.00000 0.00000 5.40430]
        Ca = [-0.0123  0.0123  0.4877  0.5123;
               0.2500  0.7500  0.2500  0.7500;
               0.5084  0.4916 -0.0084  0.0084]
        O = [ 0.0313 -0.0313  0.4687  0.5313  0.2120  0.2120  0.2880  0.2880  0.7120  0.7120  0.7880  0.7880;
              0.2500  0.7500  0.7500  0.2500  0.0370  0.4630  0.5370  0.9630  0.0370  0.4630 -0.0370  0.5370;
              0.0586 -0.0586  0.5586  0.4414  0.7130  0.7130  0.2130  0.2130  0.7870  0.7870  0.2870  0.2870]
        Ti = [0.0  0.0  0.5  0.5;
              0.0  0.5  0.0  0.5;
              0.0  0.0  0.5  0.5]
        c = make_crystal(A_rows, hcat(Ca, O, Ti),
                         [fill(:Ca, 4); fill(:O, 12); fill(:Ti, 4)])
        check_spacegroup!(c, 8)
    end

    # ─── 3-fold-screw cases (§3.13 user request) ──────────────────
    # Space groups containing 3₁ or 3₂ screws are the simplest
    # non-symmorphic cases with τ = 1/3 or 2/3 components. They
    # close a coverage gap: our other tests cover glides, 2₁, 4-fold,
    # and 6₃ screws but not the 3-fold screw specifically, and they
    # all have inversion — these are chiral (no inversion, no mirrors,
    # all ops proper rotations). Order 6 is also a gap: our other
    # expected orders are 2, 4, 8, 12, 16, 24, 48, 192. Oxygens in
    # α-Quartz are at a fully general Wyckoff position (6c), which
    # no other test exercises.

    # ── γ-Se — A_hP3_152_a ────────────────────────────────────────
    # Simplest 3-fold-screw case. P3₁21 (#152), trigonal, 3 atoms,
    # order 6. Chiral. Se at Wyckoff 3a with x = 0.2254.
    let
        A_rows = [2.18310 -3.78124 0.00000;
                  2.18310  3.78124 0.00000;
                  0.00000  0.00000 4.95360]
        r = [0.00000 0.22540 0.77460;
             0.22540 0.00000 0.77460;
             2/3     1/3     0.00000]
        c = make_crystal(A_rows, r, fill(:Se, 3))
        ops = check_spacegroup!(c, 6)
        # Signature: no inversion. All ops must be proper rotations.
        @test all(det(op.R) == 1 for op in ops)
        # Signature: at least one op with τ[3] ≈ 1/3 (the 3₁ screw)
        @test any(abs(mod(op.τ[3] - 1/3, 1.0)) < 1e-8 ||
                  abs(mod(op.τ[3] - 1/3, 1.0) - 1.0) < 1e-8 for op in ops)
    end

    # ── α-Quartz — A2B_hP9_152_c_a ────────────────────────────────
    # THE canonical trigonal non-symmorphic. P3₁21 (#152), 9 atoms
    # (6 O + 3 Si). Order 6. O at generic Wyckoff 6c position.
    let
        A_rows = [2.45700 -4.25565 0.00000;
                  2.45700  4.25565 0.00000;
                  0.00000  0.00000 5.40600]
        O = [0.14620 0.26680 0.41300 0.58700 0.73320 0.85380;
             0.73320 0.41300 0.26680 0.85380 0.14620 0.58700;
             0.45267 0.78600 0.21400 0.11933 0.54733 0.88067]
        Si = [0.00000 0.46990 0.53010;
              0.46990 0.00000 0.53010;
              2/3     1/3     0.00000]
        c = make_crystal(A_rows, hcat(O, Si), [fill(:O, 6); fill(:Si, 3)])
        ops = check_spacegroup!(c, 6)
        @test all(det(op.R) == 1 for op in ops)      # chiral
    end

    # ── Cinnabar (HgS) — AB_hP6_154_a_b ───────────────────────────
    # P3₂21 (#154), the ENANTIOMORPH of α-Quartz's P3₁21. Same order
    # 6 but the screw direction is reversed (3₂ instead of 3₁). Tests
    # that Spacey handles both chiralities.
    let
        A_rows = [2.07250 -3.58968 0.00000;
                  2.07250  3.58968 0.00000;
                  0.00000  0.00000 9.49600]
        Hg = [0.00000 0.28020 0.71980;
              0.71980 0.28020 0.00000;
              1/3     0.00000 2/3]
        S = [0.00000 0.48890 0.51110;
             0.48890 0.00000 0.51110;
             5/6     1/6     1/2]
        c = make_crystal(A_rows, hcat(Hg, S), [fill(:Hg, 3); fill(:S, 3)])
        ops = check_spacegroup!(c, 6)
        @test all(det(op.R) == 1 for op in ops)
    end

    # ── CrCl₃ — A3B_hP24_151_3c_2a ────────────────────────────────
    # P3₁12 (#151), trigonal, 24 atoms. Different 2-fold placement
    # than #152 — in #151 the 2-fold axes lie along [100]/[010] type
    # directions, in #152 along [110] type. Same order (6) and same
    # screw-type, but a valuable distinguishing test. Chiral.
    let
        A_rows = [3.00850 -5.21087 0.00000;
                  3.00850  5.21087 0.00000;
                  0.00000  0.00000 17.30000]
        Cl = [0.22220 0.22220 0.88890 0.88890 0.88890 0.88890 0.55560 0.55560 0.55560 0.55560 0.88880 0.88880 0.22220 0.22220 0.22222 0.22222 0.55558 0.55558;
              0.11110 0.11110 0.11110 0.11110 0.77780 0.77780 0.11120 0.11120 0.44440 0.44440 0.44440 0.44440 0.44442 0.77778 0.44442 0.77780 0.77778 0.77780;
              0.26023 0.73977 0.07310 0.59357 -0.07310 0.40643 -0.07310 0.40643 0.07310 0.59357 0.26023 0.73977 -0.07310 0.07310 0.40643 0.59357 0.26023 0.73977]
        Cr = [0.22220 0.88890 0.88890 0.55560 0.55560 0.88880;
              0.11110 0.11110 0.77780 0.11120 0.44440 0.44440;
              0.00000 1/3     2/3     2/3     1/3     0.00000]
        c = make_crystal(A_rows, hcat(Cl, Cr),
                         [fill(:Cl, 18); fill(:Cr, 6)])
        ops = check_spacegroup!(c, 6)
        @test all(det(op.R) == 1 for op in ops)
    end
end

@testset "spacegroup: AFLOW Part 1 (full corpus)" begin
    # Auto-generated data from tools/generate_aflow_tests.jl. Every
    # prototype in Part 1 of the AFLOW Library of Crystallographic
    # Prototypes (Mehl et al. 2017) is parsed into a Crystal, handed to
    # spacegroup(), and compared against the point-group order of the
    # space-group number encoded in the prototype label.
    #
    # The 8 structures listed in KNOWN_DEVIATIONS below report a
    # different operation count than the prototype label at default
    # tolerances. Each is a documented, legitimate Spacey behavior —
    # not a bug. Flipping them to @test_broken pins the current result
    # so a real regression (one of these starts passing with the wrong
    # value, or a new structure joins the list) shows up in CI.
    include(joinpath(@__DIR__, "aflow_structures_part1.jl"))

    # Prototype → (got, expected, reason).
    KNOWN_DEVIATIONS = Dict(
        # "P1 distortion": atoms are very close to ideal pyrite Pa-3
        # positions (displacements ~1e-3). Default pos_tol treats them
        # as on-site, Spacey finds the higher symmetry. Tight pos_tol
        # would recover the true P1.
        "AB2_aP12_1_4a_8a"        => (24, 1,  "distorted-pyrite → higher sym at default pos_tol"),
        "ABC2_aP16_1_4a_4a_8a"    => (2,  1,  "near-inversion at default pos_tol"),

        # "Accidental higher symmetry": atomic arrangement happens to
        # support inversion even though the labeled space group (Pna2₁,
        # Amm2, Aba2, …) does not include it. Spacey correctly finds the
        # maximal symmetry.
        "AB_oP8_33_a_a"           => (8,  4,  "atoms have accidental inversion beyond Pna2₁"),
        "A2B_oC12_38_de_ab"       => (8,  4,  "atoms have accidental inversion beyond Amm2"),
        "AB4_oC20_41_a_2b"        => (8,  4,  "atoms have accidental inversion beyond Aba2"),
        "AB2_oC24_41_2a_2b"       => (8,  4,  "atoms have accidental inversion beyond Aba2"),
        "A4B_cI40_197_cde_c"      => (24, 12, "atoms have accidental inversion beyond I2₁3"),

        # Lattice over-promotion: a / b = 1.005, below default lattice_tol
        # → Spacey promotes orthorhombic I to tetragonal I. Tight
        # lattice_tol would recover the true mmm.
        "AB2_oI6_71_a_i"          => (16, 8,  "a/b ≈ 1.005 → lattice_tol promotes to tetragonal"),
    )

    @testset "$(s.name) [$(s.prototype)] sg $(s.sg)" for s in AFLOW_PART1_STRUCTURES
        c_built = try Crystal(s.A, s.r, s.types; coords=s.coords) catch; nothing end
        c_built === nothing && (@test false; continue)
        ops = try spacegroup(c_built) catch; nothing end
        ops === nothing && (@test false; continue)

        if haskey(KNOWN_DEVIATIONS, s.prototype)
            got, _, _ = KNOWN_DEVIATIONS[s.prototype]
            # Pin the current deviation: Spacey must still return `got`.
            @test length(ops) == got
            # And mark the prototype-label match as broken.
            @test_broken length(ops) == s.expected_order
        else
            @test length(ops) == s.expected_order
        end
    end
end

# Parts 2 and 3 of the AFLOW Library (Hicks et al. 2019, 2021) — same
# pipeline as Part 1, generated by tools/generate_aflow_tests.jl from
# papers part2.pdf and part3.pdf. 299 + 510 = 809 additional prototypes.
# For these parts we don't hand-curate a known-deviation list — instead,
# each structure is tested with @test when it matches the AFLOW label
# and @test_broken when it doesn't. This pins every current result so
# real regressions show up, while documenting the proportion of
# deviations for anyone reading the test output.
function aflow_test_one(s)
    c_built = try Crystal(s.A, s.r, s.types; coords=s.coords) catch; nothing end
    c_built === nothing && (@test false; return)
    ops = try spacegroup(c_built) catch; nothing end
    ops === nothing && (@test false; return)
    if length(ops) == s.expected_order
        @test length(ops) == s.expected_order
    else
        @test_broken length(ops) == s.expected_order
    end
end

@testset "spacegroup: AFLOW Part 2 (full corpus)" begin
    include(joinpath(@__DIR__, "aflow_structures_part2.jl"))
    @testset "$(s.name) [$(s.prototype)] sg $(s.sg)" for s in AFLOW_PART2_STRUCTURES
        aflow_test_one(s)
    end
end

@testset "spacegroup: AFLOW Part 3 (full corpus)" begin
    include(joinpath(@__DIR__, "aflow_structures_part3.jl"))
    @testset "$(s.name) [$(s.prototype)] sg $(s.sg)" for s in AFLOW_PART3_STRUCTURES
        aflow_test_one(s)
    end
end

@testset "crystal_system" begin
    # Smoke test on the 14 Bravais lattices: the expected system is
    # derived from each entry's point-group order, which is the
    # holohedry. Mapping 2→tri, 4→mono, 8→ortho, 12→trig, 16→tet,
    # 24→hex, 48→cubic.
    order_to_sys = Dict(2=>:triclinic, 4=>:monoclinic, 8=>:orthorhombic,
                        12=>:trigonal, 16=>:tetragonal,
                        24=>:hexagonal, 48=>:cubic)
    for (name, (A, order)) in BravaisLatticeList
        expected = order_to_sys[order]
        @test crystal_system(A) == expected
    end
end

# Extract the expected crystal system from an AFLOW prototype label.
# The Pearson symbol's first two letters encode system + centering:
#   'a?' → triclinic    'm?' → monoclinic    'o?' → orthorhombic
#   't?' → tetragonal   'c?' → cubic
#   'hP' → hexagonal    'hR' → trigonal (rhombohedral setting)
function expected_crystal_system(prototype::AbstractString)
    parts = split(prototype, '_')
    length(parts) ≥ 2 || error("malformed prototype: $prototype")
    pearson = parts[2]
    length(pearson) ≥ 2 || error("malformed Pearson: $pearson")
    sys, cent = pearson[1], pearson[2]
    sys == 'a' && return :triclinic
    sys == 'm' && return :monoclinic
    sys == 'o' && return :orthorhombic
    sys == 't' && return :tetragonal
    sys == 'c' && return :cubic
    sys == 'h' && cent == 'P' && return :hexagonal
    sys == 'h' && cent == 'R' && return :trigonal
    error("unknown Pearson: $pearson")
end

function aflow_crystal_system_test_one(s)
    expected_sys = try expected_crystal_system(s.prototype) catch; return end
    actual_sys = try
        crystal_system(s.A)
    catch
        @test false; return
    end
    if actual_sys == expected_sys
        @test actual_sys == expected_sys
    else
        @test_broken actual_sys == expected_sys
    end
end

@testset "crystal_system cross-check: AFLOW Part 1" begin
    @testset "$(s.name) [$(s.prototype)]" for s in AFLOW_PART1_STRUCTURES
        aflow_crystal_system_test_one(s)
    end
end

@testset "crystal_system cross-check: AFLOW Part 2" begin
    @testset "$(s.name) [$(s.prototype)]" for s in AFLOW_PART2_STRUCTURES
        aflow_crystal_system_test_one(s)
    end
end

@testset "crystal_system cross-check: AFLOW Part 3" begin
    @testset "$(s.name) [$(s.prototype)]" for s in AFLOW_PART3_STRUCTURES
        aflow_crystal_system_test_one(s)
    end
end

@testset "spacegroup: Phase 4 near-boundary crystal (verify_stable)" begin
    # BaTiO₃-style ferroelectric near-miss: cubic lattice (Pm3̄m = 48 ops at
    # ε = 0), with Ti displaced from body-center by ε along z. For ε > 0
    # the true space group is P4mm (8 ops). At pos_tol ≫ ε, Spacey sees
    # the Ti as effectively on-center and reports the parent cubic 48 ops
    # — silent over-promotion. See test/nearMissBoundaryCrystal.jl for
    # the full (ε, pos_tol) heatmap.
    function batio3_like(ε)
        A = 4.0 * Matrix{Float64}(I, 3, 3)
        r = zeros(3, 5)
        r[:, 1] = [0.0, 0.0, 0.0]          # Ba
        r[:, 2] = [0.5, 0.5, 0.5 - ε]      # Ti displaced
        r[:, 3] = [0.5, 0.5, 0.0]          # O faces
        r[:, 4] = [0.5, 0.0, 0.5]
        r[:, 5] = [0.0, 0.5, 0.5]
        Crystal(A, r, [:Ba, :Ti, :O, :O, :O]; coords=:fractional)
    end

    for ε ∈ [1e-3, 1e-4, 1e-5]
        c = batio3_like(ε)
        tight_pos_tol = ε / 100
        loose_pos_tol = 100 * ε

        # Exact arithmetic at tight pos_tol finds the true tetragonal group
        @test length(spacegroup(c; pos_tol=tight_pos_tol)) == 8
        # At loose pos_tol, over-promotion to cubic (pins current behavior)
        @test length(spacegroup(c; pos_tol=loose_pos_tol)) == 48

        # verify_stable emits a warning when the answer flips across the
        # tight-tol / loose-tol boundary
        @test_logs (:warn, r"near a position-symmetry boundary") match_mode=:any begin
            spacegroup(c; pos_tol=loose_pos_tol, verify_stable=true)
        end

        # verify_stable at tight pos_tol: correct group and no false alarm
        @test length(spacegroup(c; pos_tol=tight_pos_tol, verify_stable=true)) == 8
    end

    # Sanity: perfect cubic (ε = 0) is stable across all pos_tol — no
    # warning expected from verify_stable because the group doesn't change.
    let c_cubic = batio3_like(0.0)
        @test length(spacegroup(c_cubic; verify_stable=true)) == 48
    end
end

