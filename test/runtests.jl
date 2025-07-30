using Spacey
using Test
using LinearAlgebra
using MinkowskiReduction
    
# Timing tests are not consistent on github virtual machines.
# Timing tests are only checked locally, see runTimingTests.jl

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
    a,b,c,ops = snapToSymmetry_SVD(u,v,w,ops)
    @test det([a b c])≈det([a1 a2 a3])
    @test norm(a)≈norm(b)≈norm(c)
    @test length(ops)==48

    # Simple tetragonal example of snap
    a1 = [1+.01,0,0]; a2 = [0.,1-.01,0]; a3 = [0,0,1.5];
    u,v,w = minkReduce(a1,a2,a3)
    ops,_ = pointGroup_robust(u,v,w;tol=5e-2)
    @test isagroup(ops)
    a,b,c,ops = snapToSymmetry_SVD(u,v,w,ops)
    @test length(ops)==16
    @test norm(a)≈norm(b)
    @test det([a b c])≈det([a1 a2 a3])
    
    # Simple orthorhombic example of snap
    println("Orthorhombic example of snap")
    a1 = [1,0.01,-.005]; a2 = [0.001,2,-0.01]; a3 = [0.02,-0.003,1.5];
    u,v,w = minkReduce(a1,a2,a3)
    ops,_ = pointGroup_robust(u,v,w;tol=5e-2)
    a,b,c,ops = snapToSymmetry_SVD(u,v,w,ops)
    @test length(ops)==8
    @test det([a b c])≈det([a1 a2 a3])
    @test isagroup(ops)

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

