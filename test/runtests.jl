using Spacey
using Test
using LinearAlgebra
using MinkowskiReduction

@testset "Spacey.jl" begin
    Δ = 0.2 # Fractional change in speedup between fast and slow algorithms that triggers an error
    #ideal hex lattice, 60° between basal plane vectors
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
    
    # Timing tests are not consistent on github virtual machines.
    # Timing tests are only checked locally, see runTimingTests.jl

    # Add a timing check against basic and "fast" make sure fast is fast
    # hexagonal case
    # slow = @belapsed pointGroup_simple($g, $h, $i)
    # fast = @belapsed pointGroup_fast($g, $h, $i)
    # @test isapprox(slow / fast, 154, rtol=Δ)
    # # cubic case
    # slow = @belapsed pointGroup_simple($d, $e, $f)
    # fast = @belapsed pointGroup_fast($d, $e, $f)
    # @test isapprox(slow / fast, 100, rtol=Δ)
    # # tetragonal case
    # d, e, f = threeDrotation([1.1, 0, 0], [0, 1, 0], [0, 0, 1], π / 3, π / 5, π / 7)
    # slow = @belapsed pointGroup_simple($d, $e, $f)
    # fast = @belapsed pointGroup_fast($d, $e, $f)
    # @test isapprox(slow / fast, 350, rtol=Δ)
    # # orthorhombic case
    # d, e, f = threeDrotation([1.1, 0, 0], [0, 0.9, 0], [0, 0, 1], π / 3, π / 5, π / 7)
    # slow = @belapsed pointGroup_simple($d, $e, $f)
    # fast = @belapsed pointGroup_fast($d, $e, $f)
    # @test isapprox(slow / fast, 20, rtol=Δ)
    # # Rhombohedral case
    u = [1, 1, 2]
    v = [1, 2, 1]
    w = [2, 1, 1]
    u, v, w = threeDrotation(u, v, w, π / 3, π / 5, π / 7)
    @test length(pointGroup_fast(u, v, w)) == 12

    # Cases of small tiny noise in input
    ϵ = 10.0e-10; u = [1, 0, ϵ]; v = [ϵ,1,0]; w = [0,0,1]
    @test length(pointGroup_fast(u,v,w))==48
    ϵ = 10.0e-8; u = [1, 0, ϵ]; v = [ϵ,1,0]; w = [0,0,1]
    @test length(pointGroup_fast(u,v,w))≠48
    # High aspect ratio cases
    a = 2^25; u = [1, 0, 0]; v = [0,a,0];  w = [0,0,1]
    @test length(pointGroup_fast(u,v,w))==16
    a = 2^26; u = [1, 0, 0]; v = [0,a,0];  w = [0,0,1]
    @test length(pointGroup_fast(u,v,w))≠16

    # Simple cubic example of snap function, ~1% noise
    a1 = [1+.01,0,0]; a2 = [0.,1-.01,0]; a3 = [0,0,1-.001];
    u,v,w = minkReduce(a1,a2,a3)
    ops,_ = pointGroup_robust(u,v,w)
    @test isagroup(ops)
    a,b,c,ops = snapToSymmetry_SVD(u,v,w,ops)
    @test det([a b c])≈det([a1 a2 a3])
    @test norm(a)≈norm(b)≈norm(c)
    @test length(ops)==48

    # Simple tetragonal example of snap
    a1 = [1+.01,0,0]; a2 = [0.,1-.01,0]; a3 = [0,0,1.5];
    u,v,w = minkReduce(a1,a2,a3)
    ops,_ = pointGroup_robust(u,v,w)
    @test isagroup(ops)
    a,b,c,ops = snapToSymmetry_SVD(u,v,w,ops)
    @test length(ops)==16
    @test norm(a)≈norm(b)
    @test det([a b c])≈det([a1 a2 a3])
    
    # Simple orthorhombic example of snap
    a1 = [1,0.01,-.05]; a2 = [0.001,2,-0.01]; a3 = [0.02,-0.003,1.5];
    u,v,w = minkReduce(a1,a2,a3)
    ops,_ = pointGroup_robust(u,v,w)
    a,b,c,ops = snapToSymmetry_SVD(u,v,w,ops)
    @test length(ops)==8
    @test det([a b c])≈det([a1 a2 a3])
    @test isagroup(ops)

    a1 = [1/16 - .0001, .001, .0001]
    a2 = [-0.001, 16-.0001, -.0001]
    a3 = [-.0001, .0001, 1.51]
    ops,_ = pointGroup_robust(a1,a2,a3)
    @test isagroup(ops)
end
