using Spacey
using Test
using BenchmarkTools

# These tests were originally part of the set that run automatically, but the difference between
# "fast" and "slow" algorithms is much smaller on the github virtual machines. So keeping these
# timing tests separate so I can compare my own machine times when updated algorithms
@testset "Spacey.jl" begin
    Δ=0.1 # Fractional change in speedup that triggers an error
    #ideal hex lattice, 60° between basal plane vectors
    u = [1, 0, 0]
    v = [0.5, √3 / 2, 0]
    w = [0, 0, √(8 / 3)]
    @test length(pointGroup_simple(u, v, w)) == 24
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
    # Add a timing check against basic and "fast"; make sure fast is fast
    # hexagonal case
    slow = @belapsed pointGroup_simple($g, $h, $i)
    fast = @belapsed pointGroup_fast($g, $h, $i)
    println("actual speed up: ", round(slow/fast,digits=2))
    @test isapprox(slow / fast, 27, rtol=Δ)
    # cubic case
    slow = @belapsed pointGroup_simple($d, $e, $f)
    fast = @belapsed pointGroup_fast($d, $e, $f)
    println("actual speed up: ", round(slow / fast,digits=2))
    @test isapprox(slow / fast, 53, rtol=Δ)
    # tetragonal case
    d, e, f = threeDrotation([1.1, 0, 0], [0, 1, 0], [0, 0, 1], π / 3, π / 5, π / 7)
    slow = @belapsed pointGroup_simple($d, $e, $f)
    fast = @belapsed pointGroup_fast($d, $e, $f)
    println("actual speed up: ", round(slow / fast,digits=2))
    @test isapprox(slow / fast, 61, rtol=Δ)
    # orthorhombic case
    d, e, f = threeDrotation([1.1, 0, 0], [0, 0.9, 0], [0, 0, 1], π / 3, π / 5, π / 7)
    slow = @belapsed pointGroup_simple($d, $e, $f)
    fast = @belapsed pointGroup_fast($d, $e, $f)
    println("actual speed up: ", round(slow / fast,digits=2))
    @test isapprox(slow / fast, 64, rtol=Δ)
    u = [1, 1, 2]
    v = [1, 2, 1]
    w = [2, 1, 1]
    u, v, w = threeDrotation(u, v, w, π / 3, π / 5, π / 7)
    @test length(pointGroup_fast(u, v, w)) == 12
    # Faster version
    slow = @belapsed pointGroup_simple($g, $h, $i)
    fast = @belapsed pointGroup_fast($g, $h, $i)
    println("actual speed up: ", round(slow/fast,digits=2))
    @test isapprox(slow / fast, 95, rtol=Δ)
    # cubic case
    slow = @belapsed pointGroup_simple($d, $e, $f)
    fast = @belapsed pointGroup_fast($d, $e, $f)
    println("actual speed up: ", round(slow / fast,digits=2))
    @test isapprox(slow / fast, 307, rtol=Δ)
    # tetragonal case
    d, e, f = threeDrotation([1.1, 0, 0], [0, 1, 0], [0, 0, 1], π / 3, π / 5, π / 7)
    slow = @belapsed pointGroup_simple($d, $e, $f)
    fast = @belapsed pointGroup_fast($d, $e, $f)
    println("actual speed up: ", round(slow / fast,digits=2))
    @test isapprox(slow / fast, 180, rtol=Δ)
    # orthorhombic case
    d, e, f = threeDrotation([1.1, 0, 0], [0, 0.9, 0], [0, 0, 1], π / 3, π / 5, π / 7)
    slow = @belapsed pointGroup_simple($d, $e, $f)
    fast = @belapsed pointGroup_fast($d, $e, $f)
    println("actual speed up: ", round(slow / fast,digits=2))
    @test isapprox(slow / fast, 307, rtol=Δ)
end
