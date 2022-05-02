using Spacey
using Test
using BenchmarkTools

@testset "Spacey.jl" begin
    #ideal hex lattice, 60° between basal plane vectors
    u = [1, 0, 0]
    v = [0.5, √3 / 2, 0]
    w = [0, 0, √(8 / 3)]
    @test length(pointGroup_basic(u, v, w)) == 24
    #ideal hex lattice, 120° between basal plane vectors
    a = -u + v
    b = v
    c = w
    @test length(pointGroup_basic(a, b, c)) == 24
    # simple cubic lattice in unique orientation
    d, e, f = threeDrotation([1, 0, 0], [0, 1, 0], [0, 0, 1], π / 3, π / 5, π / 7)
    @test length(pointGroup_basic(d, e, f)) == 48
    # hexagonal lattice in unique orientation
    g, h, i = threeDrotation([1, 0, 0], [0.5, √3 / 2, 0], [0, 0, √(8 / 3)], π / 7, π / 11, π / 3)
    @test length(pointGroup_basic(g, h, i)) == 24
    # Add a timing check against basic and "fast" make sure fast is fast
    # hexagonal case
    slow = @belapsed pointGroup_basic($g, $h, $i)
    fast = @belapsed pointGroup($g, $h, $i)
    @test slow / fast > 25
    # cubic case
    slow = @belapsed pointGroup_basic($d, $e, $f)
    fast = @belapsed pointGroup($d, $e, $f)
    @test slow / fast > 54
    # tetragonal case
    d, e, f = threeDrotation([1.1, 0, 0], [0, 1, 0], [0, 0, 1], π / 3, π / 5, π / 7)
    slow = @belapsed pointGroup_basic($d, $e, $f)
    fast = @belapsed pointGroup($d, $e, $f)
    @test slow / fast > 60
    # orthorhombic case
    d, e, f = threeDrotation([1.1, 0, 0], [0, 0.9, 0], [0, 0, 1], π / 3, π / 5, π / 7)
    slow = @belapsed pointGroup_basic($d, $e, $f)
    fast = @belapsed pointGroup($d, $e, $f)
    @test slow / fast > 62
    u = [1, 1, 2]
    v = [1, 2, 1]
    w = [2, 1, 1]
    u, v, w = threeDrotation(u, v, w, π / 3, π / 5, π / 7)
    @test length(pointGroup(u, v, w)) == 12
end
