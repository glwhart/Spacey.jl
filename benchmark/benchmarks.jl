# Benchmark suite for Spacey.jl, driven by PkgBenchmark.jl.
#
# This file is the single source of truth for which inputs to time. It is
# *not* run as part of `Pkg.test()` — execute it via `PkgBenchmark.judge`
# (see "How to use" below) so timings are compared between two commits on
# the same machine. Absolute timings on shared CI are unreliable; ratio-based
# comparison kills most of the noise.
#
# Coverage:
#   - `pointGroup`: matrix-form and three-vector-form on canonical inputs
#     (cubic, hexagonal, the FCC primitive, a high-aspect-ratio basis).
#   - `spacegroup`: representative crystal sizes and types — NaCl primitive
#     (2 atoms), diamond primitive (2 atoms, full O_h), BaTiO₃ (5 atoms,
#     ferroelectric), HCP cadmium (2 atoms, hex).
#   - `crystal_system`: wraps `pointGroup`, included as a regression check
#     for that delegation.
#   - `isSpacegroupOp`: hot inner-loop primitive used by `spacegroup`.
#
# Selection rationale: each benchmark is small enough to run in milliseconds
# (the suite finishes in seconds), and the inputs cover both the
# small-symmetry and large-symmetry regimes.
#
# How to use (manually, before a release):
#
#   julia --project=benchmark
#   julia> using PkgBenchmark
#   julia> result = judge("Spacey", "main", "HEAD")  # compare branches
#   julia> export_markdown("benchmark.md", result)   # human-readable diff
#
# Or, on a single commit:
#
#   julia> result = benchmarkpkg("Spacey")
#   julia> export_markdown("benchmark.md", result)
#
# `judge` is the load-bearing API: it runs both commits on the same machine
# and reports relative speedups/slowdowns with the noise mostly cancelled.

using BenchmarkTools, Spacey, LinearAlgebra

const SUITE = BenchmarkGroup()

# --- pointGroup -----------------------------------------------------------
SUITE["pointGroup"] = BenchmarkGroup()

let A_cubic = Matrix{Float64}(I, 3, 3)
    SUITE["pointGroup"]["matrix, simple cubic (48 ops)"] =
        @benchmarkable pointGroup($A_cubic)
end

let u = [1.0, 0, 0], v = [-0.5, sqrt(3)/2, 0], w = [0.0, 0, sqrt(8/3)]
    SUITE["pointGroup"]["3-vector, ideal HCP (24 ops)"] =
        @benchmarkable pointGroup($u, $v, $w)
end

let A_fcc = [0.0 0.5 0.5; 0.5 0.0 0.5; 0.5 0.5 0.0]
    SUITE["pointGroup"]["matrix, FCC primitive (48 ops)"] =
        @benchmarkable pointGroup($A_fcc)
end

# High aspect ratio — exercises the volume-normalization path
let A_thin = [1.0 0 0; 0 1.0 0; 0 0 100.0]
    SUITE["pointGroup"]["matrix, AR=100 tetragonal (16 ops)"] =
        @benchmarkable pointGroup($A_thin)
end

# --- crystal_system -------------------------------------------------------
SUITE["crystal_system"] = BenchmarkGroup()

let A_cubic = Matrix{Float64}(I, 3, 3)
    SUITE["crystal_system"]["cubic"] = @benchmarkable crystal_system($A_cubic)
end

let A_tet = [1.0 0 0; 0 1.0 0; 0 0 1.5]
    SUITE["crystal_system"]["tetragonal"] = @benchmarkable crystal_system($A_tet)
end

# --- spacegroup -----------------------------------------------------------
SUITE["spacegroup"] = BenchmarkGroup()

let A_fcc = [0.0 0.5 0.5; 0.5 0.0 0.5; 0.5 0.5 0.0]
    c_NaCl = Crystal(A_fcc, [0.0 0.5; 0.0 0.5; 0.0 0.5], [:Na, :Cl]; coords=:fractional)
    SUITE["spacegroup"]["NaCl primitive (48 ops)"] = @benchmarkable spacegroup($c_NaCl)

    c_diamond = Crystal(A_fcc, [0.0 0.25; 0.0 0.25; 0.0 0.25], [:C, :C]; coords=:fractional)
    SUITE["spacegroup"]["diamond primitive (48 ops)"] = @benchmarkable spacegroup($c_diamond)
end

let A_hex = [1.0 -0.5 0.0; 0.0 sqrt(3)/2 0.0; 0.0 0.0 sqrt(8/3)]
    c_HCP = Crystal(A_hex, [0.0 1/3; 0.0 2/3; 0.0 0.5], [:Cd, :Cd]; coords=:fractional)
    SUITE["spacegroup"]["HCP Cd (24 ops)"] = @benchmarkable spacegroup($c_HCP)
end

let A_BTO = 4.0 * Matrix{Float64}(I, 3, 3), δ = 0.05/4
    r_BTO = [0.0 0.5 0.5 0.5 0.0; 0.0 0.5 0.5 0.0 0.5; 0.0 0.5-δ 0.0 0.5 0.5]
    c_BaTiO3 = Crystal(A_BTO, r_BTO, [:Ba, :Ti, :O, :O, :O]; coords=:fractional)
    # Default pos_tol over-promotes to cubic (48 ops); tight pos_tol resolves
    # the displacement to 8 ops. Both are useful as benchmarks since they
    # exercise different code paths through `isSpacegroupOp`.
    SUITE["spacegroup"]["BaTiO3 default pos_tol (48 ops, over-promoted)"] =
        @benchmarkable spacegroup($c_BaTiO3)
    SUITE["spacegroup"]["BaTiO3 tight pos_tol (8 ops, resolved)"] =
        @benchmarkable spacegroup($c_BaTiO3; pos_tol=1e-4)
end

# --- isSpacegroupOp -------------------------------------------------------
SUITE["isSpacegroupOp"] = BenchmarkGroup()

let A_fcc = [0.0 0.5 0.5; 0.5 0.0 0.5; 0.5 0.5 0.0],
    c_NaCl = Crystal(A_fcc, [0.0 0.5; 0.0 0.5; 0.0 0.5], [:Na, :Cl]; coords=:fractional),
    R_I = Matrix{Int}(I, 3, 3), τ_zero = zeros(3)

    SUITE["isSpacegroupOp"]["identity on NaCl (true)"] =
        @benchmarkable isSpacegroupOp($R_I, $τ_zero, $c_NaCl)

    τ_off = [0.5, 0.0, 0.0]
    SUITE["isSpacegroupOp"]["wrong shift on NaCl (false)"] =
        @benchmarkable isSpacegroupOp($R_I, $τ_off, $c_NaCl)
end

# --- verify_stable overhead -----------------------------------------------
SUITE["verify_stable"] = BenchmarkGroup()

let A_cubic = Matrix{Float64}(I, 3, 3)
    SUITE["verify_stable"]["pointGroup off"] =
        @benchmarkable pointGroup($A_cubic; verify_stable=false)
    SUITE["verify_stable"]["pointGroup on (re-runs at tol/1000)"] =
        @benchmarkable pointGroup($A_cubic; verify_stable=true)
end

let A_fcc = [0.0 0.5 0.5; 0.5 0.0 0.5; 0.5 0.5 0.0]
    c_NaCl = Crystal(A_fcc, [0.0 0.5; 0.0 0.5; 0.0 0.5], [:Na, :Cl]; coords=:fractional)
    SUITE["verify_stable"]["spacegroup off"] =
        @benchmarkable spacegroup($c_NaCl; verify_stable=false)
    SUITE["verify_stable"]["spacegroup on (re-runs at pos_tol/1000)"] =
        @benchmarkable spacegroup($c_NaCl; verify_stable=true)
end
