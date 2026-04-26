using Spacey
using MinkowskiReduction
using LinearAlgebra
using Printf

# Crystal-level near-miss boundary diagnostic. Parallel to
# test/nearMissBoundary.jl (which did the point-group case). Shows how
# the reported space-group operation count depends on (displacement ε,
# pos_tol) for crystals that approach a higher-symmetry parent via a
# small atomic displacement.
#
# Panel 1 — BaTiO₃-style ferroelectric displacement.
#   Cubic Pm3̄m parent (space group 221, order 48) with Ba at corners,
#   Ti at body center, O at face centers. Apply a small z-displacement
#   δ to Ti: (½, ½, ½) → (½, ½, ½ − δ). For δ > 0 the true space group
#   is P4mm (order 8). Pos_tol too loose → Spacey reports 48 silently.
#
# Legend:
#   " " correct, stable
#   "~" correct at user pos_tol, but verify_stable would false-alarm
#   "*" over-promoted, verify_stable catches (tight ≠ loose)
#   "!" over-promoted, verify_stable silent (tight == loose, both wrong)
#   "?" under-reported

function mark(n_loose, n_tight, truth)
    if n_loose == truth
        return n_loose == n_tight ? " " : "~"
    elseif n_loose > truth
        return n_loose == n_tight ? "!" : "*"
    else
        return "?"
    end
end

function batio3_panel(εs, pos_tols)
    println("="^80)
    println("Panel: BaTiO₃-style ferroelectric Ti displacement")
    println("       A = diag(4, 4, 4) Å;  Ba (0,0,0);  Ti (½,½,½−ε);  O at faces")
    println("       truth: ε = 0 → cubic (48 ops);  ε > 0 → tetragonal (8 ops)")
    println("="^80)
    print(@sprintf("%12s", "ε \\ pos_tol"))
    for t ∈ pos_tols
        print(@sprintf("%11.0e", t))
    end
    println()
    println("-"^(12 + 11*length(pos_tols)))

    a = 4.0
    A = a * Matrix{Float64}(I, 3, 3)

    for ε ∈ εs
        r = zeros(3, 5)
        r[:, 1] = [0.0, 0.0, 0.0]                 # Ba
        r[:, 2] = [0.5, 0.5, 0.5 - ε]             # Ti (displaced)
        r[:, 3] = [0.5, 0.5, 0.0]                 # O face
        r[:, 4] = [0.5, 0.0, 0.5]                 # O face
        r[:, 5] = [0.0, 0.5, 0.5]                 # O face
        types = [:Ba, :Ti, :O, :O, :O]
        c = Crystal(A, r, types; coords=:fractional)
        truth = ε == 0 ? 48 : 8
        print(@sprintf("%12.0e", ε))
        for t ∈ pos_tols
            try
                n_loose = length(spacegroup(c; pos_tol=t))
                n_tight = length(spacegroup(c; pos_tol=t/1000))
                print(@sprintf("%10d%s", n_loose, mark(n_loose, n_tight, truth)))
            catch e
                print(@sprintf("%11s", "err"))
            end
        end
        println()
    end
    println()
end

# Displacement sweep (in fractional units of c = 4 Å, so ε = 0.01 = 0.04 Å)
εs = [0.0, 1e-5, 1e-4, 1e-3, 3e-3, 1e-2, 3e-2, 1e-1]
pos_tols = [1e-8, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1]

batio3_panel(εs, pos_tols)

println("Summary:  \" \"=correct  \"~\"=false-alarm  \"*\"=detected wrong  \"!\"=silent wrong  \"?\"=under-report")
