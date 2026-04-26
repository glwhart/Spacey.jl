using Spacey
using MinkowskiReduction
using LinearAlgebra
using Printf

# Near-miss boundary diagnostic. Three panels, each a (ε, tol) heatmap of the
# group size reported by pointGroup_robust, marked against the truth and against
# what verify_stable (tight_tol = tol/1000) would detect.
#
# Legend:
#   " " correct, stable
#   "~" correct at user tol, but verify_stable would false-alarm (tight differs)
#   "*" over-promoted, verify_stable fires (tight ≠ loose)
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

function panel(label, basis_fn, truth_fn, εs, tols)
    println()
    println("="^80)
    println(label)
    println("="^80)
    print(@sprintf("%12s", "ε \\ tol"))
    for t ∈ tols
        print(@sprintf("%10.0e", t))
    end
    println()
    println("-"^(12 + 10*length(tols)))
    for ε ∈ εs
        A = basis_fn(ε)
        u, v, w = minkReduce(eachcol(A)...)[1:3]
        truth = truth_fn(ε)
        print(@sprintf("%12.0e", ε))
        for t ∈ tols
            try
                n_loose = length(Spacey.pointGroup_robust(u, v, w; tol=t)[1])
                n_tight = length(Spacey.pointGroup_robust(u, v, w; tol=t/1000)[1])
                print(@sprintf("%9d%s", n_loose, mark(n_loose, n_tight, truth)))
            catch e
                print(@sprintf("%10s", "err"))
            end
        end
        println()
    end
    println()
    print("truth by ε: ")
    for ε ∈ εs
        print(@sprintf("%d ", truth_fn(ε)))
    end
    println()
end

εs   = [0.0, 1e-8, 1e-6, 1e-5, 1e-4, 1e-3, 3e-3, 1e-2, 3e-2, 1e-1]
tols = [1e-10, 1e-8, 1e-6, 1e-4, 1e-3, 1e-2, 1e-1]

# Panel 1: Tetragonal → Cubic.  c/a = 1+ε, truth = 16 for ε>0, 48 for ε=0.
panel(
    "Panel 1: Tetragonal → Cubic   A = diag(1, 1, 1+ε)",
    ε -> [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0+ε],
    ε -> ε == 0 ? 48 : 16,
    εs, tols,
)

# Panel 2: Orthorhombic → Tetragonal.  b/a = 1+ε (c unrelated). truth = 8 for ε>0, 16 for ε=0.
panel(
    "Panel 2: Orthorhombic → Tetragonal   A = diag(1, 1+ε, 1.3)",
    ε -> [1.0 0.0 0.0; 0.0 1.0+ε 0.0; 0.0 0.0 1.3],
    ε -> ε == 0 ? 16 : 8,
    εs, tols,
)

# Panel 3: Rhombohedral → Simple Cubic.
# b_i = ê_i + ε(others); three equal-length vectors with equal angles.
# At ε=0 it's cubic (48). For ε>0 it's rhombohedral (12).
panel(
    "Panel 3: Rhombohedral → Cubic   b_i = ê_i + ε·(1-ê_i)",
    ε -> [1.0 ε ε; ε 1.0 ε; ε ε 1.0],
    ε -> ε == 0 ? 48 : 12,
    εs, tols,
)

println()
println("Summary:  \" \"=correct  \"~\"=false-alarm  \"*\"=detected wrong  \"!\"=silent wrong  \"?\"=under-report")
