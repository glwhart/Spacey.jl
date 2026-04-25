# Spacey.jl

Spacey.jl finds the **point group** of a Bravais lattice and the **space group** of a crystal, with algorithms specifically designed to be robust to finite-precision floating-point error in the input.

This documentation is organised into four sections following the [Diátaxis](https://diataxis.fr/) framework. Pick the one that matches what you need right now:

| If you're …                    | Go to …                              |
|---                             |---                                   |
| New to Spacey                  | [Tutorials](tutorials/index.md)      |
| Solving a specific task        | [How-to guides](how-to/index.md)     |
| Looking up a function          | [Reference](reference/index.md)      |
| Trying to understand the why   | [Explanation](explanation/index.md)  |

## 30-second smoke test

```julia
julia> using Pkg; Pkg.add(url="https://github.com/glwhart/Spacey.jl")

julia> using Spacey, LinearAlgebra

julia> LG, G = pointGroup_robust([1.0,0,0], [0,1.0,0], [0,0,1.0]);

julia> length(G)   # 48 — the cubic point group
48
```

If that runs and prints `48`, the package is loaded and working.

!!! note "Documentation status"
    These docs are being built out in phases. The skeleton is live; content lands in subsequent phases (see `documentation_plan.md` in the repo root).
