using Documenter
using Spacey

DocMeta.setdocmeta!(Spacey, :DocTestSetup, :(using Spacey, LinearAlgebra); recursive=true)

makedocs(
    sitename = "Spacey.jl",
    authors  = "Gus Hart and contributors",
    repo     = Remotes.GitHub("glwhart", "Spacey.jl"),
    format   = Documenter.HTML(
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical  = "https://glwhart.github.io/Spacey.jl",
        assets     = String[],
    ),
    modules  = [Spacey],
    # Only enforce that *exported* symbols appear in @docs/@autodocs blocks.
    # Spacey has a handful of internal helpers with docstrings (e.g. avgVecOverOps)
    # that we don't want to surface in the user-facing API; checkdocs=:exports
    # lets them keep their docstrings without forcing them into the reference.
    checkdocs = :exports,
    pages = [
        "Home" => "index.md",
        "Tutorials" => [
            "tutorials/index.md",
            "tutorials/01-first-pointgroup.md",
            "tutorials/02-first-spacegroup.md",
        ],
        "How-to guides" => [
            "how-to/index.md",
            "how-to/construct-a-crystal.md",
            "how-to/find-pointgroup.md",
            "how-to/find-spacegroup.md",
            "how-to/classify-bravais.md",
            "how-to/handle-noisy-data.md",
            "how-to/detect-tolerance-dependence.md",
            "how-to/compose-and-apply-ops.md",
            "how-to/snap-to-symmetry.md",
        ],
        "Reference" => [
            "reference/index.md",
            "reference/crystals.md",
            "reference/point-groups.md",
            "reference/space-groups.md",
            "reference/snap.md",
            "reference/helpers.md",
        ],
        "Explanation" => [
            "explanation/index.md",
            "explanation/why-minkowski.md",
            "explanation/tolerances.md",
            "explanation/over-promotion.md",
            "explanation/canonicalising-tau.md",
            "explanation/crystal-system-vs-bravais.md",
            "explanation/algorithm-overview.md",
            "explanation/validation-strategy.md",
        ],
    ],
)

deploydocs(
    repo      = "github.com/glwhart/Spacey.jl",
    devbranch = "main",
    devurl    = "dev",
    target    = "build",
    branch    = "gh-pages",
    versions  = ["stable" => "v^"],
)
