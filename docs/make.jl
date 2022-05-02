using Documenter
using Spacey

makedocs(
    sitename = "Spacey",
    format = Documenter.HTML(),
    modules = [Spacey]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
