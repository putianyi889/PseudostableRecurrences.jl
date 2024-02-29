using Documenter
using PseudostableRecurrences

makedocs(
    sitename = "PseudostableRecurrences",
    format = Documenter.HTML(),
    modules = [PseudostableRecurrences]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "https://github.com/putianyi889/PseudostableRecurrences.jl.git",
    devurl = "dev",
    versions = ["stable" => "v^", "v#.#.#", devurl => devurl]
)
