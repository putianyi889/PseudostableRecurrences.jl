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
    repo = "github.com/putianyi889/PseudostableRecurrences.jl.git",
    devurl = "dev",
    versions = ["stable" => "v^", 
        "v#.#.#", 
        "v0.0.3" => "v0.0.3-doc1", 
        "v0.0.2" => "v0.0.2-doc2", 
        "v0.0.1" => "v0.0.1-doc2", 
        "dev" => "dev"]
)
