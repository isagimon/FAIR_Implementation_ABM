using Pkg

# Build docs in an isolated environment
Pkg.activate(@__DIR__)
Pkg.develop(PackageSpec(path=joinpath(@__DIR__, "..")))
Pkg.instantiate()

using Documenter
using FAIR_Implementation_ABM

makedocs(
    sitename = "FAIR_Implementation_ABM",
    format = Documenter.HTML(),
    modules = [FAIR_Implementation_ABM],
    pages = [
        "Home" => "index.md"
    ],
)

