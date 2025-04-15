using Documenter
using FAIR_Implementation_ABM

makedocs(
    sitename = "FAIR_Implementation_ABM",
    format = Documenter.HTML(),
    modules = [FAIR_Implementation_ABM],
    pages = [
        "Home" => "index.md"
    ]
)

