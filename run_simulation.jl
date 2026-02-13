# run_simulation.jl
# --------------------------------------------------
# Launches the FAIR_Implementation_ABM simulation.
#
# Usage:
#   julia run_simulation.jl
# --------------------------------------------------

using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

using FAIR_Implementation_ABM

# Run one simulation using the parameters in Input_Parameters.csv
FAIR_Implementation_ABM.run_simulation()
