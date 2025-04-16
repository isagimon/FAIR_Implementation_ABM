# run_simulation.jl
# --------------------------------------------------
# Launches the FAIR_Implementation_ABM simulation
# using Environment_and_Movement.jl
# --------------------------------------------------

using Pkg
Pkg.activate(".")
Pkg.instantiate()

# Run the simulation
include("src/Environment_and_Movement.jl")
