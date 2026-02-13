# 🧪 Using the FAIR_Implementation_ABM simulation model

This document provides instructions to set up, configure, run, and post-process the protein aggregation simulation.

---

## Prerequisites

- Julia (≥ 1.8; ≥ 1.10 recommended)
- `git` (optional, for cloning)

---

## Install dependencies

From the repository root:

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```

---

## Run a simulation

From the repository root:

```bash
julia run_simulation.jl
```

This will:

- activate and instantiate the pinned Julia environment (`Project.toml` + `Manifest.toml`)
- load parameters from `Input_Parameters.csv`
- run the ABM simulation
- write outputs to a timestamped run directory inside `Data_Collection/` (by default)

Example output directory:

- `Data_Collection/Simulation_2025-04-16_11-08-52/`

---

## Configure inputs (`Input_Parameters.csv`)

Edit `Input_Parameters.csv` to change:

- `Lattice_Size` — lattice width (simulation volume scales as `Lattice_Size^3`)
- `Max_NumberMonomers_Native`, `Max_NumberMonomers_AggregateProne`
- kinetic probabilities:
  - `Native_to_AggregateProne`
  - `AggregateProne_to_Native`
  - `Oligomer_Formation`
  - `Oligomer_Dissociation_rate`
  - `Fibril_Formation`
  - `Fibril_Growth`
- `Probability_of_Oligomer_Removal` — set to `0` to disable clearance
- `MAX_NumberMovements` — number of timesteps
- crowding controls:
  - `Spheres?` (TRUE/FALSE)
  - `Obstacle_Radius`
  - `Crowder_Concentration_Spheres`
- `Directory` — output root directory (recommended to keep this **relative**, default: `Data_Collection`)

### Portable overrides (optional)

You may override inputs without editing files by setting environment variables **before** starting Julia:

- `FAIR_ABM_PARAMETER_FILE` — alternate parameter CSV path
- `FAIR_ABM_OUTPUT_DIR` — alternate output directory root

---

## Output files

Each simulation run writes its outputs into a run directory:

- `Simulation_Information.csv` — key run parameters (for provenance)
- `Input_Parameters_used.csv` — exact parameter CSV copied into the run directory (for provenance)
- `Oligomer_and_Aggregate_Count_Results.csv`
- `Native_and_AggregateProne_Count_Results.csv`
- `Oligomers_Cleared.csv`
- `MSD_Data.csv`
- `Fibril_Length_Count_Results.csv`
- Optional per-timestep snapshots: `TimestepXX.csv` (if enabled in the code path)

See `DATA_DICTIONARY.md` for column definitions.

---

## Post-simulation analysis

The analysis pipeline aggregates results across multiple run directories and writes summary CSVs/figures under:

- `Data_Collection/Compare_Simulations/`

Run the full pipeline:

```bash
julia Run_All_Analysis_Scripts.jl
```

Analysis parameters are read from `Input_Parameters_Analysis.csv` (notably `Directory` and `Total_Timesteps`).

---

## Troubleshooting

- If outputs are not created, confirm `Directory` in `Input_Parameters.csv` is writable.
- If running on Windows, ensure output directories do not contain `:` characters. (This repository uses a Windows-safe timestamp format by default.)
