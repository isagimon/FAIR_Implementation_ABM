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
- `Directory` — output root directory (optional)

### Portable overrides (optional)

You may override inputs without editing files by setting environment variables **before** starting Julia:

- `FAIR_ABM_PARAMETER_FILE` — alternate parameter CSV path
- `FAIR_ABM_OUTPUT_DIR` — alternate output directory root

## Movement options

The model samples from 18 FCC nearest-neighbor movement directions plus an explicit `"None"` option, which represents remaining in place for a timestep.

This means a monomer may:
- move to one of 18 neighboring FCC lattice positions, or
- remain in place when `"None"` is sampled

The `"None"` option has a probability of 1/19 under uniform sampling.

---

## Output files

Each simulation run writes its outputs into a run directory:

- `Simulation_Information.csv` — key run parameters (for provenance)
- `Input_Parameters_used.csv` — exact copy of the parameter file used for the run
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


Analysis parameters are read from `Input_Parameters_Analysis.csv` (notably `Directory` and `Total_Timesteps`).

---

## Generating publication-quality figures (Figures 3, 4, and 5)

Publication-quality figure generation is supported through the Python scripts in `Analysis/`.

Available scripts:
- `Analysis/Plot_Fig_3_publication.py`
- `Analysis/Plot_Fig_4_publication.py`
- `Analysis/Plot_Fig_5_publication.py`

These scripts use Python dependencies pinned in:

- `Analysis/requirements.txt`

Install them with:

```bash
pip install -r Analysis/requirements.txt