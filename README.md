# 🧬 Protein Aggregation Simulation (FAIR_Implementation_ABM)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15230587.svg)](https://doi.org/10.5281/zenodo.15230587)

This repository provides a Julia implementation of an agent-based model (ABM) for protein aggregation on a 3D FCC lattice. The model captures stochastic monomer movement, conformational switching (Native ↔ AggregateProne), oligomer formation/dissociation, fibril/aggregate growth, macromolecular crowding (optional static spherical obstacles), and an optional oligomer clearance mechanism implemented as stochastic removal.

---

## Quick start

```bash
git clone https://github.com/isagimon/FAIR_Implementation_ABM.git
cd FAIR_Implementation_ABM
julia run_simulation.jl
```

Outputs are written to a timestamped run directory under `Data_Collection/` by default (e.g., `Data_Collection/Simulation_2025-04-16_11-08-52/`).

---

## Repository structure

- `src/FAIR_Implementation_ABM.jl` — Julia package entry point (public API)
- `src/Agents.jl` — lattice generation and initial agent assignment
- `src/Environment_and_Movement.jl` — movement/aggregation dynamics, data collection, and export
- `run_simulation.jl` — convenience script to activate/instantiate the Julia environment and run one simulation
- `Analysis/` — post-simulation analysis scripts
- `Input_Parameters.csv` — simulation configuration
- `Input_Parameters_Analysis.csv` — analysis configuration
- `USAGE.md` — detailed usage and output descriptions
- `DATA_DICTIONARY.md` — column definitions for output files

---

## Post-simulation analysis (`Analysis/`)

The analysis scripts operate on simulation run folders under `Data_Collection/` and write aggregated results to `Data_Collection/Compare_Simulations/`.

Run the full analysis pipeline:

```bash
julia Run_All_Analysis_Scripts.jl
```

Key scripts:
- `Run_All_Analysis_Scripts.jl` — master pipeline driver
- `Append_AggregateProne_and_Native.jl` — aggregates Native/AggregateProne counts across runs
- `Append_Aggregate_and_Oligomer.jl` — aggregates Aggregate/Oligomer counts across runs
- `Average_All_Monomers_vs_Timesteps.jl` — plots average counts over time
- `Append_Oligomer_Clearance_Data.jl` — aggregates cleared-oligomer counts across runs
- `Average_Oligomers_Cleared_vs_Timesteps.jl` — plots average cleared counts over time

---

## Configuration

Edit `Input_Parameters.csv` to change:
- lattice size and agent counts
- kinetic probabilities
- number of timesteps
- whether crowding is enabled and its parameters
- output location (`Directory`)

For portability, `Directory` should be a relative path (default: `Data_Collection`).

Optional environment overrides (useful in HPC/CI):
- `FAIR_ABM_PARAMETER_FILE` — path to an alternate parameter CSV (instead of `Input_Parameters.csv`)
- `FAIR_ABM_OUTPUT_DIR` — override the output root directory (instead of `Directory`)

---

## Dependencies

The core model uses standard Julia packages including:
- `Random`, `Dates`
- `CSV`, `DataFrames`
- `Plots`

(See `Project.toml` / `Manifest.toml` for the fully pinned environment.)

---

## License

Apache License 2.0. See `LICENSE`.

---

## Authors

- Santiago Schnell
- Conner Sandefur
- Isabella Gimón

---

## Citation

Please cite the repository (see `CITATION.cff`) if you use it in published work.
