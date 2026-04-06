# 🧬 Protein Aggregation Simulation (FAIR_Implementation_ABM)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15230586.svg)](https://doi.org/10.5281/zenodo.15230586)

This repository provides a Julia implementation of an agent-based model (ABM) for protein aggregation on a 3D face-centered cubic (FCC) lattice. The model captures stochastic monomer movement, conformational switching (Native ↔ AggregateProne), oligomer formation/dissociation, fibril/aggregate growth, macromolecular crowding (optional static spherical obstacles), and an optional oligomer clearance mechanism implemented as stochastic removal.

The repository is organized to support reproducibility and FAIR software practices, including pinned environments, per-run provenance files, analysis workflows, and archived manuscript data.

---

## Quick start

```bash
git clone https://github.com/isagimon/FAIR_Implementation_ABM.git
cd FAIR_Implementation_ABM
julia run_simulation.jl
```

Outputs are written to a timestamped run directory under `Data_Collection/` by default (e.g., `Data_Collection/Simulation_2025-04-16_11-08-52/`).

Each run directory includes provenance artifacts such as:
- `Simulation_Information.csv`
- `Input_Parameters_used.csv`

---

## Repository structure

- `src/FAIR_Implementation_ABM.jl` — Julia package entry point (public API)
- `src/Agents.jl` — lattice generation and initial agent assignment
- `src/Environment_and_Movement.jl` — movement/aggregation dynamics, data collection, and export
- `run_simulation.jl` — convenience script to activate/instantiate the Julia environment and run one simulation
- `Run_All_Analysis_Scripts.jl` — driver script for the full post-simulation analysis workflow
- `Analysis/` — post-simulation analysis scripts and separate analysis environment
- `data/` — manuscript data used to generate publication figures
- `Data_Collection/` — simulation outputs and aggregated analysis results
- `Input_Parameters.csv` — simulation configuration
- `Input_Parameters_Analysis.csv` — analysis configuration
- `USAGE.md` — detailed usage and output descriptions
- `DATA_DICTIONARY.md` — column definitions for output files and manuscript data
- `CITATION.cff` — citation metadata
- `.zenodo.json` — Zenodo metadata

---

## Post-simulation analysis (`Analysis/`)

The analysis scripts operate on simulation run folders under `Data_Collection/` and write aggregated results to `Data_Collection/Compare_Simulations/`.

Run the full analysis pipeline:

```bash
julia Run_All_Analysis_Scripts.jl

```
`Run_All_Analysis_Scripts.jl` activates the separate Julia environment defined in `Analysis/Project.toml`, which keeps plotting and analysis dependencies separate from the core simulation environment.

Key scripts:
- `Run_All_Analysis_Scripts.jl` — master pipeline driver
- `Append_AggregateProne_and_Native.jl` — aggregates Native/AggregateProne counts across runs
- `Append_Aggregate_and_Oligomer.jl` — aggregates Aggregate/Oligomer counts across runs
- `Average_All_Monomers_vs_Timesteps.jl` — plots average counts over time
- `Append_Oligomer_Clearance_Data.jl` — aggregates cleared-oligomer counts across runs
- `Average_Oligomers_Cleared_vs_Timesteps.jl` — plots average cleared counts over time

---

## Manuscript data (`data/`)

The `data/` directory contains the ensemble CSV files used to generate the manuscript figures.

- `data/Figure_3/` — monomer-state count data
- `data/Figure_4/` — aggregate-count comparisons with and without oligomer clearance
- `data/Figure_5/` — simulation runtime data

These files are distinct from the example outputs and are intended to preserve the data underlying the publication figures.

---

## Configuration

Edit `Input_Parameters.csv` to change:
- lattice size and agent counts
- kinetic probabilities
- number of timesteps
- whether crowding is enabled and its parameters
- output location (`Directory`, optional)

Optional environment overrides (useful in HPC/CI):
- `FAIR_ABM_PARAMETER_FILE` — path to an alternate parameter CSV (instead of `Input_Parameters.csv`)
- `FAIR_ABM_OUTPUT_DIR` — override the output root directory

Per-run outputs are organized automatically, and the parameter file used for a run is copied into the run directory as `Input_Parameters_used.csv` for provenance.

---

## Movement options

The model defines 18 FCC nearest-neighbor movement directions plus an explicit `"None"` option representing no movement during a timestep. This means a monomer may move to one of 18 neighboring FCC positions or remain in place when `"None"` is sampled.

---

## Dependencies

### Core simulation environment

The core model uses Julia packages defined in:
- `Project.toml`
- `Manifest.toml`

Core packages include:
- `Random`
- `Dates`
- `CSV`
- `DataFrames`

### Analysis environment

The analysis workflow uses a separate Julia environment defined in:
- `Analysis/Project.toml`

Publication-figure scripts also use Python dependencies pinned in:
- `Analysis/requirements.txt`

This separation keeps plotting and analysis dependencies out of the core simulation environment.

---

## License

Apache License 2.0. See `LICENSE`.

---

## Authors

- Isabella Gimón
- Conner Sandefur
- Santiago Schnell

---

## Citation

Please cite the repository (see `CITATION.cff`) if you use it in published work.
