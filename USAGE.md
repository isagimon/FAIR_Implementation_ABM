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

---

## Movement options

By default, the movement sampler includes the 18 FCC neighbor directions plus a `"None"` (stay-in-place) option.
This yields a 1/19 probability of no movement per movement attempt. If you prefer strictly moving steps, remove
`"None"` from `Possible_Movement_Options` in `src/Environment_and_Movement.jl`.

---

## Output files

Each simulation run writes its outputs into a run directory:

- `Simulation_Information.csv` — key run parameters (for provenance; single-row summary)
- `Input_Parameters_used.csv` — exact copy of the parameter CSV used for the run (full provenance)
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

This script automatically activates the **analysis environment** at `Analysis/Project.toml` and instantiates
its dependencies (notably `Plots`) before running the plotting steps.

Analysis parameters are read from `Input_Parameters_Analysis.csv` (notably `Directory` and `Total_Timesteps`).

---

## Notes on threading

The model is currently validated for **single-thread** execution. If you encounter instability when running with multiple threads, set:

```bash
export JULIA_NUM_THREADS=1
```

Then run:

```bash
julia run_simulation.jl
```

---

## Generating publication-quality figures (Figures 3, 4, and 5)

The manuscript figures were produced by Python scripts in `Analysis/`, not by
the Julia analysis pipeline. They require a separate Python environment.


### Install Python dependencies

```bash
pip install -r Analysis/requirements.txt
```

Exact versions are pinned in `Analysis/requirements.txt` to match the
environment used to generate the submitted figures.

### Run the scripts

From the repository root:

```bash
python Analysis/Plot_Fig_3_publication.py
python Analysis/Plot_Fig_4_publication.py
python Analysis/Plot_Fig_5_publication.py
```

Each script reads the CSV files from `data/Figure_3/`,
`data/Figure_4/`, or `data/Figure_5/` relative to the repository root 
and writes  three output files alongside the script:

| File                                             | Format | Use |
|--------------------------------------------------|---|---|
| `Figure_3.eps` / `Figure_4.eps` / `Figure_5.eps` | EPS vector | Primary journal submission format |
| `Figure_3.tif` / `Figure_4.tif` / `Figure_5.tif` | TIFF, 600 dpi, RGB | Combination artwork (journal minimum: 600 dpi) |
| `Figure_3.pdf` / `Figure_4.pdf` / `Figure_5.pdf` | PDF vector | Backup; handles transparency correctly |

### What the shaded bands represent

Figures 3 and 4: The shaded band around each line is mean ± SEM (standard 
error of the mean) across the n = 5 independent simulation runs at each 
timestep. SEM communicates the uncertainty of the estimated mean trajectory; 
the individual run values can be inspected directly in the source CSV files.

Figure 5: The dashed vertical line marks the mean runtime and the shaded 
band spans ± 1 SD. SD is used here because the histogram shows the full 
distribution of individual simulation runtimes, and ± 1 SD describes the 
spread of that distribution rather than the uncertainty of the mean.

### Note on the Julia plotting scripts

`Analysis/Plot_Fig_3.jl`, and `Analysis/Plot_Fig_4.jl`
use the Julia `Plots` package and produce equivalent line plots suitable for 
exploratory analysis. They produce PNG output with Plots.jl default styling 
and do not include SEM shading. They are retained for users who prefer to work 
entirely within the Julia environment.

---

## Troubleshooting

- If outputs are not created, confirm `Directory` in `Input_Parameters.csv` is writable.
- If running on Windows, ensure output directories do not contain `:` characters. (This repository uses a Windows-safe timestamp format by default.)

---

