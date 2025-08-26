# 🧪 Using the FAIR_Implementation_ABM Simulation Model

This document provides detailed instructions on how to set up, configure, and run the protein aggregation simulation using the Agent-Based Model implemented in Julia.

---

## 📦 Prerequisites

Before running the simulation, ensure you have the following:

- [Julia](https://julialang.org/downloads/) (version ≥ 1.8 recommended)
- `git` (if you're cloning the repository)
- A text editor (e.g., VS Code, Sublime, Atom)

To install the required Julia packages, open Julia and run:

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```

---

## 🚀 Running the Simulation

After cloning the repository:

```bash
git clone https://github.com/isagimon/FAIR_Implementation_ABM.git
cd FAIR_Implementation_ABM
julia run_simulation.jl
```

This will automatically:
- Activate the project environment
- Install any missing dependencies
- Launch the simulation using `Environment_and_Movement.jl`
- Load all parameters from `Input_Parameters.csv`

Simulation results will be saved to a timestamped folder inside `Data_Collection/`.


---

## 🛠️ Input Parameters

The parameters below can be configured directly at the top of the `Main_Simulation.jl` file:

- `MAX_NumberMovements`: Maximum number of simulation timesteps.
- `Native_to_AggregateProne`: Probability of native to AggregateProne transition.
- `AggregateProne_to_Native`: Probability of AggregateProne to native transition.
- `Oligomer_Formation`: Probability of oligomer formation.
- `Oligomer_Dissociation_rate`: Probability of oligomer dissociation.
- `Fibril_Formation`: Probability of fibril formation.
- `Fibril_Growth`: Probability of fibril growth.
- `Lattice_Size`: Size of the cubic lattice.
- `Max_NumberMonomers_Native`: Initial number of native monomers.
- `Max_NumberMonomers_AggregateProne`: Initial number of AggregateProne-prone monomers.
- `Obstacle_Radius`: Radius of spherical crowders (if used).
- `Crowder_Concentration_Spheres`: Concentration of crowders (optional).
- `Obstacle`: Boolean to enable or disable crowders.
- `Sphere_Volume`: Volume of a single crowder sphere.
- `Probability_of_Oligomer_Removal`: Probability of removing a random oligomer in a timestep

---

## 📤 Output Files

> All simulation results are saved inside a timestamped folder within the `Data_Collection/` directory. Each run creates a new subdirectory named in the format `Simulation_<timestamp>` (e.g., `Data_Collection/Simulation_2025-04-16_11-08-52_AM/`), containing all output files for that specific run. This ensures results from each simulation are cleanly organized and reproducible.

The simulation generates the following outputs:

- `Oligomer_and_Aggregate_Count_Results.csv`: Tracks oligomer and aggregate counts over time.
- `Native_and_AggregateProne_Count_Results.csv`: Counts of native vs. AggregateProne monomers.
- `Oligomers_Cleared.csv`: Number of monomers cleared from oligomeric species
- `MSD_Data.csv`: Mean squared displacement data for monomer diffusion.
- `Fibril_Length_Count_Results.csv`: Tracks fibril length distributions.
- `TimestepXXX.csv`: Snapshots of lattice states at specified timesteps.

---

## 🧾 Data Dictionary

Here’s what the columns represent in key `.csv` files:

### `Oligomer_and_Aggregate_Count_Results.csv`
| Column     | Description                           |
|------------|---------------------------------------|
| Timestep   | Simulation timepoint (integer)        |
| Oligomers  | Count of state 3 monomers (oligomers) |
| Aggregates | Count of state 4 monomers (aggregates/fibrils) |

### `Native_and_AggregateProne_Count_Results.csv`
| Column     | Description                    |
|------------|--------------------------------|
| Timestep   | Simulation timepoint           |
| Native     | Number of native monomers (N)  |
| AggregateProne  | Number of AggregateProne monomers (A) |

### `MSD_Data.csv`
| Column       | Description                                      |
|--------------|--------------------------------------------------|
| Timestep     | Simulation timepoint                             |
| MSD_Monomer  | Average Mean Squared Displacement of all monomers |
| MSD_Aggregate| Average Mean Squared Displacement of all monomers aggregates    |


### `Fibril_Length_Count_Results.csv`
| Column        | Description                                               |
|---------------|-----------------------------------------------------------|
| Fibril_Length | Number of monomers in each fibril                         |
| Count         | How many fibrils of that length were observed at the end |


### `Oligomers_Cleared.csv`
| Column                 | Description                                              |
|------------------------|----------------------------------------------------------|
| Timestep               | Simulation timepoint                                     |
| Number_Monomers_Cleared| Number of individual monomers cleared at each timestep. Since only oligomers can be cleared, this represents monomers removed from oligomeric species. |


---

## 🧠 Tips

- For faster performance on multi-core machines, consider using Julia's multithreading capabilities with `JULIA_NUM_THREADS`.
- All paths and output configurations can be modified in `Main_Simulation.jl`.

---

For more information on the model logic or code structure, consult the in-line comments in `Agents.jl` and `Main_Simulation.jl`, or reach out via the repository’s Issues page.

