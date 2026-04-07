# 📊 Data Dictionary for Protein Aggregation Simulation Output

This document describes the structure and contents of all `.csv` files generated during a simulation run of the FAIR_Implementation_ABM model.

Each simulation run creates a new output folder inside the `Data_Collection/` directory, with the following naming convention:

```
Simulation_<timestamp>/
```

Example:
```
Data_Collection/Simulation_2025-04-16_11-08-52_AM/
```

Each subfolder contains a complete set of outputs described below.

## 1. `Oligomer_and_Aggregate_Count_Results.csv`

This file contains the counts of oligomers and aggregates at each timestep of the simulation.

* **File Description:** `Oligomer_and_Aggregate_Count_Results.csv`
* **Data Type:** Comma-Separated Values (CSV)
* **Columns:**

    * `Timestep`
        * **Description:** The simulation timestep. This represents a discrete point in time within the simulation.
        * **Data Type:** Integer
        * **Units:** Timesteps (dimensionless)
        * **Example:** 1, 100, 5000

    * `Oligomers`
        * **Description:** The number of oligomers present in the simulation lattice at the given timestep. Oligomers are defined as aggregates of monomers in state 3.
        * **Data Type:** Integer
        * **Units:** Number of oligomers
        * **Example:** 0, 25, 120

    * `Aggregates`
        * **Description:** The number of larger aggregates (fibrils) present in the simulation lattice at the given timestep. Aggregates are defined as aggregates of monomers in state 4.
        * **Data Type:** Integer
        * **Units:** Number of aggregates
        * **Example:** 0, 5, 40

## 2. `Native_and_AggregateProne_Count_Results.csv`

This file contains the counts of native and AggregateProne monomers at each timestep.

* **File Description:** `Native_and_AggregateProne_Count_Results.csv`
* **Data Type:** Comma-Separated Values (CSV)
* **Columns:**

    * `Timestep`
        * **Description:** The simulation timestep.
        * **Data Type:** Integer
        * **Units:** Timesteps (dimensionless)
        * **Example:** 1, 100, 5000

    * `Native`
        * **Description:** The number of monomers in the native state (state 1) at the given timestep.
        * **Data Type:** Integer
        * **Units:** Number of monomers
        * **Example:** 450, 300, 100

    * `AggregateProne`
        * **Description:** The number of monomers in the AggregateProne state (state 2) at the given timestep.
        * **Data Type:** Integer
        * **Units:** Number of monomers
        * **Example:** 50, 200, 400

## 3. `MSD_Data.csv`

This file contains the Mean Squared Displacement (MSD) data for monomers and aggregates.

* **File Description:** `MSD_Data.csv`
* **Data Type:** Comma-Separated Values (CSV)
* **Columns:**

    * `Timestep`
        * **Description:** The simulation timestep.
        * **Data Type:** Integer
        * **Units:** Timesteps (dimensionless)
        * **Example:** 1, 100, 5000

    * `MSD_Monomer`
        * **Description:** The Mean Squared Displacement of individual monomers at the given timestep. This measures the average distance monomers have moved over time.
        * **Data Type:** Float64
        * **Units:** Lattice units squared
        * **Example:** 1.25, 5.78, 22.10

    * `MSD_Aggregate`
        * **Description:** The Mean Squared Displacement of aggregates (oligomers and fibrils) at the given timestep. This measures the average distance aggregates have moved over time.
        * **Data Type:** Float64
        * **Units:** Lattice units squared
        * **Example:** 0.50, 2.34, 10.55

## 4. `Fibril_Length_Count_Results.csv`

This file contains the counts of fibrils of different lengths over time.

* **File Description:** `Fibril_Length_Count_Results.csv`
* **Data Type:** Comma-Separated Values (CSV)
* **Columns:**

    * `Timesteps`
        * **Description:** The simulation timestep.
        * **Data Type:** Integer
        * **Units:** Timesteps (dimensionless)
        * **Example:** 1, 100, 5000

    * `1`, `2`, `3`, ..., `N` (where `N` is `max_fibril_size`)
        * **Description:** Columns representing the count of fibrils of each length. For example, column `5` contains the number of fibrils with a length of 5 monomers.
        * **Data Type:** Integer
        * **Units:** Number of fibrils
        * **Example:** 0, 10, 5 (for a given length)

## 5. `Oligomers_Cleared.csv`

This file contains the cumulative number of monomers cleared from oligomers at each simulation timestep.

* **File Description:** `Oligomers_Cleared.csv`  
* **Data Type:** Comma-Separated Values (CSV)  
* **Columns:**

  * `Timestep`  
      * **Description:** The simulation timestep.  
      * **Data Type:** Integer  
      * **Units:** Timesteps (dimensionless)  
      * **Example:** 0, 1, 100, 5000  

  * `Number_Monomers_Cleared`  
      * **Description:** The cumulative number of monomers cleared from oligomers at the given timestep.  
      * **Data Type:** Integer  
      * **Units:** Number of monomers  
      * **Example:** 0, 2, 15, 100  

## 6. `Simulation_Information.csv`

This file records the *key simulation settings* used for a specific run as a **single-row provenance snapshot**.
It is written automatically into each run folder (e.g., `Data_Collection/Simulation_<run_id>/`).

* **File Description:** `Simulation_Information.csv`
* **Data Type:** Comma-Separated Values (CSV)
* **Format:** One header row (parameter names) + one data row (values).

**Columns (header fields):**

- `Run_ID` — Run identifier used in the output directory name
- `Lattice_Size` — FCC lattice width (simulation volume scales as `Lattice_Size^3`)
- `Number_Native_Monomers`
- `Number_AggregateProne_Monomers`
- `Timesteps`
- `P_Native_to_AggregateProne`
- `P_AggregateProne_to_Native`
- `P_Oligomer_Formation`
- `P_Oligomer_Dissociation`
- `P_Fibril_Formation`
- `P_Fibril_Growth`
- `P_Oligomer_Removal`
- `Crowders_Enabled`
- `Obstacle_Radius`
- `Crowder_Concentration_Spheres`

**Example:**

```csv
Run_ID,Lattice_Size,Number_Native_Monomers,Number_AggregateProne_Monomers,Timesteps,...
2026-02-13_13-46-05,3,15,15,1000,...
```

---

### Publication Figures — Statistical Convention

Figures 3 and 4 in the manuscript display **mean ± SEM** trajectories
computed from the ensemble CSV files described above.

For each timestep *t* and each species, the plotted values are:

- **Line**: arithmetic mean across the *n* simulation columns at timestep *t*
- **Shaded band**: mean ± one standard error of the mean
  (SEM = SD / √n, where n = 5 simulations)

SEM is used rather than SD because the figures show the precision of the
estimated mean trajectory across runs, not the spread of individual
simulation outcomes. The full per-run data for every timestep are available
in the source CSV files listed above, enabling readers to compute any
alternative summary statistic.

The raw ensemble CSV files that these scripts read are stored in
`data/Figure_3/` and `data/Figure_4/` in the repository root.
The scripts that produce the figures are
`Analysis/Plot_Fig_3_publication.py` and
`Analysis/Plot_Fig_4_publication.py`. See `USAGE.md`
for instructions.


## 7. `Input_Parameters_used.csv`

This is an **exact copy** of the parameter CSV used to launch the run (either `Input_Parameters.csv`
or an alternate file passed via `FAIR_ABM_PARAMETER_FILE`). It is copied into the run directory to
support full provenance and reproducibility.

* **File Description:** `Input_Parameters_used.csv`
* **Data Type:** Comma-Separated Values (CSV)
* **Columns:**

  * `Input_Parameters` — Parameter name
  * `Values` — Value used for the run
  * `Instructions` — Human-readable guidance / notes

**Notes:**
- This file preserves the full original parameter table (including instructions), which is not
  contained in `Simulation_Information.csv`.
- Together, `Simulation_Information.csv` (machine-friendly summary) and `Input_Parameters_used.csv`
  (full provenance copy) allow reconstruction of the run configuration.


## Analysis Output Files (`Compare_Simulations/`)

The following files are automatically generated by the post-processing analysis scripts in the analysis/ folder. Each CSV file aggregates data across multiple simulation runs for comparison and visualization.

Each column (after the first) corresponds to a unique simulation, and each row corresponds to a simulation timestep.

## 8. `Appending_AggregateProne_Count.csv`

This file contains the number of AggregateProne monomers (state 2) recorded at each timestep across multiple simulations.

* **File Description:** `Appending_AggregateProne_Count.csv`
* **Data Type:** Comma-Separated Values (CSV)
* **Columns:**

    * `Timesteps`
        * **Description:** The simulation timestep. This represents a discrete point in time shared across all simulations.
        * **Data Type:** Integer
        * **Units:** Timesteps (dimensionless)
        * **Example:** 0, 1, 100, 5000

    * `Simulation_<timestamp>`
        * **Description:** Number of AggregateProne monomers in a specific simulation at each timestep. One column per simulation.
        * **Data Type:** Integer
        * **Units:** Number of monomers
        * **Example:** 50, 200, 400


## 9. `Appending_Native_Count.csv`

This file contains the number of native monomers (state 1) recorded at each timestep across multiple simulations.

* **File Description:** `Appending_Native_Count.csv`
* **Data Type:** Comma-Separated Values (CSV)
* **Columns:**

    * `Timesteps`
        * **Description:** The simulation timestep.
        * **Data Type:** Integer
        * **Units:** Timesteps (dimensionless)
        * **Example:** 0, 1, 100, 5000

    * `Simulation_<timestamp>`
        * **Description:** Number of native monomers in a specific simulation at each timestep. One column per simulation.
        * **Data Type:** Integer
        * **Units:** Number of monomers
        * **Example:** 300, 200, 100


## 10. `Appending_Oligomer_Count.csv`

This file contains the number of oligomers (state 3) recorded at each timestep across multiple simulations.

* **File Description:** `Appending_Oligomer_Count.csv`
* **Data Type:** Comma-Separated Values (CSV)
* **Columns:**

    * `Timesteps`
        * **Description:** The simulation timestep.
        * **Data Type:** Integer
        * **Units:** Timesteps (dimensionless)
        * **Example:** 0, 1, 100, 5000

    * `Simulation_<timestamp>`
        * **Description:** Number of oligomers in a specific simulation at each timestep. One column per simulation.
        * **Data Type:** Integer
        * **Units:** Number of oligomers
        * **Example:** 0, 25, 120


## 11. `Appending_Aggregate_Count.csv`

This file contains the number of aggregates (state 4, e.g., fibrils) recorded at each timestep across multiple simulations.

* **File Description:** `Appending_Aggregate_Count.csv`
* **Data Type:** Comma-Separated Values (CSV)
* **Columns:**

    * `Timesteps`
        * **Description:** The simulation timestep.
        * **Data Type:** Integer
        * **Units:** Timesteps (dimensionless)
        * **Example:** 0, 1, 100, 5000

    * `Simulation_<timestamp>`
        * **Description:** Number of fibrils or large aggregates in a specific simulation at each timestep. One column per simulation.
        * **Data Type:** Integer
        * **Units:** Number of aggregates
        * **Example:** 0, 5, 40

## 12. `Appending_Oligomers_Clearance_Count.csv`

This file compiles the cumulative number of monomers cleared from oligomers across multiple simulation runs, aligned by timestep.

* **File Description:** `Appending_Oligomers_Clearance_Count.csv`  
* **Data Type:** Comma-Separated Values (CSV)  
* **Columns:**

  * `Timesteps`  
      * **Description:** The simulation timestep shared across all runs.  
      * **Data Type:** Integer  
      * **Units:** Timesteps (dimensionless)  
      * **Example:** 0, 1, 100, 5000  

  * `Simulation_<timestamp>`  
      * **Description:** The cumulative number of monomers cleared from oligomers at each timestep for a specific simulation run.  
      * **Data Type:** Integer  
      * **Units:** Number of monomers  
      * **Example:** 0, 2, 4, 12  
      * **Notes:**  
        - Each simulation column is named using its folder timestamp (e.g., `Simulation_8:26:2025 Time 12-36-33.387 PM`).  
        - Additional simulations will appear as new columns in this file.  


## 13. `Average_All_Monomers_States_vs_Timesteps_<timestamp>.png`

This image shows the average number of monomers in each state—Native, AggregateProne, Oligomer, and Aggregate—across all simulations over time.

* **File Description:** `Average_All_Monomers_States_vs_Timesteps_<timestamp>.png`
* **Data Type:** Portable Network Graphics (PNG)
* **Location:** `Compare_Simulations/`
* **Plot Contents:**

    * **X-axis:** `Timesteps`
        * **Description:** Simulation timestep index.
        * **Units:** Timesteps (dimensionless)
        * **Example:** 0, 1, 25, 50, 100

    * **Y-axis:** `Average Count`
        * **Description:** Average number of monomers in each state across all simulations at each timestep.
        * **Units:** Count of monomers
        * **Example:** 0.0, 1.5, 7.0, 9.0

    * **Line Colors and Labels:**
        * `Orange` — AggregateProne monomers (State 2)
        * `Green` — Native monomers (State 1)
        * `Red` — Oligomers (State 3)
        * `Blue` — Aggregates (State 4)

* **Notes:**
    * This plot aggregates data from:
        * `Appending_AggregateProne_Count.csv`
        * `Appending_Native_Count.csv`
        * `Appending_Oligomer_Count.csv`
        * `Appending_Aggregate_Count.csv`
    * It is generated using the script: `Average_All_Monomers_vs_Timesteps.jl`


## 14. `Oligomers_Cleared_vs_Timesteps_<timestamp>.png`

This image shows the average number of monomers cleared from oligomers across all simulations over time.

* **File Description:** `Oligomers_Cleared_vs_Timesteps_<timestamp>.png`  
* **Data Type:** Portable Network Graphics (PNG)  
* **Location:** `Compare_Simulations/`  
* **Plot Contents:**

    * **X-axis:** `Timesteps`  
        * **Description:** Simulation timestep index.  
        * **Units:** Timesteps (dimensionless)  
        * **Example:** 0, 1, 100, 5000  

    * **Y-axis:** `Average Number of Monomers Cleared`  
        * **Description:** Average cumulative number of monomers cleared from oligomers across all simulations at each timestep.  
        * **Units:** Number of monomers  
        * **Example:** 0.0, 2.5, 10.0, 50.0  

    * **Line Style and Label:**  
        * `Blue` — Average number of monomers cleared  

* **Notes:**  
    * This plot aggregates data from:  
        * `Appending_Oligomers_Clearance_Count.csv`  
    * It is generated using the script: `Average_Oligomers_Cleared_vs_Timesteps.jl`  

 
