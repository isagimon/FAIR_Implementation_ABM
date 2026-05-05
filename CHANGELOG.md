# Changelog

All notable changes to `FAIR_Implementation_ABM` are documented here.

This project is archived on Zenodo. The concept DOI
[10.5281/zenodo.15230586](https://doi.org/10.5281/zenodo.15230586) always
resolves to the latest version. Individual version DOIs are listed below.

---
## v2.2.3 - Figure 5 Runtime Data Update

**Zenodo DOI:** [10.5281/zenodo.20032424](https://doi.org/10.5281/zenodo.20032424)

### Changed
- Updated Figure 5 runtime benchmark data and associated repository files for consistency.

## [v2.2.2] — 2026-04-02

**Zenodo DOI:** [10.5281/zenodo.15230586](https://doi.org/10.5281/zenodo.15230586)

This release is the final, publication-ready version of the codebase submitted
with the manuscript. It adds the manuscript data and publication-quality figure
scripts, separates the analysis environment from the core simulation
environment, and hardens output-directory resolution and thread-safety in
oligomer clearance. No scientific or algorithmic logic was changed relative to
v2.1.1.

### Added

- **Publication figure scripts (Python).** Three Python scripts in `Analysis/`
  reproduce the submitted manuscript figures at journal-required resolution:
  - `Plot_Fig_3_publication.py` — monomer state counts (Figure 3; line plot
    with mean ± SEM shading across *n* = 5 runs).
  - `Plot_Fig_4_publication.py` — aggregate count with vs. without oligomer
    clearance (Figure 4; same shading convention).
  - `Plot_Fig_5_publication.py` — simulation runtime distribution (Figure 5;
    histogram with mean ± 1 SD band).
  - Each script writes three output formats: EPS (primary journal submission),
    TIFF at 600 dpi (combination artwork), and PDF (backup).
  - `Analysis/requirements.txt` pins the exact Python package versions
    (matplotlib, numpy, pandas, Pillow) used to generate the submitted figures
    (Python 3.12.3 on Ubuntu 24.04).

- **Exploratory Julia figure scripts.** `Analysis/Plot_Fig_3.jl` and
  `Analysis/Plot_Fig_4.jl` produce equivalent PNG line plots using `Plots.jl`
  for users who prefer to work entirely within Julia. They do not include SEM
  shading and are retained for exploratory use.

- **Manuscript data directory.** `data/` contains the ensemble CSV files used
  to generate each figure, organized by figure number:
  - `data/Figure_3/` — four CSVs of monomer-state counts aggregated across
    five independent runs (no crowding; no oligomer clearance).
  - `data/Figure_4/No_Oligomer_Clearance/` and `data/Figure_4/Oligomer_Clearance/`
    — aggregate counts comparing baseline against clearance enabled (*p* = 0.05).
  - `data/Figure_5/Simulation_Runtimes_Minutes.csv` — wall-clock runtimes (in
    minutes) for 300 benchmark runs executed on the Notre Dame Center for
    Research Computing (1,000 monomers; 5,000 timesteps; 30 × 30 × 30 FCC
    lattice; 15 parallel jobs).
  - `data/README.md` documents the contents, provenance, and plotting scripts
    for each subfolder.

- **Separate analysis Julia environment.** `Analysis/Project.toml` defines an
  independent Julia environment for the post-processing pipeline, keeping
  plotting dependencies (e.g., `Plots`) out of the core simulation environment.

- **`Data_Collection/Compare_Simulations/.gitkeep`** — placeholder file to
  track the output directory for aggregated analysis results in version control.

- **Root-level `zenodo.json`** added alongside the existing `.zenodo.json` for
  compatibility with Zenodo's deposit workflow.

- **DOI badge** added to `README.md`.

### Changed

- **Analysis environment activation.** `Run_All_Analysis_Scripts.jl` now
  activates `Analysis/Project.toml` (`Pkg.activate(joinpath(@__DIR__,
  "Analysis"))`) instead of the repository root environment, so that `Plots`
  and other plotting dependencies are isolated from the core simulation
  environment.

- **Output-directory resolution hardened.** The precedence logic in
  `Environment_and_Movement.jl` (both at load time and inside
  `Resolve_Output_Directory`) was rewritten to guard against the edge case
  where a relative path resolves to the repository root itself, ensuring
  outputs always land under `Data_Collection/` rather than at the repo root.
  Relative paths in `FAIR_ABM_OUTPUT_DIR` and `Input_Parameters.csv` are now
  resolved consistently against the repository root rather than the Julia
  working directory.

- **`run_simulation` signature.** The `output_root` keyword argument changed
  from `AbstractString` with a default of `OUTPUT_ROOT` to
  `Union{Nothing,AbstractString}` with a default of `nothing`, so that
  directory resolution is fully delegated to `Resolve_Output_Directory` at
  call time rather than at module load time.

- **`Make_Directory` provenance copy.** The copy of the parameter CSV into the
  run directory (`Input_Parameters_used.csv`) is now guarded by an `isfile()`
  check before attempting the copy, preventing errors in test environments
  where a real parameter file may not be present.

- **Thread-safe oligomer clearance.** The oligomer-clearance routine in
  `Environment_and_Movement.jl` now collects all coordinates belonging to a
  chosen oligomer inside a `lock(dict_lock)` block before removing them,
  eliminating a potential race condition in multi-threaded execution. The
  removed-count tracking was simplified to `length(coords_to_clear)`.

- **Core dependencies slimmed further.** `Plots` and `Profile` were removed
  from the root `Project.toml` (they are no longer needed in the core
  simulation environment). `Manifest.toml` was regenerated and reduced from
  ~1,245 to ~312 lines as a result.

- **`Agents.jl` parameter validation.** `"Directory"` was removed from the
  list of required parameters passed to `validate_parameters!`, since the
  output directory is now resolved independently of the parameter CSV and the
  field is optional.

- **`DATA_DICTIONARY.md`** expanded from 331 to 392 lines to document the
  additional output file `Input_Parameters_used.csv` and the manuscript data
  CSVs under `data/`.

- **`USAGE.md`** expanded to include: a dedicated section on generating
  publication-quality figures (Figures 3, 4, and 5) with the Python scripts; a
  note on movement options (18 FCC neighbors plus optional `"None"`
  stay-in-place); a threading caveat recommending single-thread execution; and
  documentation of `Input_Parameters_used.csv` as a provenance artifact per
  run.

- **`README.md`** updated to reflect: the separate analysis environment; the
  new `data/` directory; the additional provenance file
  (`Input_Parameters_used.csv`); the movement-options note; and split
  dependency documentation (core vs. analysis).

- **Metadata.** `CITATION.cff` updated to version `v2.2.2` and DOI
  `10.5281/zenodo.15230586`; author order corrected to Gimón, Sandefur,
  Schnell throughout `CITATION.cff`, `Project.toml`, `README.md`, and
  `zenodo.json`.

- **`docs/Project.toml`** simplified: the self-referential
  `FAIR_Implementation_ABM` dependency was removed, leaving only `Documenter`.

- **`.gitignore`** — `Manifest.toml` removed from the ignore list (the
  manifest is now committed to pin the reproducible environment); trailing
  newline added.

### Removed

- `COAUTHOR_CHANGE_SUMMARY.md` — internal working document from v2.1.1
  summarizing publication-readiness changes; not relevant to the public
  release.

---

## [v2.1.1] — 2026-03-31

**Zenodo DOI:** [10.5281/zenodo.19353711](https://doi.org/10.5281/zenodo.19353711)

*v2.2.3 is the authoritative public release and supersedes v2.2.2.*

This release covered publication-readiness and repository-quality improvements
made between the initial paper-submission snapshot (v1.2.0) and the public
Zenodo deposit. The core aggregation rules (oligomer formation, fibril
elongation, dissociation) and the lattice geometry were preserved. However,
several algorithmic changes were introduced alongside the infrastructure work;
these are documented explicitly below.

### Added

- **Oligomer clearance mechanism.** A new optional stochastic removal process
  was introduced: at each timestep, with probability
  `Probability_of_Oligomer_Removal`, a randomly chosen oligomer and all its
  constituent monomers are removed from the lattice. This mechanism is
  implemented in `Randomly_Remove_Oligomer()`, `Track_Cleared_Monomers()`,
  `Save_Number_Monomers_Cleared()`, and `Export_Monomers_Cleared_Data()`, and
  is controlled by the new `Probability_of_Oligomer_Removal` parameter in
  `Input_Parameters.csv`. Setting this parameter to `0` disables clearance.
  This mechanism was absent from v1.2.0.

- **Standard Julia package entry point.** `src/FAIR_Implementation_ABM.jl`
  was added so the repository can be activated as a proper Julia package with
  `Pkg.activate` and `using FAIR_Implementation_ABM`.

- **`run_simulation.jl` convenience script.** Activates and instantiates the
  pinned Julia environment, then calls `run_simulation()`, giving users a
  single command to execute a simulation without manually entering the Julia
  REPL.

- **Analysis pipeline.** The `Analysis/` directory was introduced with five
  post-processing scripts that aggregate simulation output across multiple run
  directories and produce summary CSV files and plots:
  - `Append_AggregateProne_and_Native.jl`
  - `Append_Aggregate_and_Oligomer.jl`
  - `Average_All_Monomers_vs_Timesteps.jl`
  - `Append_Oligomer_Clearance_Data.jl`
  - `Average_Oligomers_Cleared_vs_Timesteps.jl`

- **`Run_All_Analysis_Scripts.jl`** — master driver that sequentially invokes
  all five analysis scripts.

- **`Input_Parameters_Analysis.csv`** — configuration file for the analysis
  pipeline (specifies `Directory` and `Total_Timesteps`).

- **Minimal test suite.** `test/runtests.jl` was added, covering: filesystem-
  safe timestamp generation; movement-dispatch consistency (all movement keys
  map to a function); output-directory precedence (`ENV` > function argument >
  CSV parameter > default); and per-run provenance-file creation
  (`Simulation_Information.csv` and `Input_Parameters_used.csv`).

- **GitHub Actions CI workflow.** `.github/workflows/ci.yml` runs tests and
  builds documentation on every push and pull request.

- **Per-run provenance output.** Each simulation run now writes
  `Simulation_Information.csv` (a single-row summary of key input parameters)
  and `Input_Parameters_used.csv` (an exact copy of the parameter CSV) into
  the run directory.

- **`CITATION.cff`** — machine-readable citation metadata (CFF format 1.2.0)
  including title, version, DOI, date-released, abstract, author list, license,
  and repository URL.

- **`.zenodo.json`** — machine-readable Zenodo deposit metadata including
  title, version, creators, description, license, upload type, and keywords.

- **`.gitignore`** — rules covering Julia build artifacts, generated output
  directories, macOS `.DS_Store` files, editor files, and log files.

- **`CHANGELOG.md`** — this file (introduced at v2.1.1).

- **`COAUTHOR_CHANGE_SUMMARY.md`** — internal document summarizing
  publication-readiness updates for co-author review.

- **`example_outputs/`** — two demonstration run directories with pre-computed
  CSV outputs and PNG figures, allowing users to verify their installation
  without running a full simulation:
  - `example_outputs/Simulation_Demo_Run/` — simulation CSVs.
  - `example_outputs/Analysis_Demo_Run/` — analysis CSVs and plots.

### Changed

- **Movement sampling corrected from 13 to 18+1 directions.** Both v1.2.0 and
  v2.1.1 define 19 entries in `Possible_Movement_Options` ("One" through
  "Eighteen" plus "None"). However, v1.2.0 sampled with
  `Possible_Movement_Options[rand(1:13)]`, which silently restricted movement
  to only the first 13 directions; directions "Fourteen" through "Eighteen" and
  "None" were dead code. v2.1.1 corrected this to `rand(Possible_Movement_Options)`,
  sampling uniformly from all 18 FCC neighbor directions plus the explicit
  "None" stay-in-place option (giving a 1/19 probability of no movement per
  timestep). This changes the effective diffusion geometry of the model.

- **Sphere volume calculation replaced.** v1.2.0 read a `Sphere_Volume`
  parameter directly from `Input_Parameters.csv` and divided the target
  crowder concentration by it to compute the number of spheres to place.
  v2.1.1 removed `Sphere_Volume` from the parameter CSV and introduced
  `Return_Sphere_Volume()`, a lookup function that maps `Obstacle_Radius`
  (integers 0–6) to pre-computed lattice site counts (1, 20, 130, 400, 920,
  1850, 3300 sites respectively). Users can no longer specify sphere volume
  independently of radius.

- **Thread-safety restructuring in sphere and aggregate movement.** The
  `dictionary_emptied` flag in the sphere-coordinate and aggregate-movement
  functions was changed from a plain `Bool` (subject to data races) to a
  `Threads.Atomic{Bool}`. `Empty_Possible_Coordinates_Movement_Dict()` was
  moved from inside the threaded loop body to after thread completion.
  `Append_Possible_Coordinate_Movements` was given its own `lock(dict_lock)`
  guard. These changes affect simulation behavior under multi-threaded
  execution.

- **Source files moved to `src/`.** `Agents.jl` and `Environment_and_Movement.jl`
  were relocated from the repository root to `src/`, following standard Julia
  package layout.

- **Terminology standardized.** The state name `Amyloid` / `AmyloidProne` was
  renamed to `AggregateProne` throughout source code, documentation, parameter
  files, and data dictionaries, to match the terminology used in the
  manuscript.

- **Hard-coded absolute paths eliminated.** All occurrences of absolute paths
  (e.g., `/Users/…`) were replaced with portable relative path construction
  using `@__DIR__` and `joinpath`, making the codebase runnable on any machine
  without modification.

- **Side-effect execution removed.** In v1.2.0, `include`-ing the source files
  triggered lattice initialization and simulation execution immediately.
  v2.1.1 separates initialization (`initialize_simulation!()`) from execution
  (`run_simulation()`); neither runs automatically at import time.

- **Filesystem-safe timestamp format.** Run-directory names use a
  Windows-compatible timestamp (colons replaced with hyphens), e.g.,
  `Simulation_2025-04-16_11-08-52`.

- **Environment variable overrides.** Two optional environment variables were
  introduced for HPC and CI use:
  - `FAIR_ABM_PARAMETER_FILE` — alternate path to the parameter CSV.
  - `FAIR_ABM_OUTPUT_DIR` — override for the output root directory.

- **`Project.toml` updated.** Version bumped from `0.1.0` to `2.1.1`; UUID
  updated; all three authors added; `[compat]` bounds added for Julia ≥ 1.10,
  CSV 0.10, DataFrames 1, and Plots 1; `[extras]` and `[targets]` sections
  added for the test suite.

- **Dependencies rationalized.** `Documenter`, `ThreadsX`, `XLSX`, and
  `Profile` were removed from the core `Project.toml`. `Manifest.toml` was
  regenerated after the dependency cleanup.

- **`DATA_DICTIONARY.md` substantially expanded** — from 101 to 331 lines —
  to document all output columns, corrected terminology (`AggregateProne`),
  and the new provenance files.

- **`README.md` substantially rewritten** to describe the new package
  structure, provide an updated quick-start example using `run_simulation.jl`,
  document the analysis pipeline and configuration options, and list all
  dependencies with version context.

- **`USAGE.md` updated** to reflect renamed source files, the new
  `run_simulation.jl` entry point, `Input_Parameters.csv` configuration
  options (including environment variable overrides), and the full list of
  output files.

- **`docs/src/index.md`** and `docs/make.jl` updated to reflect the new
  package structure; `docs/Project.toml` added to manage the documentation
  build environment.

- **License.** The placeholder MIT license text in v1.2.0 was replaced with
  the full Apache License 2.0 text.

### Removed

- **`Fibril_No_Growth` parameter** — present in v1.2.0's `Input_Parameters`
  documentation and recorded in `Simulation_Information.csv`, but never used
  in any conditional logic in the simulation. Its removal has no effect on
  dynamics.

- **`Sphere_Volume` parameter** — removed from `Input_Parameters.csv`;
  replaced by the `Return_Sphere_Volume()` lookup function (see Changed above).

- **`docs/build/`** — generated documentation build artifacts removed from
  version control; the build directory is now listed in `.gitignore`.

- **macOS `.DS_Store` files** — removed from version control and added to
  `.gitignore`.

---

## [v1.2.0] — 2025-04-16

**Zenodo DOI:** [10.5281/zenodo.15230587](https://doi.org/10.5281/zenodo.15230587)

Initial snapshot of the codebase archived at the time of paper submission.
