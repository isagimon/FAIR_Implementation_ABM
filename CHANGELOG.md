# Changelog

## 1.2.0 (publication-ready snapshot)

Repository hygiene and portability improvements (no scientific/algorithmic intent changed):

- Removed macOS `.DS_Store` files and removed `docs/build/` artifacts.
- Added `.gitignore` and placeholder `.gitkeep` files for generated-output directories.
- Eliminated hard-coded absolute paths (e.g., `/Users/...`) and made paths portable using `@__DIR__`/`joinpath`.
- Introduced a standard Julia package entry point (`src/FAIR_Implementation_ABM.jl`).
- Removed side-effect execution at import time:
  - lattice initialization now happens via `initialize_simulation!()` called by `run_simulation()`
  - simulations no longer auto-run when source files are included
- Implemented a filesystem-safe timestamp format for run directories (Windows compatible).
- Added provenance outputs per run:
  - `Input_Parameters_used.csv`
  - `Simulation_Information.csv`
- Added a minimal test suite (`test/runtests.jl`) and a GitHub Actions workflow (`.github/workflows/ci.yml`).
- Updated documentation (`README.md`, `USAGE.md`) and metadata (`.zenodo.json`, `Project.toml`) for internal consistency.

