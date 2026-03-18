# Coauthor Change Summary

This document summarizes the publication-readiness updates made to the `FAIR_Implementation_ABM` repository. These changes were made to improve repository portability, reproducibility, documentation consistency, and metadata completeness. No scientific or algorithmic intent of the model was changed.

## Summary of updates

### Repository hygiene and portability
- Removed unnecessary system files and generated documentation artifacts
- Added `.gitignore` rules for generated outputs, build artifacts, and editor files
- Replaced hard-coded absolute paths with portable relative path handling
- Ensured timestamped output directories use filesystem-safe naming

### Package structure and reproducibility
- Added a standard Julia package entry point
- Removed side-effect execution at import time
- Updated `Project.toml` metadata
- Regenerated `Manifest.toml` after dependency cleanup
- Added provenance outputs for each run:
  - `Input_Parameters_used.csv`
  - `Simulation_Information.csv`

### Documentation and metadata consistency
- Updated documentation to reflect current parameter names and valid ranges
- Corrected terminology from `Amyloid` to `AggregateProne` where needed
- Updated `CITATION.cff` metadata, author order, and abstract
- Added `.zenodo.json` for machine-readable Zenodo deposit metadata
- Corrected and expanded the data dictionary

### Testing and continuous integration
- Added a minimal test suite
- Added GitHub Actions CI workflow to run tests and build documentation on pushes and pull requests

## Scientific scope
These updates do **not** change the scientific purpose of the model, the aggregation rules, or the intended interpretation of simulation results. They are repository-quality and reproducibility improvements intended to support sharing, reuse, citation, and publication.