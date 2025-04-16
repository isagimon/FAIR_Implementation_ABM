# 🧬 Protein Aggregation Simulation

This Julia project simulates the dynamics of protein aggregation on a 3D lattice. It models the stochastic movement and interaction of monomers, including state transitions (Native, Amyloid, Oligomer, Fibril) and the formation of larger aggregates under crowding and kinetic constraints.

---
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15230587.svg)](https://doi.org/10.5281/zenodo.15230587)

## 📂 Files

- `Agents.jl`: Functions for lattice generation and monomer initialization.
- `Main_Simulation.jl`: Core logic for monomer movement, state transitions, aggregation rules, and data collection.
- `USAGE.md`: Complete instructions for running the model and interpreting output.
- `LICENSE`: License file (MIT).

---

## ⚙️ Dependencies

This project uses the following Julia packages:

- `Random`
- `Plots`
- `DataFrames`
- `CSV`
- `Dates`
- `XLSX`
- `Profile`
- `Base.Threads`

---
## ⚠️ Assumptions and Limitations

- The simulation is based on a simplified kinetic model of protein aggregation.
- Monomer movement is stochastic and may not capture all biophysical constraints.
- Periodic boundary conditions are applied to simulate a continuous space.
- Crowders (if enabled) are modeled as static, spherical obstacles with fixed radii.
---

## 📜 License

This project is licensed under the [Apache License 2.0](LICENSE).  
You may not use this file except in compliance with the License.  
See the `LICENSE` file or visit [http://www.apache.org/licenses/LICENSE-2.0](http://www.apache.org/licenses/LICENSE-2.0) for full terms.


---

## 👩‍🔬 Authors

- Santiago Schnell  
- Conner Sandefur  
- Isabella Gimon
---

## 📣 Citation


This repository is licensed under the [Apache License 2.0](LICENSE), which permits reuse, modification, and redistribution under certain conditions.

You are welcome to:
- Use the code for research or educational purposes
- Modify parameters or reaction rules for your own models
- Extend the codebase to include new features

Please cite this repository if you use it in published work:
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15230587.svg)](https://doi.org/10.5281/zenodo.15230587)

For citation details, see the [`CITATION.cff`](CITATION.cff) file.


---

## 📬 Contact
---

If you have questions or would like to collaborate, please open an issue or contact the maintainer:
📧 **Isabella Gimon** (igimon@nd.edu)

