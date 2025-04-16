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

## 🚀 Quick Start

To run the simulation:

```bash
git clone https://github.com/isagimon/FAIR_Implementation_ABM.git
cd FAIR_Implementation_ABM
julia Main_Simulation.jl
```

For full setup instructions, input parameters, and output file descriptions, refer to the [USAGE.md](USAGE.md) file.

---

## 📜 License

This project is licensed under the [Apache License 2.0](LICENSE).  
You may not use this file except in compliance with the License.  
See the `LICENSE` file or visit [http://www.apache.org/licenses/LICENSE-2.0](http://www.apache.org/licenses/LICENSE-2.0) for full terms.


---

## 👩‍🔬 Authors

- Isabella Gimon  
- Santiago Schnell  
- Conner Sandefur  

---

## 📣 Citation

If you use this code in your research, please cite:

> [Authors], "Protein Aggregation Simulation Model", [Year], [Repository URL]

---

## 📬 Contact

For questions or suggestions, feel free to reach out to the authors or open an issue on the GitHub repository.

