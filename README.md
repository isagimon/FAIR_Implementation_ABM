# Protein Aggregation Simulation

## Description

This Julia project simulates the dynamics of protein aggregation on a 3D lattice. It models the movement and interactions of monomers, including state transitions (Native, Amyloid, Oligomer, Fibril) and aggregation processes.

## Files

* `Agents.jl`: Contains functions for generating the simulation lattice and assigning initial states to monomers.
* `Main_Simulation.jl`: Contains the main simulation logic, including monomer movement, state changes, aggregation rules, and data collection.

## Dependencies

* Julia (version X.X or later)
* Packages:
    * `Random`
    * `Plots`
    * `DataFrames`
    * `CSV`
    * `Dates`
    * `XLSX`
    * `Profile`
    * `Base.Threads`

## Installation

1.  Ensure you have Julia installed. You can download it from [https://julialang.org/downloads/](https://julialang.org/downloads/).
2.  Clone this repository to your local machine:
    ```bash
    git clone [repository_url]
    cd [repository_directory]
    ```
3.  Install the required Julia packages. Open the Julia REPL and run:
    ```julia
    using Pkg
    Pkg.activate(".")  # Activate the project environment (if any)
    Pkg.instantiate() # Install all dependencies in the Project.toml file
    ```

## Usage

To run the simulation:

1.  Navigate to the project directory in your Julia REPL or terminal.
2.  Execute the `Main_Simulation.jl` script:
    ```julia
    julia Main_Simulation.jl
    ```

The simulation will generate output files in the directory specified in `Main_Simulation.jl`.

## Input Parameters

The simulation parameters are defined at the beginning of `Main_Simulation.jl`:

* `MAX_NumberMovements`: Maximum number of simulation timesteps.
* `Native_to_Amyloid`: Probability of native to amyloid transition.
* `Amyloid_to_Native`: Probability of amyloid to native transition.
* `Oligomer_Formation`: Probability of oligomer formation.
* `Oligomer_Dissociation_rate`: Probability of oligomer dissociation.
* `Fibril_Formation`: Probability of fibril formation.
* `Fibril_Growth`: Probability of fibril growth.
* `Fibril_No_Growth`: Probability of fibril no growth.
* `Lattice_Size`: Size of the cubic lattice.
* `Max_NumberMonomers_Native`: Maximum number of native monomers.
* `Max_NumberMonomers_Amyloid`: Maximum number of amyloid-prone monomers.
* `Obstacle_Radius`: Radius of spherical crowders (if applicable).
* `Crowder_Concentration_Spheres`: Concentration of spherical crowders (if applicable).
* `Obstacle`: Boolean to enable/disable spherical crowders.
* `Sphere_Volume`: Volume of a single sphere (in lattice units).

## Output

The simulation generates several output files (the exact names and locations are defined in `Main_Simulation.jl`):

* `Simulation_Results.csv`: Contains the counts of oligomers and aggregates over time.
* `Native_and_Amyloid_Count_Results.csv`: Contains the counts of native and amyloid monomers over time.
* `MSD_Data.csv`: Contains the Mean Squared Displacement (MSD) data.
* `Fibril_Length_Count_Results.csv`: Contains the counts of fibrils of different lengths over time.
* `TimestepXXX.csv`: Files containing the state of the lattice at specific timesteps.

A more detailed description of the output data columns can be found in a separate `DATA_DICTIONARY.md` file (or similar).

## Data Dictionary

* (Create a separate file or section to describe the columns of your output CSV files)
    * Example:
        * `Simulation_Results.csv`
            * `Timestep`: The simulation timestep (integer).
            * `Oligomers`: The number of oligomers at this timestep (integer).
            * `Aggregates`: The number of aggregates at this timestep (integer).
        * ...

## Assumptions and Limitations

* The simulation uses a simplified model of protein aggregation.
* The movement of monomers is stochastic and does not consider complex interactions.
* The lattice has periodic boundary conditions.

## Contributing

(Optional: If you want others to contribute)

1.  Fork the repository.
2.  Create a new branch for your feature or bug fix.
3.  Make your changes and commit them.
4.  Push to the branch.
5.  Create a pull request.

## License

(Specify the license under which your code is released, e.g., MIT License)

This project is licensed under the [MIT License](LICENSE). See the `LICENSE` file for details.

## Authors

* [Your Name(s)]

## Acknowledgements

(Optional: If you want to acknowledge anyone)

* (Any funding sources, collaborators, etc.)

## Citation

If you use this code in your research, please cite it as:

* [Authors], "[Project Title]", [Year], [Repository URL]

## Contact

* [Your Email Address] (Optional)
