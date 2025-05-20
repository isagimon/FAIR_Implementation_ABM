using CSV
using Statistics
using Plots
using DataFrames
using Dates

"""
    run_plot_all_monomer_states(directory::String, Total_Timesteps::Int)

Generates a combined line plot showing the average count of each monomer state 
(Native, Amyloid, Oligomer, Aggregate) across all simulations over time. This 
function is part of the FAIR post-analysis pipeline and uses pre-aggregated 
CSV data stored in the `Compare_Simulations` folder.

# Arguments
- `directory::String`: Path to the root directory containing the `Compare_Simulations` folder.
- `Total_Timesteps::Int`: Total number of timesteps included in the simulation.

# Behavior
- Reads four input CSVs: 
  - `Appending_Amyloid_Count.csv`
  - `Appending_Native_Count.csv`
  - `Appending_Oligomer_Count.csv`
  - `Appending_Aggregate_Count.csv`
- Each file contains simulation results across multiple runs. The function computes
  the mean value at each timestep by averaging across all simulations.
- All four species are plotted together on the same graph against time (0 to Total_Timesteps).
- The final figure is saved as a timestamped `.png` in the `Compare_Simulations/` folder.

# Output
- `Average_All_Monomers_States_vs_Timesteps_<timestamp>.png`: Line plot visualizing
  the time evolution of the average monomer counts.

# Notes
- This function assumes that each CSV file has the same number of timesteps and that
  the first column is a shared `Timesteps` vector.
- Missing values in the simulation CSVs are automatically skipped when calculating means.

# FAIR Principles
- **Findable**: Output files are named with timestamps and saved in a centralized analysis directory.
- **Accessible**: Outputs use standard formats (CSV and PNG) that are easy to open and interpret.
- **Interoperable**: The function works directly on cleanly structured tabular data.
- **Reusable**: Modular structure allows the function to be reused or extended to other monomer species.
"""


function run_plot_all_monomer_states(directory::String, Total_Timesteps::Int)
    function process_data(file_path::String)
        data = CSV.read(file_path, DataFrame)
        mean_values = Float64[]

        for row in 1:nrow(data)
            row_data = data[row, 2:end]  # Exclude the first column (Timesteps)
            row_mean = mean(skipmissing(row_data))
            push!(mean_values, row_mean)
        end

        return mean_values
    end

    # Construct full paths to input files
    amyloid_file = joinpath(directory, "Compare_Simulations", "Appending_Amyloid_Count.csv")
    native_file = joinpath(directory, "Compare_Simulations", "Appending_Native_Count.csv")
    oligomer_file = joinpath(directory, "Compare_Simulations", "Appending_Oligomer_Count.csv")
    aggregate_file = joinpath(directory, "Compare_Simulations", "Appending_Aggregate_Count.csv")

    # Process all monomer state files
    println("Processing Amyloid data...")
    Amyloid_Mean = process_data(amyloid_file)

    println("Processing Native data...")
    Native_Mean = process_data(native_file)

    println("Processing Oligomer data...")
    Oligomer_Mean = process_data(oligomer_file)

    println("Processing Aggregate data...")
    Aggregate_Mean = process_data(aggregate_file)

    # X-axis: timesteps
    Timesteps = 0:Total_Timesteps

    # Create the plot
    plot(Timesteps, Amyloid_Mean, label="Amyloid", linewidth=2, linecolor=:orange)
    plot!(Timesteps, Native_Mean, label="Native", linewidth=2, linecolor=:green)
    plot!(Timesteps, Oligomer_Mean, label="Oligomer", linewidth=2, linecolor=:red)
    plot!(Timesteps, Aggregate_Mean, label="Aggregate", linewidth=2, linecolor=:blue)

    # Plot styling
    title!("Average Each Monomer vs Timesteps")
    xlabel!("Timesteps")
    ylabel!("Average Count")

    # Save output
    timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
    file_name = joinpath(directory, "Compare_Simulations", "Average_All_Monomers_States_vs_Timesteps_$timestamp.png")
    savefig(file_name)
    println("Combined graph saved as: $file_name")
end
