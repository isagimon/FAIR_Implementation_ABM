using CSV
using Statistics
using Plots
using DataFrames
using Dates

"""
    run_plot_all_monomer_states(directory::String, Total_Timesteps::Int)

Generates a line plot showing the average number of monomers in each state (Amyloid, Native, Oligomer, Aggregate)
over time, across multiple simulation runs. Reads data from summary CSVs stored in the `Compare_Simulations` folder 
and computes row-wise averages to produce time-series curves.

# Arguments
- `directory::String`: Path to the base directory containing the `Compare_Simulations` folder.
- `Total_Timesteps::Int`: Total number of timesteps to plot on the X-axis.

# Behavior
- Reads four files: `Appending_Amyloid_Count.csv`, `Appending_Native_Count.csv`,
  `Appending_Oligomer_Count.csv`, and `Appending_Fibril_Count.csv`.
- Calculates the mean value at each timestep across all simulations for each monomer type.
- Produces a single plot displaying the average value of each species over time.
- Saves the resulting plot as a PNG file with a timestamped filename.

# Output
- A time-resolved PNG plot saved to `Compare_Simulations/` in the specified directory.
"""
function run_plot_all_monomer_states(directory::String, Total_Timesteps::Int)

    """
        process_data(file_path::String) -> Vector{Float64}

    Reads a CSV file and computes the row-wise average across all simulation columns,
    skipping any missing values.

    # Arguments
    - `file_path::String`: Full path to the CSV file to process.

    # Returns
    - A vector of Float64 values representing the average count per timestep.
    """
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
    aggregate_file = joinpath(directory, "Compare_Simulations", "Appending_Fibril_Count.csv")

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
    Timesteps = 1:Total_Timesteps

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
