using CSV
using DataFrames

# Absolute path to your project folder
basic_directory = "/Users/isabellagimon/Desktop/FAIR_Implementation_ABM"

# Correct path to the parameter file in the root directory
csv_file = "$basic_directory/Input_Parameters_Analysis.csv"

# Read parameters from the CSV
parameters_df = CSV.read(csv_file, DataFrame)

# Extract values from the table
directory = parameters_df[parameters_df.Input_Parameters .== "Directory", :Values][1]
Number_Timesteps = parse(Int, parameters_df[parameters_df.Input_Parameters .== "Total_Timesteps", :Values][1])
Total_Number_Monomers = parse(Int, parameters_df[parameters_df.Input_Parameters .== "Total_Number_Monomers", :Values][1])

println("Directory: ", directory)
println("Number of Timesteps: ", Number_Timesteps)
println("Total Number of Monomers: ", Total_Number_Monomers)

# Include analysis scripts located inside the `analysis/` folder
include("$basic_directory/Analysis/Append_Amyloid_and_Native.jl")
include("$basic_directory/Analysis/Append_Aggregate_and_Oligomer.jl")
include("$basic_directory/Analysis/Average_All_Monomers_vs_Timesteps.jl")

# Run the analysis
run_append_amyloid_and_native(directory, Number_Timesteps)
run_process_aggregate_excel_sheets(directory, Number_Timesteps)
run_plot_all_monomer_states(directory, Number_Timesteps)
