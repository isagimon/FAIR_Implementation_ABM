# Load required libraries
using CSV             # For reading and writing CSV files
using DataFrames      # For handling tabular data in a structured format

# ------------------------------------------------------------------------
# Set up project paths and load parameter values from a central CSV file
# ------------------------------------------------------------------------

# Define the absolute path to the root of your project directory
basic_directory = "/Users/isabellagimon/Desktop/FAIR_Implementation_ABM"

# Construct the full path to the input parameter file located in the root directory
csv_file = "$basic_directory/Input_Parameters_Analysis.csv"

# Read the CSV file into a DataFrame where each row defines a key-value pair of parameters
parameters_df = CSV.read(csv_file, DataFrame)

# Extract the simulation data directory path from the parameter table
directory = parameters_df[parameters_df.Input_Parameters .== "Directory", :Values][1]

# Extract the number of timesteps to analyze and convert it from string to integer
Number_Timesteps = parse(Int, parameters_df[parameters_df.Input_Parameters .== "Total_Timesteps", :Values][1])

# Extract the total number of monomers (used optionally in some analysis) and convert to integer
Total_Number_Monomers = parse(Int, parameters_df[parameters_df.Input_Parameters .== "Total_Number_Monomers", :Values][1])

# Print the loaded parameters for verification and debugging
println("Directory: ", directory)
println("Number of Timesteps: ", Number_Timesteps)
println("Total Number of Monomers: ", Total_Number_Monomers)

# ------------------------------------------------------------------------
# Load analysis scripts that define and implement the core analysis logic
# ------------------------------------------------------------------------

# Load the script that extracts native and amyloid monomer counts
include("$basic_directory/Analysis/Append_Amyloid_and_Native.jl")

# Load the script that extracts aggregate and oligomer counts
include("$basic_directory/Analysis/Append_Aggregate_and_Oligomer.jl")

# Load the script that plots average monomer state counts over time
include("$basic_directory/Analysis/Average_All_Monomers_vs_Timesteps.jl")

# ------------------------------------------------------------------------
# Run all analysis functions sequentially
# ------------------------------------------------------------------------

# Append native and amyloid monomer count data across all simulations
run_append_amyloid_and_native(directory, Number_Timesteps)

# Append aggregate and oligomer count data into combined CSVs
run_process_aggregate_excel_sheets(directory, Number_Timesteps)

# Generate and save a line plot of average monomer state counts over time
run_plot_all_monomer_states(directory, Number_Timesteps)

