# Run_All_Analysis_Scripts.jl
#
# Runs the post-processing / analysis pipeline over a set of simulation output folders
# located under `Data_Collection/` (or another directory configured in
# Input_Parameters_Analysis.csv).
#
# Usage:
#   julia Run_All_Analysis_Scripts.jl
#

using Pkg
Pkg.activate(joinpath(@__DIR__, "Analysis"))
Pkg.instantiate()

using CSV
using DataFrames

# ------------------------------------------------------------------------
# Project paths and analysis parameters
# ------------------------------------------------------------------------

basic_directory = @__DIR__
csv_file = joinpath(basic_directory, "Input_Parameters_Analysis.csv")

parameters_df = CSV.read(csv_file, DataFrame)

clean_names = [strip(replace(String(name), '\ufeff' => "")) for name in names(parameters_df)]
rename!(parameters_df, Symbol.(clean_names))

# Directory containing simulation folders
directory_value = string(parameters_df[parameters_df.Input_Parameters .== "Directory", :Values][1])
directory = isabspath(directory_value) ? directory_value : joinpath(basic_directory, directory_value)

# Timesteps and total monomers
Number_Timesteps = parse(Int, string(parameters_df[parameters_df.Input_Parameters .== "Total_Timesteps", :Values][1]))
Total_Number_Monomers = parse(Int, string(parameters_df[parameters_df.Input_Parameters .== "Total_Number_Monomers", :Values][1]))

# ------------------------------------------------------------------------
# Run analysis scripts
# ------------------------------------------------------------------------

include(joinpath(basic_directory, "Analysis", "Append_AggregateProne_and_Native.jl"))
run_append_AggregateProne_and_native(directory, Number_Timesteps)

include(joinpath(basic_directory, "Analysis", "Append_Aggregate_and_Oligomer.jl"))
run_process_aggregate_excel_sheets(directory, Number_Timesteps)

include(joinpath(basic_directory, "Analysis", "Average_All_Monomers_vs_Timesteps.jl"))
run_plot_all_monomer_states(directory, Number_Timesteps)

include(joinpath(basic_directory, "Analysis", "Append_Oligomer_Clearance_Data.jl"))
Extract_Oligomers_Cleared_Count_Excel_Sheets(directory, Number_Timesteps)

include(joinpath(basic_directory, "Analysis", "Average_Oligomers_Cleared_vs_Timesteps.jl"))
Import_Data(directory, Number_Timesteps)
