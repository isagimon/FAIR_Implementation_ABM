using XLSX
using DataFrames
using Dates
using FilePaths
using CSV

"""
Module: Append_Aggregate_and_Oligomer.jl

Processes and compiles oligomer and fibril count data across multiple simulation runs,
storing results into centralized summary CSV files for downstream analysis.

This script automates the collection of time-series aggregation data by reading individual
`Oligomer_and_Aggregate_Count_Results.csv` files and extracting relevant columns for comparison.

Implements FAIR principles:

- **Findable**: Aggregated outputs are stored in clearly named CSV files within a standardized directory structure (`Compare_Simulations`).
- **Accessible**: Compatible with CSV and DataFrame formats for easy access and sharing.
- **Interoperable**: Uses standard Julia packages and file formats.
- **Reusable**: Modular design allows easy modification for different metrics or output formats.

Key Features:
- Automatically scans and processes all simulation folders starting with “Simulation”.
- Extracts and appends fibril (amorphous aggregate) and oligomer counts to separate summary files.
- Verifies and creates summary CSV files if not already present.
- Supports user-defined number of timesteps to ensure synchronized data collection.

Authors: Santiago Schnell; Conner Sandefur; Isabella Gimon  
Dependencies:
- CSV  
- DataFrames  
- Dates  
- XLSX  
- FilePaths  

License: http://www.apache.org/licenses/LICENSE-2.0
"""


"""
    run_process_aggregate_excel_sheets(directory::String, Number_Timesteps::Int)

Processes all simulation folders under `directory` and appends oligomer and fibril
(aggregate) counts into centralized CSV files.

# Arguments
- `directory`: Root directory containing all simulation subfolders.
- `Number_Timesteps`: Total number of simulation timesteps.

# Behavior
- Scans all folders starting with "Simulation".
- Reads `Oligomer_and_Aggregate_Count_Results.csv` from each.
- Extracts and appends fibril and oligomer columns to:
  - `Appending_Aggregate_Count.csv`
  - `Appending_Oligomer_Count.csv`

# Calls
- `Save_Amorphous_Count_Column`
- `Save_Oligomer_Count_Column`
"""

function run_process_aggregate_excel_sheets(directory::String, Number_Timesteps::Int)
    Amorphous_Count_Excel(directory, Number_Timesteps)
    Oligomer_Count_Excel(directory, Number_Timesteps)
    All_Folders = readdir(directory, join=true)
    Simulation_Folders = filter(folder -> isdir(folder) && startswith(basename(folder), "Simulation"), All_Folders)
    
    for Folder in Simulation_Folders
        #Line_Graph_Folder = joinpath(Folder, "Line_Graphs")
        Simulation_Results_File = joinpath(Folder, "Oligomer_and_Aggregate_Count_Results.csv")
        println("This is Folder: $Folder")
        Save_Amorphous_Count_Column(Simulation_Results_File, Folder)
        #println("This is the Amorphous_Count_File: ",Amorphous_Count_File)
        Save_Oligomer_Count_Column(Simulation_Results_File, Folder)
        #println("This is the Oligomer_Count_File: ",Oligomer_Count_File)
    end
end

"""
    Save_Amorphous_Count_Column(Amorphous_Count_File::String, Full_Path::String)

Extracts the third column (fibril/aggregate count) from a simulation output CSV
and prepares it for appending to the fibril summary CSV.

# Arguments
- `Amorphous_Count_File`: Full path to the simulation results CSV.
- `Full_Path`: Path to the simulation folder (used for naming).

# Calls
- `Append_Amorphous_Column`

# Notes
- Assumes the third column contains fibril counts.
"""

function Save_Amorphous_Count_Column(Amorphous_Count_File, Full_Path)
    Read_File = CSV.read(Amorphous_Count_File, DataFrame)
    Aggregate_Column_Name = names(Read_File)[3]
    Aggregate_Column_Data = Read_File[!, Aggregate_Column_Name]
    Aggregate_Column_DataFrame = DataFrame(Aggregate_Column = Aggregate_Column_Data)
    Append_Amorphous_Column(Aggregate_Column_DataFrame, Full_Path)
    #println("This is Second Column: $Aggregate_Column_DataFrame")
end

"""
    Append_Amorphous_Column(Aggregate_Column::DataFrame, Full_Path::String)

Appends a fibril count column for one simulation to the aggregated fibril CSV.

# Arguments
- `Aggregate_Column`: A DataFrame with one column of fibril counts.
- `Full_Path`: Path to the simulation folder for retrieving the simulation name.

# Uses
- `Location_Fibril_CSV_File`
- `Retrieve_Simulation_Name`

# Notes
- Renames the column to the simulation folder name before appending.
"""


function Append_Amorphous_Column(Aggregate_Column, Full_Path)
    Simulation_Name = Retrieve_Simulation_Name(Full_Path)
    println("This is Simulation_Name: $Simulation_Name")
    Amorphous_CSV_File = Location_Fibril_CSV_File(directory)
    Read_Amorphous_CSV_File = CSV.read(Amorphous_CSV_File, DataFrame)
    rename!(Aggregate_Column, names(Aggregate_Column)[1] => Simulation_Name)
    Read_Amorphous_CSV_File[!, Simulation_Name] = Aggregate_Column[:, Simulation_Name]
    CSV.write(Amorphous_CSV_File, Read_Amorphous_CSV_File)
end

"""
    Append_Oligomer_Column(Oligomer_Column::Vector, Full_Path::String)

Appends oligomer data from a single simulation to `Appending_Oligomer_Count.csv`.

# Arguments
- `Oligomer_Column`: A vector of oligomer count values.
- `Full_Path`: Path used to derive the simulation folder name.

# Uses
- `Location_Oligomer_CSV_File`
- `Retrieve_Simulation_Name`

# Notes
- Appends the column using the simulation folder name as header.
"""


function Append_Oligomer_Column(Oligomer_Column, Full_Path)
    Simulation_Name = Retrieve_Simulation_Name(Full_Path)
    println("This is Simulation_Name: $Simulation_Name")
    Oligomer_CSV_File = Location_Oligomer_CSV_File(directory)
    Read_Oligomer_CSV_File = CSV.read(Oligomer_CSV_File, DataFrame)
    Read_Oligomer_CSV_File[!, Simulation_Name] = Oligomer_Column
    CSV.write(Oligomer_CSV_File, Read_Oligomer_CSV_File)

end

"""
    Save_Oligomer_Count_Column(Oligomer_Count_File::String, Full_Path::String)

Extracts the second column (oligomer count) from a simulation results CSV and
prepares it for appending to the oligomer summary CSV.

# Arguments
- `Oligomer_Count_File`: Path to the results CSV.
- `Full_Path`: Path to the simulation folder (used for naming).

# Calls
- `Append_Oligomer_Column`

# Notes
- Assumes the second column contains oligomer count data.
"""

function Save_Oligomer_Count_Column(Oligomer_Count_File, Full_Path)
    Read_File = CSV.read(Oligomer_Count_File, DataFrame)
    Oligomer_Column_Name = names(Read_File)[2]
    Oligomer_Column_Data = Read_File[!, Oligomer_Column_Name]
    Oligomer_Column_DataFrame = DataFrame(Oligomer_Column = Oligomer_Column_Data)
    Append_Oligomer_Column(Oligomer_Column_Data, Full_Path)
end

"""
    Amorphous_Count_Excel(directory::String, Number_Timesteps::Int)

Initializes the fibril summary file `Appending_Aggregate_Count.csv` if it does not exist.

# Arguments
- `directory`: Root directory containing simulation data.
- `Number_Timesteps`: Total number of simulation timesteps.

# Calls
- `Checks_Amorphous_Excel_Present`
- `Append_Number_Timesteps`
"""

function Amorphous_Count_Excel(directory, Number_Timesteps)
    if Checks_Amorphous_Excel_Present(directory) == false
        Compare_Simulation_Directory = directory * "/Compare_Simulations"
        File_Name = "Appending_Aggregate_Count.csv"
        Complete_File_Path = joinpath(Compare_Simulation_Directory, File_Name)
        Timesteps = Append_Number_Timesteps(Number_Timesteps)
        CSV.write(Complete_File_Path, Timesteps)
    end

end

"""
    Checks_Amorphous_Excel_Present(directory::String) -> Bool

Checks whether the fibril summary file (`Appending_Aggregate_Count.csv`) exists.

# Arguments
- `directory`: Root directory to check.

# Returns
- `false` if file is missing; `nothing` otherwise.

# Notes
- Prints a message if file is not found.
"""


function Checks_Amorphous_Excel_Present(directory)
    Compare_Simulation_Directory = directory * "/Compare_Simulations"
    CSV_File = "Appending_Aggregate_Count.csv"
    File_Path = joinpath(Compare_Simulation_Directory, CSV_File)
    if !isfile(File_Path)
        println("File is not found: $Compare_Simulation_Directory")
        return false
    end

end

"""
    Oligomer_Count_Excel(directory::String, Number_Timesteps::Int)

Initializes `Appending_Oligomer_Count.csv` if it doesn’t already exist.

# Arguments
- `directory`: Root folder for simulation results.
- `Number_Timesteps`: Total number of simulation timesteps.

# Calls
- `Checks_Oligomer_Excel_Present`
- `Append_Number_Timesteps`
"""


function Oligomer_Count_Excel(directory, Number_Timesteps)
    if Checks_Oligomer_Excel_Present(directory) == false
        Compare_Simulation_Directory = directory * "/Compare_Simulations"
        File_Name = "Appending_Oligomer_Count.csv"
        Complete_File_Path = joinpath(Compare_Simulation_Directory, File_Name)
        Timesteps = Append_Number_Timesteps(Number_Timesteps)
        CSV.write(Complete_File_Path, Timesteps)
    end
end

"""
    Checks_Oligomer_Excel_Present(directory::String) -> Bool

Checks for the existence of `Appending_Oligomer_Count.csv`.

# Arguments
- `directory`: Root directory of simulation outputs.

# Returns
- `false` if file is not found; `nothing` otherwise.

# Notes
- Helps avoid overwriting the summary file if it already exists.
"""

function Checks_Oligomer_Excel_Present(directory)
    Compare_Simulation_Directory = directory * "/Compare_Simulations"
    CSV_File = "Appending_Oligomer_Count.csv"
    File_Path = joinpath(Compare_Simulation_Directory, CSV_File)
    if !isfile(File_Path)
        println("File is not found: $Compare_Simulation_Directory")
        return false
    end

end

"""
    Append_Number_Timesteps(Number_Timesteps::Int) -> DataFrame

Creates a DataFrame with one column labeled `Timesteps`, ranging from 0 to `Number_Timesteps`.

# Arguments
- `Number_Timesteps`: The number of time steps to record.

# Returns
- A DataFrame of the form:
    Timesteps
    0
    1
    ...
    N
"""


function Append_Number_Timesteps(Number_Timesteps) 
    return timesteps = DataFrame(Timesteps = 0:Number_Timesteps)
end

"""
    Location_Fibril_CSV_File(directory::String) -> String

Returns the full path to the central `Appending_Aggregate_Count.csv` file.

# Arguments
- `directory`: Root folder of the simulation project.

# Returns
- A string path to the fibril summary file.
"""

function Location_Fibril_CSV_File(directory) 
    return joinpath(directory, "Compare_Simulations", "Appending_Aggregate_Count.csv")
end

"""
    Location_Oligomer_CSV_File(directory::String) -> String

Returns the full path to the central `Appending_Oligomer_Count.csv` file.

# Arguments
- `directory`: Root folder of the simulation project.

# Returns
- A string path to the oligomer summary file.
"""

function Location_Oligomer_CSV_File(directory) 
    return joinpath(directory, "Compare_Simulations", "Appending_Oligomer_Count.csv")
end

"""
    Retrieve_Simulation_Name(Path::String) -> String

Extracts the simulation folder name from a full directory path using regex.

# Arguments
- `Path`: Full path to the simulation directory.

# Returns
- The simulation folder name (e.g., "Simulation_2025-05-18_14-00-00").
"""

function Retrieve_Simulation_Name(Path)
    Simulation_Name = match(r"(Simulation.*)$", Path)
    return Simulation_Name.match
end