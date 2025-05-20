using XLSX
using DataFrames
using Dates
using FilePaths
using CSV

"""
Module: Append_Native_Amyloid_Counts.jl

Processes and aggregates native and amyloid-prone monomer count data from multiple simulation folders, 
compiling results into centralized summary CSV files for cross-simulation comparison.

This script automates the extraction of relevant time-series data from `Native_and_Amyloid_Count_Results.csv` 
files generated during simulations. It isolates the columns corresponding to native and amyloid-prone monomers 
and appends them into collective CSVs housed in the `Compare_Simulations` subdirectory.

Implements FAIR principles:

- **Findable**: Output files are stored in clearly labeled, centralized CSV files for each species type.
- **Accessible**: Data is stored using standard CSV format and can be easily opened and analyzed with common tools.
- **Interoperable**: Compatible with the Julia ecosystem via CSV.jl and DataFrames.jl.
- **Reusable**: Modular functions allow for easy adjustment of aggregation targets or simulation structure.

Key Features:
- Dynamically scans all simulation folders with names starting in "Simulation".
- Initializes and appends to `Appending_Amyloid_Count.csv` and `Appending_Native_Count.csv`.
- Ensures output CSVs are synchronized with the number of simulation timesteps.
- Useful for tracking and comparing native and amyloid monomer counts across many independent runs.

Authors: Santiago Schnell, Conner Sandefur, Isabella Gimon  
Dependencies:
- CSV  
- DataFrames  
- XLSX  
- Dates  
- FilePaths  

License: http://www.apache.org/licenses/LICENSE-2.0
"""

"""
    run_append_amyloid_and_native(directory::String, Number_Timesteps::Int)

Main function that processes simulation results for amyloid-prone and native monomer counts. 

# Arguments
- `directory::String`: Path to the top-level directory containing simulation folders.
- `Number_Timesteps::Int`: Number of timesteps per simulation.

# Behavior
- Scans all folders beginning with "Simulation".
- Extracts amyloid and native monomer counts from each simulation's result file.
- Appends each result as a new column into their respective summary CSVs.

# Notes
- The expected input file per simulation is `Native_and_Amyloid_Count_Results.csv`.
"""
# Main execution function
function run_append_amyloid_and_native(directory::String, Number_Timesteps::Int)
    # Initialize arrays to store averages from timestep 0 to Number_Timesteps (inclusive)
    Avg_Oligomer_Count = zeros(Int, Number_Timesteps + 1, 2) #TIMESTEP HAS TO BEGIN AT ZERO 
    Avg_Fibril_Count = zeros(Int, Number_Timesteps + 1, 2) #TIMESTEP HAS TO BEGIN AT ZERO 

    timestamp = string(now(), dateformat"mm-dd-yyyy HH:MM:SS.sss")

    # Core function to process Excel sheets
    function Extract_Aggregate_Count_Excel_Sheets()
        Amyloid_Count_Excel(directory)
        Native_Count_Excel(directory)
        All_Folders = readdir(directory, join=true)
        Simulation_Folders = filter(folder -> isdir(folder) && startswith(basename(folder), "Simulation"), All_Folders)

        for Folder in Simulation_Folders
            Simulation_Results_File = joinpath(Folder, "Native_and_Amyloid_Count_Results.csv")
            println("This is Folder: $Folder")
            Save_Amyloid_Count_Column(Simulation_Results_File, Folder)
            Save_Native_Count_Column(Simulation_Results_File, Folder)
        end
    end

    # Run the extraction process
    Extract_Aggregate_Count_Excel_Sheets()
end


"""
    Save_Amyloid_Count_Column(Amorphous_Count_File, Full_Path)

Extracts the amyloid monomer count column (3rd column) from a simulation results file
and formats it as a DataFrame.

# Arguments
- `Amorphous_Count_File`: Path to the CSV file `Native_and_Amyloid_Count_Results.csv`.
- `Full_Path`: Full path to the simulation folder.

# Calls
- `Append_Amyloid_Column`
"""

# All helper functions remain the same but now use `directory` and `Number_Timesteps` dynamically
function Save_Amyloid_Count_Column(Amorphous_Count_File, Full_Path)
    Read_File = CSV.read(Amorphous_Count_File, DataFrame)
    Aggregate_Column_Name = names(Read_File)[3]
    Aggregate_Column_Data = Read_File[!, Aggregate_Column_Name]
    Aggregate_Column_DataFrame = DataFrame(Aggregate_Column = Aggregate_Column_Data)
    Append_Amyloid_Column(Aggregate_Column_DataFrame, Full_Path)
end

"""
    Append_Amyloid_Column(Aggregate_Column::DataFrame, Full_Path::String)

Appends amyloid monomer count data as a new column to the central summary CSV.

# Arguments
- `Aggregate_Column`: DataFrame containing amyloid monomer counts.
- `Full_Path`: Path to the simulation folder used to name the column.

# Notes
- Uses the folder name (e.g., Simulation1) as the new column header.
"""

function Append_Amyloid_Column(Aggregate_Column, Full_Path)
    Simulation_Name = Retrieve_Simulation_Name(Full_Path)
    println("This is Simulation_Name: $Simulation_Name")
    Amorphous_CSV_File = Location_Amyloid_CSV_File(directory)
    Read_Amorphous_CSV_File = CSV.read(Amorphous_CSV_File, DataFrame)
    rename!(Aggregate_Column, names(Aggregate_Column)[1] => Simulation_Name)
    Read_Amorphous_CSV_File[!, Simulation_Name] = Aggregate_Column[:, Simulation_Name]
    CSV.write(Amorphous_CSV_File, Read_Amorphous_CSV_File)
end

"""
    Append_Native_Column(Oligomer_Column::Vector, Full_Path::String)

Appends native monomer count data to the central summary CSV.

# Arguments
- `Oligomer_Column`: Vector of native monomer counts.
- `Full_Path`: Path to the simulation folder used to label the column.
"""

function Append_Native_Column(Oligomer_Column, Full_Path)
    Simulation_Name = Retrieve_Simulation_Name(Full_Path)
    println("This is Simulation_Name: $Simulation_Name")
    Oligomer_CSV_File = Location_Native_CSV_File(directory)
    Read_Oligomer_CSV_File = CSV.read(Oligomer_CSV_File, DataFrame)
    Read_Oligomer_CSV_File[!, Simulation_Name] = Oligomer_Column
    CSV.write(Oligomer_CSV_File, Read_Oligomer_CSV_File)
end

"""
    Save_Native_Count_Column(Oligomer_Count_File::String, Full_Path::String)

Extracts the native monomer count column (2nd column) from the simulation result CSV file 
(`Native_and_Amyloid_Count_Results.csv`) and prepares it for appending to the summary CSV.

# Arguments
- `Oligomer_Count_File`: Path to the simulation results CSV.
- `Full_Path`: Full path to the simulation folder used to name the new column.

# Calls
- `Append_Native_Column`
"""
function Save_Native_Count_Column(Oligomer_Count_File, Full_Path)
    Read_File = CSV.read(Oligomer_Count_File, DataFrame)
    Oligomer_Column_Name = names(Read_File)[2]
    Oligomer_Column_Data = Read_File[!, Oligomer_Column_Name]
    Oligomer_Column_DataFrame = DataFrame(Oligomer_Column = Oligomer_Column_Data)
    Append_Native_Column(Oligomer_Column_Data, Full_Path)
end

"""
    Amyloid_Count_Excel(directory::String)

Initializes the amyloid summary CSV (`Appending_Amyloid_Count.csv`) with a "Timesteps" column 
if the file does not already exist.

# Arguments
- `directory`: Path to the directory containing `Compare_Simulations`.

# Global Dependencies
- Assumes `Number_Timesteps` is accessible in scope.

# Calls
- `Append_Number_Timesteps`
"""
function Amyloid_Count_Excel(directory)
    if !Checks_Amorphous_Excel_Present(directory)
        Compare_Simulation_Directory = directory * "/Compare_Simulations"
        File_Name = "Appending_Amyloid_Count.csv"
        Complete_File_Path = joinpath(Compare_Simulation_Directory, File_Name)
        Timesteps = Append_Number_Timesteps(Number_Timesteps)
        CSV.write(Complete_File_Path, Timesteps)
    end
end

"""
    Checks_Amorphous_Excel_Present(directory::String) -> Bool

Checks whether the amyloid summary file (`Appending_Amyloid_Count.csv`) exists.

# Arguments
- `directory`: Directory path where the file should be located.

# Returns
- `true` if file exists; otherwise `false` with a printed warning.
"""
function Checks_Amorphous_Excel_Present(directory)
    Compare_Simulation_Directory = directory * "/Compare_Simulations"
    CSV_File = "Appending_Amyloid_Count.csv"
    File_Path = joinpath(Compare_Simulation_Directory, CSV_File)
    if !isfile(File_Path)
        println("File is not found: $Compare_Simulation_Directory")
        return false
    end
end

"""
    Native_Count_Excel(directory::String)

Initializes the native monomer summary CSV (`Appending_Native_Count.csv`) with a "Timesteps" column 
if the file does not already exist.

# Arguments
- `directory`: Path to the directory containing `Compare_Simulations`.

# Global Dependencies
- Assumes `Number_Timesteps` is accessible in scope.

# Calls
- `Append_Number_Timesteps`
"""
function Native_Count_Excel(directory)
    if !Checks_Oligomer_Excel_Present(directory)
        Compare_Simulation_Directory = directory * "/Compare_Simulations"
        File_Name = "Appending_Native_Count.csv"
        Complete_File_Path = joinpath(Compare_Simulation_Directory, File_Name)
        Timesteps = Append_Number_Timesteps(Number_Timesteps)
        CSV.write(Complete_File_Path, Timesteps)
    end
end

"""
    Checks_Oligomer_Excel_Present(directory::String) -> Bool

Checks whether the native monomer summary file (`Appending_Native_Count.csv`) exists.

# Arguments
- `directory`: Directory path where the file should be located.

# Returns
- `true` if file exists; otherwise `false` with a printed warning.
"""
function Checks_Oligomer_Excel_Present(directory)
    Compare_Simulation_Directory = directory * "/Compare_Simulations"
    CSV_File = "Appending_Native_Count.csv"
    File_Path = joinpath(Compare_Simulation_Directory, CSV_File)
    if !isfile(File_Path)
        println("File is not found: $Compare_Simulation_Directory")
        return false
    end
end

"""
    Append_Number_Timesteps(Number_Timesteps::Int) -> DataFrame

Creates a DataFrame with a "Timesteps" column ranging from 1 to `Number_Timesteps`.

# Arguments
- `Number_Timesteps`: Total number of simulation steps.

# Returns
- DataFrame with a single "Timesteps" column.
"""
function Append_Number_Timesteps(Number_Timesteps)
    return DataFrame(Timesteps = 1:Number_Timesteps)
end

"""
    Location_Amyloid_CSV_File(directory::String) -> String

Constructs the full file path to the amyloid monomer summary CSV file.

# Arguments
- `directory`: Root path of the simulation folder.

# Returns
- Full file path to `Appending_Amyloid_Count.csv`.
"""
function Location_Amyloid_CSV_File(directory)
    return joinpath(directory, "Compare_Simulations", "Appending_Amyloid_Count.csv")
end

"""
    Location_Native_CSV_File(directory::String) -> String

Constructs the full file path to the native monomer summary CSV file.

# Arguments
- `directory`: Root path of the simulation folder.

# Returns
- Full file path to `Appending_Native_Count.csv`.
"""
function Location_Native_CSV_File(directory)
    return joinpath(directory, "Compare_Simulations", "Appending_Native_Count.csv")
end

"""
    Retrieve_Simulation_Name(Path::String) -> String

Extracts the simulation folder name (e.g., "Simulation1") from a full file path.

# Arguments
- `Path`: Full file path including the simulation folder.

# Returns
- Name of the simulation folder as a string.
"""
function Retrieve_Simulation_Name(Path)
    Simulation_Name = match(r"(Simulation.*)$", Path)
    return Simulation_Name.match
end

function Append_Number_Timesteps(Number_Timesteps)
    return DataFrame(Timesteps = 0:Number_Timesteps)
end
