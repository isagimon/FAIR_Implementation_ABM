using DataFrames
using Dates
using CSV

"""
Module: Append_Oligomer_Clearance_Counts.jl

Aggregates per-timestep **oligomer clearance counts** across multiple simulation runs
into centralized CSV files for downstream comparison and analysis.

This script scans all `Simulation*` folders under `directory`, reads each run’s
`Oligomers_Cleared.csv`, and appends the oligomer-cleared column into a single
summary file: `Compare_Simulations/Appending_Oligomers_Clearance_Count.csv`.

Implements FAIR principles:

- **Findable**: Aggregated outputs are stored with clear, consistent filenames in `Compare_Simulations/`.
- **Accessible**: Uses CSV and DataFrames for broad tool compatibility.
- **Interoperable**: Relies on standard Julia packages and file formats.
- **Reusable**: Decomposed into small utility functions for easy modification.

Key Features:
- Creates the summary CSV (with a `Timesteps` column) if it doesn’t exist.
- Appends one column per simulation (named by its folder).
- Expects inputs produced by your ABM pipeline (`Oligomers_Cleared.csv`).

Authors: Santiago Schnell; Conner Sandefur; Isabella Gimon
Dependencies:
- CSV
- DataFrames
- Dates

License: http://www.apache.org/licenses/LICENSE-2.0
"""


timestamp = string(now(), dateformat"mm-dd-yyyy HH:MM:SS.sss")

"""
    Extract_Oligomers_Cleared_Excel_Sheets()

Scans all simulation folders under `directory` and appends their **oligomer clearance**
time series into the central summary CSV.

# Behavior
- Ensures the summary file is initialized (`Oligomer_Count_Excel`).
- Finds folders whose basename starts with `"Simulation"`.
- For each run, reads `Oligomers_Cleared.csv` and appends its oligomer-cleared column.

# Calls
- `Oligomer_Count_Excel`
- `Save_Oligomer_Count_Column`

# Global variables read
- `directory`
"""


function Extract_Oligomers_Cleared_Count_Excel_Sheets(directory, Number_Timesteps)
    Oligomer_Count_Excel(directory, Number_Timesteps)
    All_Folders = readdir(directory, join=true)
    Simulation_Folders = filter(folder -> isdir(folder) && startswith(basename(folder), "Simulation"), All_Folders)
    
    for Folder in Simulation_Folders
        Simulation_Results_File = joinpath(Folder, "Oligomers_Cleared.csv")
        println("This is Folder: $Folder")
        Save_Oligomer_Count_Column(Simulation_Results_File, Folder)
    end
    
end

"""
    Append_Oligomer_Column(Oligomer_Column, Full_Path::String)

Appends a single simulation’s **oligomer clearance** column to
`Compare_Simulations/Appending_Oligomers_Clearance_Count.csv`.

# Arguments
- `Oligomer_Column`: A vector or single-column data representing oligomer clearance counts.
- `Full_Path`: Full path to the simulation folder (used to derive the column name).

# Uses
- `Location_Oligomer_CSV_File`
- `Retrieve_Simulation_Name`

# Notes
- The appended column header is the simulation folder name (e.g., `Simulation_YYYY-MM-DD_HH-MM-SS`).

# Global variables read
- `directory`
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

Reads `Oligomers_Cleared.csv` for a given simulation and prepares the
**oligomer clearance** column for appending to the summary CSV.

# Arguments
- `Oligomer_Count_File`: Path to the simulation’s `Oligomers_Cleared.csv`.
- `Full_Path`: Full path to the simulation folder (for naming the appended column).

# Behavior
- Reads the CSV as a DataFrame.
- Extracts the **second column** (assumed to be oligomer clearance counts).
- Passes the column to `Append_Oligomer_Column`.

# Calls
- `Append_Oligomer_Column`
"""

function Save_Oligomer_Count_Column(Oligomer_Count_File, Full_Path)
    Read_File = CSV.read(Oligomer_Count_File, DataFrame)
    Oligomer_Column_Name = names(Read_File)[2]
    Oligomer_Column_Data = Read_File[!, Oligomer_Column_Name]
    Oligomer_Column_DataFrame = DataFrame(Oligomer_Column = Oligomer_Column_Data)
    Append_Oligomer_Column(Oligomer_Column_Data, Full_Path)
end

"""
    Oligomer_Count_Excel(directory::String)

Initializes the central summary file `Compare_Simulations/Appending_Oligomers_Clearance_Count.csv`
if it does not already exist.

# Arguments
- `directory`: Root directory containing all simulation subfolders.

# Behavior
- Checks for the presence of the summary CSV.
- If missing, creates it with a `Timesteps` column from `1:Number_Timesteps`.

# Calls
- `Checks_Oligomer_Excel_Present`
- `Append_Number_Timesteps`

# Global variables read
- `Number_Timesteps`
"""

function Oligomer_Count_Excel(directory, Number_Timesteps)
    if Checks_Oligomer_Excel_Present(directory) == false
        Compare_Simulation_Directory = directory * "/Compare_Simulations"
        File_Name = "Appending_Oligomers_Clearance_Count.csv"
        Complete_File_Path = joinpath(Compare_Simulation_Directory, File_Name)
        Timesteps = Append_Number_Timesteps(Number_Timesteps)
        CSV.write(Complete_File_Path, Timesteps)
    end
end

"""
    Checks_Oligomer_Excel_Present(directory::String) -> Bool

Checks whether `Compare_Simulations/Appending_Oligomers_Clearance_Count.csv` exists.

# Arguments
- `directory`: Root directory to check.

# Returns
- `false` if the file is not found; `nothing` otherwise.

# Notes
- Prints a message if the file is missing to aid debugging.
"""

function Checks_Oligomer_Excel_Present(directory)
    Compare_Simulation_Directory = directory * "/Compare_Simulations"
    CSV_File = "Appending_Oligomers_Clearance_Count.csv"
    File_Path = joinpath(Compare_Simulation_Directory, CSV_File)
    if !isfile(File_Path)
        println("File is not found: $Compare_Simulation_Directory")
        return false
    end

end

"""
    Append_Number_Timesteps() -> DataFrame

Creates a DataFrame with a single column `Timesteps` from `1:Number_Timesteps`.

# Returns
- `DataFrame` of the form:
    Timesteps
    0
    1
    …
    N

# Global variables read
- `Number_Timesteps`
"""


function Append_Number_Timesteps(Number_Timesteps) 
    return timesteps = DataFrame(Timesteps = 0:Number_Timesteps)
end

"""
    Location_Oligomer_CSV_File(directory::String) -> String

Returns the absolute path to `Compare_Simulations/Appending_Oligomers_Clearance_Count.csv`.

# Arguments
- `directory`: Root folder of the simulation project.

# Returns
- `String` path to the oligomer clearance summary CSV.
"""


function Location_Oligomer_CSV_File(directory) 
    return joinpath(directory, "Compare_Simulations", "Appending_Oligomers_Clearance_Count.csv")
end

"""
    Retrieve_Simulation_Name(Path::String) -> String

Extracts the simulation folder name from a full directory path using regex.

# Arguments
- `Path`: Full path to the simulation directory.

# Returns
- The simulation folder name (e.g., `"Simulation_2025-05-18_14-00-00"`).

# Notes
- Matches the substring starting at `"Simulation"` to the end of `Path`.
"""

function Retrieve_Simulation_Name(Path)
    Simulation_Name = match(r"(Simulation.*)$", Path)
    return Simulation_Name.match
end

