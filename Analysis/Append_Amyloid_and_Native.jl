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
    Amyloid_Count_Excel(directory, Number_Timesteps)
    Native_Count_Excel(directory, Number_Timesteps)
    All_Folders = readdir(directory, join=true)
    Simulation_Folders = filter(folder -> isdir(folder) && startswith(basename(folder), "Simulation"), All_Folders)
    
    for Folder in Simulation_Folders
        #Line_Graph_Folder = joinpath(Folder, "Line_Graphs")
        Simulation_Results_File = joinpath(Folder, "Native_and_Amyloid_Count_Results.csv")
        println("This is Folder: $Folder")
        Save_Amyloid_Count_Column(Simulation_Results_File, Folder)
        #println("This is the Amorphous_Count_File: ",Amorphous_Count_File)
        Save_Native_Count_Column(Simulation_Results_File, Folder)
        #println("This is the Oligomer_Count_File: ",Oligomer_Count_File)
    end
end

"""
    Save_Amyloid_Count_Column(Amorphous_Count_File, Full_Path)

Extracts the amyloid monomer count column from a simulation result CSV file and 
prepares it for aggregation by wrapping it in a standardized DataFrame.

# Arguments
- `Amorphous_Count_File`: The full path to the CSV file containing simulation results.
- `Full_Path`: Full path to the current simulation folder.

# Calls
- `Append_Amyloid_Column()` to insert the extracted column into the aggregated results.

# Notes
- The function assumes the amyloid count is in the third column of the result CSV.
"""
function Save_Amyloid_Count_Column(Amorphous_Count_File, Full_Path)
    Read_File = CSV.read(Amorphous_Count_File, DataFrame)
    Aggregate_Column_Name = names(Read_File)[3]
    Aggregate_Column_Data = Read_File[!, Aggregate_Column_Name]
    Aggregate_Column_DataFrame = DataFrame(Aggregate_Column = Aggregate_Column_Data)
    Append_Amyloid_Column(Aggregate_Column_DataFrame, Full_Path)
end

"""
    Append_Amyloid_Column(Aggregate_Column, Full_Path)

Appends amyloid monomer count data for a single simulation to the aggregated 
summary CSV file `Appending_Amyloid_Count.csv`.

# Arguments
- `Aggregate_Column`: A DataFrame containing a single column of amyloid counts.
- `Full_Path`: Path to the simulation folder, used to extract the simulation name.

# Global variables used
- `directory` (must be defined externally)

# Notes
- If the output CSV already contains a column for the simulation, it will be overwritten.
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
    Append_Native_Column(Oligomer_Column, Full_Path)

Appends a single simulation’s native monomer data to the 
`Appending_Native_Count.csv` file in `Compare_Simulations/`.

# Arguments
- `Oligomer_Column`: Vector or column of native monomer counts.
- `Full_Path`: Full path to the simulation folder, used to retrieve simulation name.

# Global variables used
- `directory`

# Notes
- The appended column is not renamed. It is added with the simulation name as its header.
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
    Save_Native_Count_Column(Oligomer_Count_File, Full_Path)

Extracts the native monomer column from a simulation result CSV file and passes it
to be appended to the native monomer master CSV.

# Arguments
- `Oligomer_Count_File`: Path to CSV with native and amyloid monomer counts.
- `Full_Path`: Path to simulation folder used for naming.

# Calls
- `Append_Native_Column()`

# Notes
- Assumes native monomer count is in the second column of the file.
"""

function Save_Native_Count_Column(Oligomer_Count_File, Full_Path) #NATIVE COUNT
    Read_File = CSV.read(Oligomer_Count_File, DataFrame)
    Oligomer_Column_Name = names(Read_File)[2]
    Oligomer_Column_Data = Read_File[!, Oligomer_Column_Name]
    Oligomer_Column_DataFrame = DataFrame(Oligomer_Column = Oligomer_Column_Data)
    Append_Native_Column(Oligomer_Column_Data, Full_Path)
end

"""
    Amyloid_Count_Excel(directory, Number_Timesteps)

Creates `Appending_Amyloid_Count.csv` if it does not already exist.

# Arguments
- `directory`: Path to the root simulation directory.
- `Number_Timesteps`: Number of simulation steps to include in the file.

# Calls
- `Checks_Amorphous_Excel_Present()`
- `Append_Number_Timesteps()`

# Notes
- The CSV is initialized with a `Timesteps` column from 0 to `Number_Timesteps`.
"""

function Amyloid_Count_Excel(directory, Number_Timesteps) #AMYLOID
    if Checks_Amorphous_Excel_Present(directory) == false
        Compare_Simulation_Directory = directory * "/Compare_Simulations"
        File_Name = "Appending_Amyloid_Count.csv"
        Complete_File_Path = joinpath(Compare_Simulation_Directory, File_Name)
        Timesteps = Append_Number_Timesteps(Number_Timesteps)
        CSV.write(Complete_File_Path, Timesteps)
    end

end

"""
    Checks_Amorphous_Excel_Present(directory) -> Bool

Checks whether `Appending_Amyloid_Count.csv` exists in `Compare_Simulations/`.

# Arguments
- `directory`: Root path to simulation data.

# Returns
- `false` if the file is missing; otherwise returns `nothing`.

# Notes
- Prints a message if the file is not found.
"""

function Checks_Amorphous_Excel_Present(directory) #AMYLOID
    Compare_Simulation_Directory = directory * "/Compare_Simulations"
    CSV_File = "Appending_Amyloid_Count.csv"
    File_Path = joinpath(Compare_Simulation_Directory, CSV_File)
    if !isfile(File_Path)
        println("File is not found: $Compare_Simulation_Directory")
        return false
    end

end

"""
    Native_Count_Excel(directory, Number_Timesteps)

Creates `Appending_Native_Count.csv` if it does not already exist.

# Arguments
- `directory`: Path to the root simulation directory.
- `Number_Timesteps`: Number of simulation steps to include in the file.

# Calls
- `Checks_Oligomer_Excel_Present()`
- `Append_Number_Timesteps()`

# Notes
- The CSV is initialized with a single `Timesteps` column.
"""


function Native_Count_Excel(directory, Number_Timesteps) #NATIVE
    if Checks_Oligomer_Excel_Present(directory) == false
        Compare_Simulation_Directory = directory * "/Compare_Simulations"
        File_Name = "Appending_Native_Count.csv"
        Complete_File_Path = joinpath(Compare_Simulation_Directory, File_Name)
        Timesteps = Append_Number_Timesteps(Number_Timesteps)
        CSV.write(Complete_File_Path, Timesteps)
    end
end

"""
    Checks_Oligomer_Excel_Present(directory) -> Bool

Checks whether `Appending_Native_Count.csv` exists in `Compare_Simulations/`.

# Arguments
- `directory`: Root path to simulation data.

# Returns
- `false` if the file is missing; otherwise returns `nothing`.

# Notes
- Used to avoid overwriting existing summary files.
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
    Append_Number_Timesteps(Number_Timesteps) -> DataFrame

Generates a `DataFrame` with a single `Timesteps` column ranging from 0 to `Number_Timesteps`.

# Arguments
- `Number_Timesteps`: The total number of simulation steps.

# Returns
- A DataFrame like: Timesteps = [0, 1, 2, ..., Number_Timesteps]
"""

function Append_Number_Timesteps(Number_Timesteps) 
    return timesteps = DataFrame(Timesteps = 0:Number_Timesteps)
end

"""
    Location_Amyloid_CSV_File(directory) -> String

Returns the full path to `Appending_Amyloid_Count.csv` within the `Compare_Simulations/` folder.

# Arguments
- `directory`: Root directory where simulation results are stored.
"""


function Location_Amyloid_CSV_File(directory) 
    return joinpath(directory, "Compare_Simulations", "Appending_Amyloid_Count.csv")
end

"""
    Location_Native_CSV_File(directory) -> String

Returns the full path to `Appending_Native_Count.csv` in `Compare_Simulations/`.

# Arguments
- `directory`: Root directory where simulation results are stored.
"""

function Location_Native_CSV_File(directory) #NATIVE COUNT
    return joinpath(directory, "Compare_Simulations", "Appending_Native_Count.csv")
end

"""
    Retrieve_Simulation_Name(Path) -> String

Extracts the simulation folder name from a full path using a regular expression.

# Arguments
- `Path`: Full file path to a simulation folder.

# Returns
- A string like `"Simulation_2025-05-18_14-00-00"`, used for naming output columns.
"""


function Retrieve_Simulation_Name(Path)
    Simulation_Name = match(r"(Simulation.*)$", Path)
    return Simulation_Name.match
end