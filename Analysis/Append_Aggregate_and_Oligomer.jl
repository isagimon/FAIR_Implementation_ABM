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
`Simulation_Results.csv` files and extracting relevant columns for comparison.

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

Main driver function that processes aggregate count data across all simulation folders.

# Arguments
- `directory::String`: Path to the top-level directory containing simulation data.
- `Number_Timesteps::Int`: Total number of simulation timesteps.

# Behavior
- Initializes data structures.
- Extracts and appends aggregate data to summary CSVs.

# Calls
- `Extract_Aggregate_Count_Excel_Sheets`
"""
function run_process_aggregate_excel_sheets(directory::String, Number_Timesteps::Int)
    # Initialize arrays
    Avg_Oligomer_Count = zeros(Int, Number_Timesteps, 2) #TIMESTEP HAS TO BEGIN AT ZERO 
    Avg_Fibril_Count = zeros(Int, Number_Timesteps, 2) #TIMESTEP HAS TO BEGIN AT ZERO 

    timestamp = string(now(), dateformat"mm-dd-yyyy HH:MM:SS.sss")

    """
        Extract_Aggregate_Count_Excel_Sheets()

    Collects and appends aggregate data (amorphous and oligomer) across all simulation folders
    in the specified directory.

    # Behavior
    - Generates CSVs if not already present.
    - Iterates over folders starting with "Simulation".
    - Appends relevant columns to central data files.
    """
    function Extract_Aggregate_Count_Excel_Sheets()
        Amorphous_Count_Excel(directory)
        Oligomer_Count_Excel(directory)
        All_Folders = readdir(directory, join=true)
        Simulation_Folders = filter(folder -> isdir(folder) && startswith(basename(folder), "Simulation"), All_Folders)

        for Folder in Simulation_Folders
            Simulation_Results_File = joinpath(Folder, "Simulation_Results.csv")
            println("This is Folder: $Folder")
            Save_Amorphous_Count_Column(Simulation_Results_File, Folder)
            Save_Oligomer_Count_Column(Simulation_Results_File, Folder)
        end
    end

    Extract_Aggregate_Count_Excel_Sheets()
end

"""
    Save_Amorphous_Count_Column(Amorphous_Count_File, Full_Path)

Reads and isolates the third column (aggregate count) from a CSV file and appends it
to the summary fibril CSV.

# Arguments
- `Amorphous_Count_File`: Path to the Simulation_Results.csv file.
- `Full_Path`: Full path to the specific simulation folder.

# Calls
- `Append_Amorphous_Column`
"""
function Save_Amorphous_Count_Column(Amorphous_Count_File, Full_Path)
    Read_File = CSV.read(Amorphous_Count_File, DataFrame)
    Aggregate_Column_Name = names(Read_File)[3]
    Aggregate_Column_Data = Read_File[!, Aggregate_Column_Name]
    Aggregate_Column_DataFrame = DataFrame(Aggregate_Column = Aggregate_Column_Data)
    Append_Amorphous_Column(Aggregate_Column_DataFrame, Full_Path)
end

"""
    Append_Amorphous_Column(Aggregate_Column, Full_Path)

Appends the provided aggregate count data column to the centralized fibril CSV file.

# Arguments
- `Aggregate_Column`: DataFrame with one column containing aggregate counts.
- `Full_Path`: Full path to the simulation folder used to generate the column label.

# Global dependencies
- `directory` (used in CSV file path construction)
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
    Save_Oligomer_Count_Column(Oligomer_Count_File, Full_Path)

Reads and isolates the second column (oligomer count) from a CSV file and appends it
to the summary oligomer CSV.

# Arguments
- `Oligomer_Count_File`: Path to the Simulation_Results.csv file.
- `Full_Path`: Full path to the simulation folder.
"""
function Save_Oligomer_Count_Column(Oligomer_Count_File, Full_Path)
    Read_File = CSV.read(Oligomer_Count_File, DataFrame)
    Oligomer_Column_Name = names(Read_File)[2]
    Oligomer_Column_Data = Read_File[!, Oligomer_Column_Name]
    Oligomer_Column_DataFrame = DataFrame(Oligomer_Column = Oligomer_Column_Data)
    Append_Oligomer_Column(Oligomer_Column_Data, Full_Path)
end

"""
    Append_Oligomer_Column(Oligomer_Column, Full_Path)

Appends the provided oligomer count data column to the centralized oligomer CSV file.

# Arguments
- `Oligomer_Column`: Array of oligomer count values.
- `Full_Path`: Full path to the simulation folder used to generate the column label.
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
    Amorphous_Count_Excel(directory)

Creates a new fibril count CSV file if it does not already exist.

# Arguments
- `directory`: Path to the directory containing Compare_Simulations folder.

# Calls
- `Append_Number_Timesteps`
"""
function Amorphous_Count_Excel(directory)
    if !Checks_Amorphous_Excel_Present(directory)
        Compare_Simulation_Directory = directory * "/Compare_Simulations"
        File_Name = "Appending_Fibril_Count.csv"
        Complete_File_Path = joinpath(Compare_Simulation_Directory, File_Name)
        Timesteps = Append_Number_Timesteps(Number_Timesteps)
        CSV.write(Complete_File_Path, Timesteps)
    end
end

"""
    Oligomer_Count_Excel(directory)

Creates a new oligomer count CSV file if it does not already exist.

# Arguments
- `directory`: Path to the directory containing Compare_Simulations folder.

# Calls
- `Append_Number_Timesteps`
"""
function Oligomer_Count_Excel(directory)
    if !Checks_Oligomer_Excel_Present(directory)
        Compare_Simulation_Directory = directory * "/Compare_Simulations"
        File_Name = "Appending_Oligomer_Count.csv"
        Complete_File_Path = joinpath(Compare_Simulation_Directory, File_Name)
        Timesteps = Append_Number_Timesteps(Number_Timesteps)
        CSV.write(Complete_File_Path, Timesteps)
    end
end

"""
    Checks_Amorphous_Excel_Present(directory)

Returns `true` if the fibril count CSV already exists, `false` otherwise.
Also prints a message if the file is not found.

# Arguments
- `directory`: Path to the directory to check.
"""
function Checks_Amorphous_Excel_Present(directory)
    Compare_Simulation_Directory = directory * "/Compare_Simulations"
    CSV_File = "Appending_Fibril_Count.csv"
    File_Path = joinpath(Compare_Simulation_Directory, CSV_File)
    if !isfile(File_Path)
        println("File is not found: $Compare_Simulation_Directory")
        return false
    end
end

"""
    Checks_Oligomer_Excel_Present(directory)

Returns `true` if the oligomer count CSV already exists, `false` otherwise.
Also prints a message if the file is not found.

# Arguments
- `directory`: Path to the directory to check.
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

Creates a DataFrame with a single column "Timesteps" ranging from 1 to `Number_Timesteps`.

# Arguments
- `Number_Timesteps`: Total number of timesteps to include.

# Returns
- DataFrame with "Timesteps" column.
"""
function Append_Number_Timesteps(Number_Timesteps)
    return DataFrame(Timesteps = 1:Number_Timesteps)
end

"""
    Location_Fibril_CSV_File(directory::String) -> String

Constructs the file path for the fibril count CSV.

# Arguments
- `directory`: Path to the base directory.

# Returns
- Full path to "Appending_Fibril_Count.csv"
"""
function Location_Fibril_CSV_File(directory)
    return joinpath(directory, "Compare_Simulations", "Appending_Fibril_Count.csv")
end

"""
    Location_Oligomer_CSV_File(directory::String) -> String

Constructs the file path for the oligomer count CSV.

# Arguments
- `directory`: Path to the base directory.

# Returns
- Full path to "Appending_Oligomer_Count.csv"
"""
function Location_Oligomer_CSV_File(directory)
    return joinpath(directory, "Compare_Simulations", "Appending_Oligomer_Count.csv")
end

"""
    Retrieve_Simulation_Name(Path::String) -> String

Extracts the simulation name (e.g., "Simulation1") from a file path.

# Arguments
- `Path`: File path from which to extract the simulation folder name.

# Returns
- The simulation folder name as a string.
"""
function Retrieve_Simulation_Name(Path)
    Simulation_Name = match(r"(Simulation.*)$", Path)
    return Simulation_Name.match
end
