using Random
using DataFrames
using CSV
using Dates
using Base.Threads

"""
Module: Environment_and_Movement

Performs the simulation of protein aggregation dynamics using a 3D lattice.
Includes monomer diffusion, conformational transitions (Native ↔ aggregateProne), and aggregation (Oligomer, Fibril).
Implements FAIR principles:

- **Findable**: Structured naming, exported per-timestep state CSVs.
- **Accessible**: Results stored in interoperable formats (CSV).
- **Interoperable**: Uses Julia-native types and standard packages.
- **Reusable**: Modular code, consistent global state, reproducible runs.

Authors: Santiago Schnell, Conner Sandefur, Isabella Gimón
Creation Date: [Date]
Last Modified: [Date]
Dependencies: Random, Plots, DataFrames, CSV, Dates, Profile, Base.Threads, Agents.jl
License: http://www.apache.org/licenses/LICENSE-2.0
"""


include(joinpath(@__DIR__, "Agents.jl"))


##########################
# INPUT PARAMETERS
##########################

##########################
# INPUT PARAMETERS (initialized later)
##########################

# These are set by `apply_parameters!()` after Parameters is loaded from CSV.
global MAX_NumberMovements::Int = 0
global Native_to_AggregateProne::Float64 = 0.0
global AggregateProne_to_Native::Float64 = 0.0
global Oligomer_Formation::Float64 = 0.0
global Oligomer_Dissociation_rate::Float64 = 0.0
global Fibril_Formation::Float64 = 0.0
global Fibril_Growth::Float64 = 0.0
global Probability_of_Oligomer_Removal::Float64 = 0.0

global OUTPUT_ROOT::String = ""
global Directory::String = ""  # legacy alias

# Derived parameter (depends on counts loaded in Agents.jl)
global max_fibril_size::Int = 0

"""
apply_parameters!()

Reads values from the global `Parameters` dictionary and sets derived globals used by the simulation.
Call this *after* `Parameters = load_csv_parameters(PARAMETER_FILE)` is executed.
"""
function apply_parameters!()
    # --- Simulation controls ---
    global MAX_NumberMovements        = Int(Parameters["MAX_NumberMovements"])
    global Native_to_AggregateProne   = Float64(Parameters["Native_to_AggregateProne"])
    global AggregateProne_to_Native   = Float64(Parameters["AggregateProne_to_Native"])
    global Oligomer_Formation         = Float64(Parameters["Oligomer_Formation"])
    global Oligomer_Dissociation_rate = Float64(Parameters["Oligomer_Dissociation_rate"])
    global Fibril_Formation           = Float64(Parameters["Fibril_Formation"])
    global Fibril_Growth              = Float64(Parameters["Fibril_Growth"])
    global Probability_of_Oligomer_Removal = Float64(Parameters["Probability_of_Oligomer_Removal"])

    global Max_NumberMonomers_Native         = Int(Parameters["Max_NumberMonomers_Native"])
    global Max_NumberMonomers_AggregateProne = Int(Parameters["Max_NumberMonomers_AggregateProne"])

   # --- Output root default (can later be overridden by ENV or run_simulation output_root) ---
local data_root = joinpath(@__DIR__, "..", "Data_Collection")
mkpath(data_root)

if haskey(Parameters, "Directory") && !isempty(strip(String(Parameters["Directory"])))
    global OUTPUT_ROOT = abspath(String(Parameters["Directory"]))
else
    global OUTPUT_ROOT = abspath(data_root)
end

global Directory = OUTPUT_ROOT

    # --- Derived parameter ---
    global max_fibril_size = Int(Max_NumberMonomers_AggregateProne + Max_NumberMonomers_Native)

    return nothing
end






# Run-specific output directory and identifier (set when Make_Directory() is called)
global timestamp = ""
global directory = ""



##########################
# GLOBAL STATE VARIABLES
##########################

global timesteps = 0            # Current timestep counter
global CurrentTimeStep = 0      # Counter for exporting timestep information
global Fibril_Length_Count      # DataFrame to store fibril length counts
const Total_Cleared_Monomers = Ref(0) #Variable keeps track of monomers cleared

# AUTHOR NOTE:
# "None" is included as an explicit stay-in-place option in Possible_Movement_Options.
# Because movement is sampled uniformly, this gives each monomer a 1/19 probability
# of remaining in place during a timestep before any movement-function dispatch occurs.
#
# Keep "None" if an explicit no-move event is biologically intended.
# Remove "None" if monomers should always attempt one of the 18 geometric moves,
# and only stay in place when movement is blocked or fails by rule.
# Movement directions
const Possible_Movement_Options = [
    "One", "Two", "Three", "Four", "Five", "Six", "Seven", "Eight", "Nine", "Ten",
    "Eleven", "Twelve", "Thirteen", "Fourteen", "Fifteen", "Sixteen", "Seventeen",
    "Eighteen", "None"
]                               # Array of possible movement options

# Movement and tracking
global Possible_Coordinate_Movements_Dict = Dict{Tuple{Float64, Float64, Float64}, Tuple{Float64, Float64, Float64}}() # Dictionary to store possible coordinate movements
global available_numbers = collect(1:200000)    # Array of available unique numbers

##########################
# DATA COLLECTION STRUCTURES
##########################

results_df       = DataFrame(Timestep = Int[], Oligomers = Int[], Aggregates = Int[])                # DataFrame to store simulation results
results_df_two   = DataFrame(Timestep = Int[], Native = Int[], AggregateProne = Int[])                      # DataFrame to store native and aggregateProne counts
msd_data         = DataFrame(Timestep = Int[], MSD_Monomer = Float64[], MSD_Aggregate = Float64[])   # DataFrame to store MSD data
Number_Monomers_Cleared = DataFrame(Timestep = Int[], Number_Monomers_Cleared = Int[]) #DataFrame to store cleared monomers 

##########################
# THREAD SAFETY
##########################

const dict_lock = ReentrantLock()

##########################
# TIMESTAMPING
##########################

"""
safe_timestamp()

Returns a filesystem-safe timestamp string (safe on Windows/macOS/Linux).
"""
function safe_timestamp()
    return Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
end

"""
    Resolve_Output_Directory(output_root, repo_root, data_root) -> String

Determines the base output directory for a simulation run.

Priority order:
1. `ENV["FAIR_ABM_OUTPUT_DIR"]` if defined and nonempty
2. Explicit `output_root` passed to `Make_Directory`
3. `Initial_Conditions["Directory"]` if available
4. Fallback to `data_root`

If a candidate path resolves to the repository root, it is redirected
to `data_root` to avoid writing run folders into the repo root.

# Arguments
- `output_root`: Optional user-specified output directory.
- `repo_root`: Absolute path to the repository root.
- `data_root`: Absolute path to the `Data_Collection` directory.

# Returns
- `String`: Absolute path to the chosen output directory.
"""

function Resolve_Output_Directory(output_root, repo_root, data_root)
    # 1. Environment variable wins
    env_out = get(ENV, "FAIR_ABM_OUTPUT_DIR", "")
    if !isempty(env_out)
        return abspath(env_out)
    end

    # 2. Explicit function argument wins if ENV is not set
    if output_root !== nothing && !isempty(strip(String(output_root)))
        return abspath(String(output_root))
    end

    # 3. Use Directory from the loaded parameter CSV if present
    if @isdefined(Parameters) && haskey(Parameters, "Directory")
        param_dir = String(Parameters["Directory"])
        if !isempty(strip(param_dir))
            return abspath(param_dir)
        end
    end

    # 4. Fallback to repo/Data_Collection
    return abspath(data_root)
end
"""
Make_Directory()

Creates the directory structure for saving simulation results and calls Input_Parameters().

# Global variables modified
- `directory`

# Calls
- `Input_Parameters`
"""
function Make_Directory(; output_root::Union{Nothing,AbstractString}=nothing,
                        run_id::AbstractString=safe_timestamp(),
                        parameter_file::Union{Nothing,AbstractString}=nothing)

    local repo_root = abspath(joinpath(@__DIR__, ".."))
    local data_root = abspath(joinpath(repo_root, "Data_Collection"))
    mkpath(data_root)

    output_root = Resolve_Output_Directory(output_root, repo_root, data_root)
    mkpath(output_root)

    global timestamp = run_id
    global directory = abspath(joinpath(output_root, "Simulation_$timestamp"))
    mkpath(directory)

    if parameter_file !== nothing
        param_copy_path = joinpath(directory, "Input_Parameters_used.csv")
        cp(parameter_file, param_copy_path; force=true)
    end

    Input_Parameters()

    return directory
end

"""
Input_Parameters()

Saves the input parameters of the simulation to a CSV file.

# Global variables read
- `timestamp`
- `Lattice_Size`
- `Max_NumberMonomers_Native`
- `Max_NumberMonomers_AggregateProne`
- `MAX_NumberMovements`
- `Native_to_AggregateProne`
- `AggregateProne_to_Native`
- `Oligomer_Formation`
- `Oligomer_Dissociation_rate`
- `Fibril_Formation`
- `Fibril_Growth`

# Calls
- `Append_Input_Parameters`
"""



function Input_Parameters()
    # Store a run-local copy of key parameters for provenance
    File_Path = joinpath(directory, "Simulation_Information.csv")

    Data = DataFrame(
        Run_ID = [timestamp],
        Lattice_Size = [Lattice_Size],
        Number_Native_Monomers = [Max_NumberMonomers_Native],
        Number_AggregateProne_Monomers = [Max_NumberMonomers_AggregateProne],
        Timesteps = [MAX_NumberMovements],
        P_Native_to_AggregateProne = [Native_to_AggregateProne],
        P_AggregateProne_to_Native = [AggregateProne_to_Native],
        P_Oligomer_Formation = [Oligomer_Formation],
        P_Oligomer_Dissociation = [Oligomer_Dissociation_rate],
        P_Fibril_Formation = [Fibril_Formation],
        P_Fibril_Growth = [Fibril_Growth],
        P_Oligomer_Removal = [Probability_of_Oligomer_Removal],
        Crowders_Enabled = [Obstacle],
        Obstacle_Radius = [Obstacle_Radius],
        Crowder_Concentration_Spheres = [Crowder_Concentration_Spheres],
    )

    Append_Input_Parameters(File_Path, Data)
end

"""
Append_Input_Parameters(File_Path, Data)

Appends the input parameters to an existing CSV file, creating the file if it does not exist.

# Arguments
- `File_Path::String`: The path to the CSV file.
- `Data::DataFrame`: The DataFrame containing the input parameters.
"""
function Append_Input_Parameters(File_Path, Data)
    if isfile(File_Path)
        existing_data = CSV.read(File_Path, DataFrame)
        appended_data = vcat(existing_data, Data)
        CSV.write(File_Path, appended_data)
    else
        CSV.write(File_Path, Data)
    end
end

"""
Counting_Timesteps()

Increments the global timestep counter and prints the current timestep.

# Global variables modified
- `timesteps`

# Returns
- `Int`: The incremented timestep count.
"""

function Counting_Timesteps()
    global timesteps = timesteps + 1
    println("----------------------------------")
    println("This is what timesteps looks like: ",timesteps)
    return timesteps
end

"""
Randomly_Chooses_Monomer()

Randomly selects a monomer from the lattice.

# Global variables read
- `Locations_and_States_Dict`

# Returns
- `Tuple{Float64, Float64, Float64}`: The coordinates of the randomly chosen monomer.
"""

function Randomly_Chooses_Monomer()

    while true
        # Choose a random key (coordinate) from the dictionary
        RandomLocation = rand(keys(Locations_and_States_Dict))
        State, _ = Locations_and_States_Dict[RandomLocation]
        
        # Check if the state is non-zero (indicating it's a monomer)
        if State != 0
            return RandomLocation
        end
    end
end

"""
Current_Time()

Generates a timestamp string.

# Returns
- `String`: The timestamp string.
"""

function Current_Time()

current_time = now()
hour_now = hour(current_time)
minute_now = minute(current_time)
second_now = second(current_time)
millisecond_now = Dates.millisecond(current_time)
am_pm = ifelse(hour_now < 12, "AM", "PM")
adjusted_hour = ifelse(hour_now > 12, hour_now - 12, ifelse(hour_now == 0, 12, hour_now))

timestamp_for_timestep = string(month(current_time), ":", day(current_time), ":", year(current_time), " Time ", adjusted_hour, "-", minute_now, "-", second_now, ".", millisecond_now, " ", am_pm)
return timestamp_for_timestep

end

"""
Initial_Conditions()

Initializes DataFrames to store simulation results.

# Global variables modified
- `results_df_two`
- `results_df`
- 'Number_Monomers_Cleared'
"""

function Initial_Conditions()
    global results_df_two
    global results_df
    global Number_Monomers_Cleared
    push!(results_df_two, (0,Max_NumberMonomers_Native, Max_NumberMonomers_AggregateProne ))
    push!(results_df, (0, 0, 0))
    push!(Number_Monomers_Cleared, (0, 0))
end

"""
Save_MSD_Data(timesteps)

Saves the Mean Squared Displacement (MSD) data for monomers and aggregates.

# Arguments
- `timesteps::Int`: The current timestep.

# Global variables modified
- `msd_data`

# Calls
- `compute_MSD`
- `compute_MSD_aggregates`
"""

function Save_MSD_Data(timesteps)
    global msd_data
    msd_calculated_monomers = compute_MSD() #This only meant for native or AggregateProne monomers
    msd_calculated_aggregates = compute_MSD_aggregates()
    push!(msd_data, (timesteps, msd_calculated_monomers, msd_calculated_aggregates))
end

"""
Save_Number_Monomers_Cleared(timesteps)

Appends a record of the current `timesteps` and the total number of cleared monomers to `Number_Monomers_Cleared`.

# Global variables read
- `Number_Monomers_Cleared`
- `Total_Cleared_Monomers`
"""

function Save_Number_Monomers_Cleared(timesteps)
    global Number_Monomers_Cleared
    push!(Number_Monomers_Cleared, (timesteps, Total_Cleared_Monomers[]))
end

"""
Export_MSD_Data()

Exports the Mean Squared Displacement (MSD) data to a CSV file.

# Global variables read
- `directory`
- `msd_data`
"""

function Export_MSD_Data()
    file_path = joinpath(directory, "MSD_Data.csv")
    CSV.write(file_path, msd_data)
end


"""
Randomly_Chooses_Movement()

Randomly selects a movement option for a monomer.

# Returns
- `String`: The randomly chosen movement option.
"""

function Randomly_Chooses_Movement() 
    return rand(Possible_Movement_Options)
end

"""
Export_Timestep_Information()

Exports the state of the lattice at the current timestep to a CSV file.

# Global variables read
- `Locations_and_States_Dict`

# Calls
- `Current_TimeStep`
- `Export_DataFrame`
"""

function Export_Timestep_Information()
    #global Locations_and_States_Dict
    df = DataFrame()

    # Extract values from each key-value pair in Locations_and_States_Dict
    df.X_Coordinates = [coord[1] for coord in keys(Locations_and_States_Dict)]
    df.Y_Coordinates = [coord[2] for coord in keys(Locations_and_States_Dict)]
    df.Z_Coordinates = [coord[3] for coord in keys(Locations_and_States_Dict)]
    df.States = [value[1] for value in values(Locations_and_States_Dict)]
    df.Unique_Number = [value[2] for value in values(Locations_and_States_Dict)]

    Timestep = Current_TimeStep()
    Export_DataFrame(df, Timestep)
end

"""
Current_TimeStep()

Increments the global CurrentTimeStep counter.

# Global variables modified
- `CurrentTimeStep`

# Returns
- `Int`: The incremented CurrentTimeStep count.
"""

function Current_TimeStep()
    global CurrentTimeStep = CurrentTimeStep + 1
    return CurrentTimeStep
end

"""
Export_DataFrame(df, CurrentTimestep)

Exports a DataFrame to a CSV file, with filename based on the timestep.

# Arguments
- `df::DataFrame`: The DataFrame to export.
- `CurrentTimestep::Int`: The current timestep.

# Global variables read
- `directory`
"""


function Export_DataFrame(df, CurrentTimestep)
    # Mirror the original naming convention (Timestep00X / Timestep0XX / TimestepXXX)
    file_name = if CurrentTimestep < 10
        "Timestep00$(CurrentTimestep).csv"
    elseif CurrentTimestep < 100
        "Timestep0$(CurrentTimestep).csv"
    else
        "Timestep$(CurrentTimestep).csv"
    end

    CSV.write(joinpath(directory, file_name), df)
end

"""
Count_Oligomers_Aggregates()

Counts the number of oligomers and aggregates in the lattice.

# Global variables read
- `Locations_and_States_Dict`

# Returns
- `Tuple{Int, Int}`: The count of oligomers and aggregates.
"""

function Count_Oligomers_Aggregates()
    oligomers_count = 0
    aggregates_count = 0

    # Loop through the dictionary to count monomers in state 3 and 4
    for (_, (state, _)) in Locations_and_States_Dict
        if state == 3
            oligomers_count += 1
        elseif state == 4
            aggregates_count += 1
        end
    end

    return oligomers_count, aggregates_count
end

"""
Count_Native_AggregateProne()

Counts the number of native and AggregateProne monomers in the lattice.

# Global variables read
- `Locations_and_States_Dict`

# Returns
- `Tuple{Int, Int}`: The count of native and AggregateProne monomers.
"""

function Count_Native_AggregateProne()
    native_count = 0
    AggregateProne_count = 0

    # Loop through the dictionary to count monomers in state 3 and 4
    for (_, (state, _)) in Locations_and_States_Dict
        if state == 1
            native_count += 1
        elseif state == 2
            AggregateProne_count += 1
        end
    end

    return native_count, AggregateProne_count

end

"""
Save_Data(timestep, oligomers, aggregates)

Saves the counts of oligomers and aggregates for a given timestep.

# Arguments
- `timestep::Int`: The current timestep.
- `oligomers::Int`: The count of oligomers.
- `aggregates::Int`: The count of aggregates.

# Global variables modified
- `results_df`
"""


function Save_Data(timestep, oligomers, aggregates)
    # Append the results to the DataFrame in memory
    global results_df
    push!(results_df, (timestep, oligomers, aggregates))
end

"""
Save_Data_Two(timesteps, native, AggregateProne)

Saves the counts of native and AggregateProne monomers for a given timestep.

# Arguments
- `timesteps::Int`: The current timestep.
- `native::Int`: The count of native monomers.
- `AggregateProne::Int`: The count of AggregateProne monomers.

# Global variables modified
- `results_df_two`
"""

function Save_Data_Two(timesteps, native, AggregateProne)
    global results_df_two
    push!(results_df_two, (timesteps, native, AggregateProne))
end

"""
Create_Fibril_Length_DataFrame()

Creates a DataFrame to store the counts of fibrils of different lengths.

# Global variables modified
- `Fibril_Length_Count`

# Calls
- `Format_Fibril_Length_DataFrame`
"""

function Create_Fibril_Length_DataFrame()

    # 🔴 SAFETY CHECK
    if max_fibril_size == 0
        error("max_fibril_size == 0. apply_parameters!() was not executed correctly before creating the fibril dataframe.")
    end

    # Get the data and column names from the formatter
    data, column_names = Format_Fibril_Length_DataFrame()

    # Initialize the global Fibril_Length_Count DataFrame
    global Fibril_Length_Count = DataFrame(data, Symbol.(column_names))
end


"""
Format_Fibril_Length_DataFrame()

Formats the data for the fibril length count DataFrame.

# Global variables read
- `MAX_NumberMovements`
- `max_fibril_size`

# Returns
- `Tuple{Matrix{Int64}, Vector{String}}`: The data matrix and column names.
"""

function Format_Fibril_Length_DataFrame()
    # Ensure MAX_NumberMovements is an integer
    max_movements = Int(MAX_NumberMovements)

    # Create the Timesteps column (row numbers)
    timesteps = collect(1:max_movements)

    # Create columns for fibril sizes (1 to max_fibril_size)
    column_names = ["Timesteps"; string.(1:max_fibril_size)]  # Column names as strings
    data = hcat(timesteps, zeros(Int, max_movements, max_fibril_size))  # Initialize data with zeros

    return data, column_names

end

"""
Export_Final_Results()

Exports the final simulation results to a CSV file.

# Global variables read
- `directory`
- `results_df`
"""

function Export_Final_Results()
    file_path = joinpath(directory, "Oligomer_and_Aggregate_Count_Results.csv")
    CSV.write(file_path, results_df)
end

"""
Export_Final_Results_Two()

Exports the final counts of native and AggregateProne monomers to a CSV file.

# Global variables read
- `directory`
- `results_df_two`
"""

function Export_Final_Results_Two()
    file_path = joinpath(directory, "Native_and_AggregateProne_Count_Results.csv")
    CSV.write(file_path, results_df_two)
end

"""
Export_Fibril_Length_Count()

Exports the fibril length count data to a CSV file.

# Global variables read
- `directory`
- `Fibril_Length_Count`
"""


function Export_Fibril_Length_Count()
    #global Directory = directory * "Simulation_$timestamp"
    file_path = joinpath(directory, "Fibril_Length_Count_Results.csv")
    CSV.write(file_path, Fibril_Length_Count)
end

"""
Export_Monomers_Cleared_Data()

Exports the recorded number of cleared monomers over time to a CSV file.

# Global variables read
- `Directory`
- `timestamp`
- `Number_Monomers_Cleared`
"""

function Export_Monomers_Cleared_Data()
    file_path = joinpath(directory, "Oligomers_Cleared.csv")
    CSV.write(file_path, Number_Monomers_Cleared)
end

"""
Total_Monomers()

Calculates the total number of monomers in the simulation.

# Global variables read
- `Max_NumberMonomers_AggregateProne`
- `Max_NumberMonomers_Native`

# Returns
- `Int`: The total number of monomers.
"""

function Total_Monomers()
    Total_Monomers = Max_NumberMonomers_AggregateProne + Max_NumberMonomers_Native
    return Total_Monomers
end

"""
Count_Fibril_Length(timestep)

Counts the number of fibrils of each length at a given timestep and updates the Fibril_Length_Count DataFrame.

# Arguments
- `timestep::Int`: The current timestep.

# Global variables read
- `Locations_and_States_Dict`
- `max_fibril_size`

# Global variables modified
- `Fibril_Length_Count`

# Calls
- `Filter_Aggregate`
- `Update_Fibril_Length_DataFrame`
"""


function Count_Fibril_Length(timestep)
       # Filter coordinates where state == 4 (fibrils)
       Fibril_Keys = Filter_Aggregate()

       # Dictionary to store counts of fibril sizes
       fibril_counts = Dict{Int, Int}()
   
       # Iterate through filtered keys to count fibril sizes
       for key in Fibril_Keys
           _, unique_number = Locations_and_States_Dict[key]
           fibril_counts[unique_number] = get(fibril_counts, unique_number, 0) + 1
       end
   
       # Tally the sizes and update the DataFrame
       fibril_size_counts = Dict{Int, Int}()
       for (_, count) in fibril_counts
           if count <= max_fibril_size
               fibril_size_counts[count] = get(fibril_size_counts, count, 0) + 1
           end
       end
   
       # Update the DataFrame
       Update_Fibril_Length_DataFrame(timestep, fibril_size_counts, max_fibril_size)

end

"""
Update_Fibril_Length_DataFrame(timestep, fibril_lengths, max_fibril_size)

Updates the Fibril_Length_Count DataFrame with the counts of fibrils of different lengths.

# Arguments
- `timestep::Int`: The current timestep.
- `fibril_lengths::Dict{Int, Int}`: A dictionary mapping fibril length to count.
- `max_fibril_size::Int`: The maximum possible fibril size.

# Global variables modified
- `Fibril_Length_Count`
"""


function Update_Fibril_Length_DataFrame(timestep, fibril_lengths, max_fibril_size)
    global Fibril_Length_Count

    # Iterate over all possible fibril sizes (1 to max_fibril_size)
    for size in 1:max_fibril_size
        column_name = Symbol(string(size))
        # Set the count for this size at the current timestep
        Fibril_Length_Count[timestep, column_name] = get(fibril_lengths, size, 0)
    end
end

"""
Movement()

Runs the main simulation loop for protein aggregation dynamics.  
Iterates over timesteps, moves monomers according to movement rules,  
updates states, collects data, and exports results at the end.

# Global variables read
- `MAX_NumberMovements`
- `Locations_and_States_Dict`
- `MovementFunctions`
- `Probability_of_Oligomer_Removal`

# Functions called
- `Create_Fibril_Length_DataFrame`
- `Initial_Conditions`
- `Counting_Timesteps`
- `Current_Time`
- `Randomly_Chooses_Movement`
- `Retrieve_X_Coordinate`
- `Retrieve_Y_Coordinate`
- `Retrieve_Z_Coordinate`
- `Distinguishing_Monomers`
- `Filter_Oligomer`
- `Filter_Aggregate`
- `Filter_Native`
- `Filter_AggregateProne`
- `Save_Data`
- `Save_Data_Two`
- `Save_MSD_Data`
- `Save_Number_Monomers_Cleared`
- `Count_Fibril_Length`
- `Random_Dice`
- `Randomly_Remove_Oligomer`
- `Export_Final_Results`
- `Export_Final_Results_Two`
- `Export_MSD_Data`
- `Export_Fibril_Length_Count`
- `Export_Monomers_Cleared_Data`
"""



function Movement()
    Create_Fibril_Length_DataFrame()
    Initial_Conditions()
    #Export_Timestep_Information()
    # While loop for timesteps
    while MAX_NumberMovements >= Counting_Timesteps()
        # Record the timestamp
        Current_Time_Timestamp = Current_Time()
        Current_Time_Raw = now()
        println("Timestamp: $Current_Time_Timestamp")

        # Loop through all monomers in Locations_and_States_Dict
        for Monomer in keys(Locations_and_States_Dict)
            # Retrieve the state and unique number
            state, unique_number = Locations_and_States_Dict[Monomer]
            
            # Skip inactive monomers (state == 0)
            if state == 0
                continue
            end
            
            # Randomly choose a movement type (e.g., "One" to "Eighteen")
            Movement = Randomly_Chooses_Movement()

            # Guard for stay-in-place option
            if Movement == "None"
                 continue
            end


            # Retrieve the coordinates of the current monomer
            X_Coordinate_Monomer = Retrieve_X_Coordinate(Monomer)
            Y_Coordinate_Monomer = Retrieve_Y_Coordinate(Monomer)
            Z_Coordinate_Monomer = Retrieve_Z_Coordinate(Monomer)
            
            # Get movement options for the chosen movement type
            movement_options = MovementFunctions[Movement]

            # Check conditions and apply the first valid movement
            for (condition, movement_func) in movement_options
                if condition(X_Coordinate_Monomer, Y_Coordinate_Monomer, Z_Coordinate_Monomer)
                    # Compute the new location
                    Location_Movement = movement_func(Monomer)
                    
                    # Move the monomer to the new location
                    Distinguishing_Monomers(Monomer, Location_Movement, Movement)
                    
                    # Stop checking conditions once a movement is applied
                    break
                end
            end
        end

        #Export_Timestep_Information()

           # Collect data after the timestep
           oligomer_keys = Filter_Oligomer()
           aggregate_keys = Filter_Aggregate()
           native_keys = Filter_Native()
           AggregateProne_keys = Filter_AggregateProne()
           
           oligomers, aggregates = length(oligomer_keys), length(aggregate_keys)
           native, AggregateProne = length(native_keys), length(AggregateProne_keys)

          # Save data for this timestep
          Save_Data(timesteps, oligomers, aggregates)
          Save_Data_Two(timesteps, native, AggregateProne)
          Save_MSD_Data(timesteps)
          Save_Number_Monomers_Cleared(timesteps)
          Count_Fibril_Length(timesteps)
        # Timing for the timestep
        Past_Time_Raw = Current_Time_Raw
        Current_Time_Raw = now()
        Time_Taken_For_This_Timestep = Current_Time_Raw - Past_Time_Raw
        println("Time taken for timestep $timesteps: $Time_Taken_For_This_Timestep")

        if Random_Dice() < Probability_of_Oligomer_Removal
            #println("I am going to randomly remove an Oligomer")
            Randomly_Remove_Oligomer()
        end

        
    end
    
    # Export results after all timesteps
    Export_Final_Results()
    Export_Final_Results_Two()
    Export_MSD_Data()
    Export_Fibril_Length_Count()
    Export_Monomers_Cleared_Data()
end

"""
Randomly_Remove_Oligomer()

Selects a random oligomer from the system and removes all its associated monomers.  
Updates the lattice states, clears the oligomer's center of mass information,  
and tracks the number of monomers removed.

# Global variables read
- `Locations_and_States_Dict`

# Functions called
- `Filter_Oligomer`
- `Update_Locations_States`
- `Remove_Center_of_Mass_Info`
- `Track_Cleared_Monomers`
"""


function Randomly_Remove_Oligomer()
    All_Oligomers = Filter_Oligomer()

    if isempty(All_Oligomers)
        println("No oligomers available to remove.")
        return
    end

    Random_Coordinate = rand(All_Oligomers)
    _, Unique_Number = Locations_and_States_Dict[Random_Coordinate]

    coords_to_remove = Vector{typeof(Random_Coordinate)}()

    for (coord, (_, unique_num)) in Locations_and_States_Dict
        if unique_num == Unique_Number
            push!(coords_to_remove, coord)
        end
    end

    removed_count = 0

    for coord in coords_to_remove
        Update_Locations_States(coord, 0, 0)
        removed_count += 1
    end

    Remove_Center_of_Mass_Info(Unique_Number)
    Track_Cleared_Monomers(removed_count)
end

"""
Track_Cleared_Monomers(n::Int)

Increments the total count of cleared monomers by `n`.

# Global variables read
- `Total_Cleared_Monomers`
"""

function Track_Cleared_Monomers(n::Int)
    Total_Cleared_Monomers[] += n
end

"""
Retrieve_X_Coordinate(Monomer)

Retrieves the X coordinate of a monomer.

# Arguments
- `Monomer::Tuple{Float64, Float64, Float64}`: The coordinates of the monomer.

# Returns
- `Float64`: The X coordinate of the monomer.
"""

function Retrieve_X_Coordinate(Monomer)
    return Monomer[1]  # Return the X coordinate
end

"""
Retrieve_Y_Coordinate(Monomer)

Retrieves the Y coordinate of a monomer.

# Arguments
- `Monomer::Tuple{Float64, Float64, Float64}`: The coordinates of the monomer.

# Returns
- `Float64`: The Y coordinate of the monomer.
"""

function Retrieve_Y_Coordinate(Monomer)
    return Monomer[2]  # Return the Y coordinate
end

"""
Retrieve_Z_Coordinate(Monomer)

Retrieves the Z coordinate of a monomer.

# Arguments
- `Monomer::Tuple{Float64, Float64, Float64}`: The coordinates of the monomer.

# Returns
- `Float64`: The Z coordinate of the monomer.
"""

function Retrieve_Z_Coordinate(Monomer)
    return Monomer[3]  # Return the Z coordinate
end

"""
Movement_One_Coordinate(Monomer)

Calculates the new coordinates after a movement of +1 in the X direction.

# Arguments
- `Monomer::Tuple{Float64, Float64, Float64}`: The current coordinates of the monomer.

# Returns
- `Tuple{Float64, Float64, Float64}`: The new coordinates after the movement.
"""

function Movement_One_Coordinate(Monomer) 
    X_Coordinate = Retrieve_X_Coordinate(Monomer) + 1
    Y_Coordinate = Retrieve_Y_Coordinate(Monomer)
    Z_Coordinate = Retrieve_Z_Coordinate(Monomer) 
    return X_Coordinate , Y_Coordinate , Z_Coordinate
end

"""
Movement_One_Coordinate_Exception(Monomer)

Calculates the new coordinates after a movement of +1 in the X direction, with an exception for the boundary.

# Arguments
- `Monomer::Tuple{Float64, Float64, Float64}`: The current coordinates of the monomer.

# Returns
- `Tuple{Float64, Float64, Float64}`: The new coordinates after the movement.
"""

function Movement_One_Coordinate_Exception(Monomer) 
    X_Coordinate = 0
    Y_Coordinate = Retrieve_Y_Coordinate(Monomer)
    Z_Coordinate = Retrieve_Z_Coordinate(Monomer) 
    return X_Coordinate , Y_Coordinate , Z_Coordinate
end

"""
Movement_One_Coordinate_Exception_Second(Monomer)

Calculates the new coordinates after a movement of +1 in the X direction, with a second exception for the boundary.

# Arguments
- `Monomer::Tuple{Float64, Float64, Float64}`: The current coordinates of the monomer.

# Global variables read
- `Lattice_Size`

# Returns
- `Tuple{Float64, Float64, Float64}`: The new coordinates after the movement.
"""

function Movement_One_Coordinate_Exception_Second(Monomer) 
    X_Coordinate = (Retrieve_X_Coordinate(Monomer) + 1 ) - Lattice_Size
    Y_Coordinate = Retrieve_Y_Coordinate(Monomer)
    Z_Coordinate = Retrieve_Z_Coordinate(Monomer) 
    return X_Coordinate , Y_Coordinate , Z_Coordinate
end

"""
Movement_Two_Coordinate(Monomer)

Calculates the new coordinates after a movement of +1 in the Y direction.

# Arguments
- `Monomer::Tuple{Float64, Float64, Float64}`: The current coordinates of the monomer.

# Returns
- `Tuple{Float64, Float64, Float64}`: The updated coordinates after the movement in the +Y direction.
"""

function Movement_Two_Coordinate(Monomer) 
    X_Coordinate = Retrieve_X_Coordinate(Monomer)
    Y_Coordinate = Retrieve_Y_Coordinate(Monomer) + 1
    Z_Coordinate = Retrieve_Z_Coordinate(Monomer) 
    return X_Coordinate , Y_Coordinate , Z_Coordinate
end

"""
Movement_Two_Coordinate_Exception(Monomer)

Handles boundary wrapping when a monomer moves in the +Y direction and exceeds the lattice boundary.

# Arguments
- `Monomer::Tuple{Float64, Float64, Float64}`: The current coordinates of the monomer.

# Global variables read
- `Lattice_Size`

# Returns
- `Tuple{Float64, Float64, Float64}`: The updated coordinates after applying the boundary condition in Y.
"""

function Movement_Two_Coordinate_Exception(Monomer) 
    X_Coordinate = Retrieve_X_Coordinate(Monomer)
    Y_Coordinate = 0
    Z_Coordinate = Retrieve_Z_Coordinate(Monomer) 
    return X_Coordinate , Y_Coordinate , Z_Coordinate
end

"""
Movement_Two_Coordinate_Exception_Second(Monomer)

Handles movement of +1 in Y with second exception boundary wrapping.

# Arguments
- `Monomer::Tuple{Float64, Float64, Float64}`

# Global variables read
- `Lattice_Size`

# Returns
- `Tuple{Float64, Float64, Float64}`
"""


function Movement_Two_Coordinate_Exception_Second(Monomer) 
    X_Coordinate = Retrieve_X_Coordinate(Monomer)
    Y_Coordinate = (Retrieve_Y_Coordinate(Monomer) + 1 ) - Lattice_Size
    Z_Coordinate = Retrieve_Z_Coordinate(Monomer) 
    return X_Coordinate , Y_Coordinate , Z_Coordinate
end

"""
Movement_Three_Coordinate(Monomer)

Calculates the new coordinates after a movement of +1 in the Z direction.

# Arguments
- `Monomer::Tuple{Float64, Float64, Float64}`: The current coordinates of the monomer.

# Returns
- `Tuple{Float64, Float64, Float64}`: The updated coordinates after the movement in the +Z direction.
"""

function Movement_Three_Coordinate(Monomer) 
    X_Coordinate = Retrieve_X_Coordinate(Monomer)
    Y_Coordinate = Retrieve_Y_Coordinate(Monomer)
    Z_Coordinate = Retrieve_Z_Coordinate(Monomer) + 1
    return X_Coordinate , Y_Coordinate , Z_Coordinate
end

"""
Movement_Three_Coordinate_Exception(Monomer)

Handles boundary wrapping when a monomer moves in the +Z direction and exceeds the lattice boundary.

# Arguments
- `Monomer::Tuple{Float64, Float64, Float64}`: The current coordinates of the monomer.

# Global variables read
- `Lattice_Size`

# Returns
- `Tuple{Float64, Float64, Float64}`: The updated coordinates after applying the boundary condition in Z.
"""

function Movement_Three_Coordinate_Exception(Monomer) 
    X_Coordinate = Retrieve_X_Coordinate(Monomer)
    Y_Coordinate = Retrieve_Y_Coordinate(Monomer)
    Z_Coordinate = 0
    return X_Coordinate , Y_Coordinate , Z_Coordinate
end

"""
Movement_Three_Coordinate_Exception_Second(Monomer)

Handles movement of +1 in Z with second exception boundary wrapping.

# Arguments
- `Monomer::Tuple{Float64, Float64, Float64}`

# Global variables read
- `Lattice_Size`

# Returns
- `Tuple{Float64, Float64, Float64}`
"""


function Movement_Three_Coordinate_Exception_Second(Monomer) 
    X_Coordinate = Retrieve_X_Coordinate(Monomer)
    Y_Coordinate = Retrieve_Y_Coordinate(Monomer)
    Z_Coordinate = (Retrieve_Z_Coordinate(Monomer) + 1 ) - Lattice_Size
    return X_Coordinate , Y_Coordinate , Z_Coordinate
end

"""
Movement_Four_Coordinate(Monomer)

Calculates the new coordinates after a movement of -1 in the X direction.

# Arguments
- `Monomer::Tuple{Float64, Float64, Float64}`: The current coordinates of the monomer.

# Returns
- `Tuple{Float64, Float64, Float64}`: The updated coordinates after the movement in the -X direction.
"""


function Movement_Four_Coordinate(Monomer) 
    X_Coordinate = Retrieve_X_Coordinate(Monomer) - 1
    Y_Coordinate = Retrieve_Y_Coordinate(Monomer)
    Z_Coordinate = Retrieve_Z_Coordinate(Monomer) 
    return X_Coordinate , Y_Coordinate , Z_Coordinate
end

"""
Movement_Four_Coordinate_Exception(Monomer)

Handles boundary wrapping when a monomer moves in the -X direction and crosses the minimum boundary (X = 0).

# Arguments
- `Monomer::Tuple{Float64, Float64, Float64}`: The current coordinates of the monomer.

# Global variables read
- `Lattice_Size`

# Returns
- `Tuple{Float64, Float64, Float64}`: The updated coordinates after applying the boundary condition in X.
"""


function Movement_Four_Coordinate_Exception(Monomer) 
    X_Coordinate = Lattice_Size
    Y_Coordinate = Retrieve_Y_Coordinate(Monomer)
    Z_Coordinate = Retrieve_Z_Coordinate(Monomer) 
    return X_Coordinate , Y_Coordinate , Z_Coordinate
end

"""
Movement_Four_Coordinate_Exception_Second(Monomer)

Handles movement of -1 in X with second exception from 0.5 boundary position.

# Arguments
- `Monomer::Tuple{Float64, Float64, Float64}`

# Global variables read
- `Lattice_Size`

# Returns
- `Tuple{Float64, Float64, Float64}`
"""


function Movement_Four_Coordinate_Exception_Second(Monomer) 
    X_Coordinate = Lattice_Size - .5
    Y_Coordinate = Retrieve_Y_Coordinate(Monomer)
    Z_Coordinate = Retrieve_Z_Coordinate(Monomer) 
    return X_Coordinate , Y_Coordinate , Z_Coordinate
end

"""
Movement_Five_Coordinate(Monomer)

Calculates the new coordinates after a movement of -1 in the Y direction.

# Arguments
- `Monomer::Tuple{Float64, Float64, Float64}`: The current coordinates of the monomer.

# Returns
- `Tuple{Float64, Float64, Float64}`: The updated coordinates after the movement in the -Y direction.
"""


function Movement_Five_Coordinate(Monomer) 
    X_Coordinate = Retrieve_X_Coordinate(Monomer)
    Y_Coordinate = Retrieve_Y_Coordinate(Monomer) - 1
    Z_Coordinate = Retrieve_Z_Coordinate(Monomer) 
    return X_Coordinate , Y_Coordinate , Z_Coordinate
end

"""
Movement_Five_Coordinate_Exception(Monomer)

Handles boundary wrapping when a monomer moves in the -Y direction and crosses the minimum boundary (Y = 0).

# Arguments
- `Monomer::Tuple{Float64, Float64, Float64}`: The current coordinates of the monomer.

# Global variables read
- `Lattice_Size`

# Returns
- `Tuple{Float64, Float64, Float64}`: The updated coordinates after applying the boundary condition in Y.
"""


function Movement_Five_Coordinate_Exception(Monomer) 
    X_Coordinate = Retrieve_X_Coordinate(Monomer)
    Y_Coordinate = Lattice_Size
    Z_Coordinate = Retrieve_Z_Coordinate(Monomer) 
    return X_Coordinate , Y_Coordinate , Z_Coordinate
end

"""
Movement_Five_Coordinate_Exception_Second(Monomer)

Handles movement of -1 in Y with second exception from 0.5 boundary position.

# Arguments
- `Monomer::Tuple{Float64, Float64, Float64}`

# Global variables read
- `Lattice_Size`

# Returns
- `Tuple{Float64, Float64, Float64}`
"""


function Movement_Five_Coordinate_Exception_Second(Monomer) 
    X_Coordinate = Retrieve_X_Coordinate(Monomer)
    Y_Coordinate = Lattice_Size - .5
    Z_Coordinate = Retrieve_Z_Coordinate(Monomer) 
    return X_Coordinate , Y_Coordinate , Z_Coordinate
end

"""
Movement_Six_Coordinate(Monomer)

Calculates the new coordinates after a movement of -1 in the Z direction.

# Arguments
- `Monomer::Tuple{Float64, Float64, Float64}`: The current coordinates of the monomer.

# Returns
- `Tuple{Float64, Float64, Float64}`: The updated coordinates after the movement in the -Z direction.
"""


function Movement_Six_Coordinate(Monomer) 
    X_Coordinate = Retrieve_X_Coordinate(Monomer)
    Y_Coordinate = Retrieve_Y_Coordinate(Monomer)
    Z_Coordinate = Retrieve_Z_Coordinate(Monomer) - 1
    return X_Coordinate , Y_Coordinate , Z_Coordinate
end

"""
Movement_Six_Coordinate_Exception(Monomer)

Handles boundary wrapping when a monomer moves in the -Z direction and crosses the minimum boundary (Z = 0).

# Arguments
- `Monomer::Tuple{Float64, Float64, Float64}`: The current coordinates of the monomer.

# Global variables read
- `Lattice_Size`

# Returns
- `Tuple{Float64, Float64, Float64}`: The updated coordinates after applying the boundary condition in Z.
"""


function Movement_Six_Coordinate_Exception(Monomer) 
    X_Coordinate = Retrieve_X_Coordinate(Monomer)
    Y_Coordinate = Retrieve_Y_Coordinate(Monomer)
    Z_Coordinate = Lattice_Size
    return X_Coordinate , Y_Coordinate , Z_Coordinate
end

"""
Movement_Six_Coordinate_Exception_Second(Monomer)

Handles movement of -1 in Z with second exception from 0.5 boundary position.

# Arguments
- `Monomer::Tuple{Float64, Float64, Float64}`

# Global variables read
- `Lattice_Size`

# Returns
- `Tuple{Float64, Float64, Float64}`
"""


function Movement_Six_Coordinate_Exception_Second(Monomer) 
    X_Coordinate = Retrieve_X_Coordinate(Monomer)
    Y_Coordinate = Retrieve_Y_Coordinate(Monomer)
    Z_Coordinate = Lattice_Size - .5
    return X_Coordinate , Y_Coordinate , Z_Coordinate
end

"""
Movement_Seven_Coordinate(Monomer)

Calculates the new coordinates after a diagonal movement of +0.5 in X and -0.5 in Y.

# Arguments
- `Monomer::Tuple{Float64, Float64, Float64}`: The current coordinates of the monomer.

# Returns
- `Tuple{Float64, Float64, Float64}`: The updated coordinates after the diagonal movement in the XY plane.
"""

function Movement_Seven_Coordinate(Monomer) 
    X_Coordinate = Retrieve_X_Coordinate(Monomer) + .5
    Y_Coordinate = Retrieve_Y_Coordinate(Monomer) - .5
    Z_Coordinate = Retrieve_Z_Coordinate(Monomer)
    return X_Coordinate , Y_Coordinate , Z_Coordinate 
end
 
"""
Movement_Seven_Coordinate_Exception(Monomer)

Handles special boundary cases for diagonal movement in the XY plane when at corners or edges.

# Arguments
- `Monomer::Tuple{Float64, Float64, Float64}`: The current coordinates of the monomer.

# Returns
- `Tuple{Float64, Float64, Float64}`: The coordinates adjusted for edge or corner boundary conditions.
"""


function Movement_Seven_Coordinate_Exception(Monomer) 
    X_Coordinate = Retrieve_X_Coordinate(Monomer)
    Y_Coordinate = Retrieve_Y_Coordinate(Monomer)
    Z_Coordinate = Retrieve_Z_Coordinate(Monomer)
    return X_Coordinate , Y_Coordinate , Z_Coordinate 
end

"""
Movement_Seven_Coordinate_Exception_Second(Monomer)

Handles diagonal XY movement with intermediate boundary exception logic.

# Arguments
- `Monomer::Tuple{Float64, Float64, Float64}`

# Returns
- `Tuple{Float64, Float64, Float64}`
"""


function Movement_Seven_Coordinate_Exception_Second(Monomer) 
    X_Coordinate = Retrieve_Y_Coordinate(Monomer)
    Y_Coordinate = Retrieve_X_Coordinate(Monomer)
    Z_Coordinate = Retrieve_Z_Coordinate(Monomer)
    return  X_Coordinate , Y_Coordinate, Z_Coordinate 
end

"""
Movement_Seven_Coordinate_Exception_Third(Monomer)

Handles edge case for movement on X = Lattice_Size and Y = 0.

# Arguments
- `Monomer::Tuple{Float64, Float64, Float64}`

# Returns
- `Tuple{Float64, Float64, Float64}`
"""

function Movement_Seven_Coordinate_Exception_Third(Monomer) 
    X_Coordinate = 0
    Y_Coordinate = Lattice_Size
    Z_Coordinate = Retrieve_Z_Coordinate(Monomer)
    return X_Coordinate , Y_Coordinate , Z_Coordinate 
end

"""
Movement_Eight_Coordinate(Monomer)

Calculates the new coordinates after a diagonal movement of -0.5 in X and +0.5 in Y.

# Arguments
- `Monomer::Tuple{Float64, Float64, Float64}`: The current coordinates of the monomer.

# Returns
- `Tuple{Float64, Float64, Float64}`: The updated coordinates after the diagonal movement in the XY plane.
"""


function Movement_Eight_Coordinate(Monomer) 
    X_Coordinate = Retrieve_X_Coordinate(Monomer) - .5
    Y_Coordinate = Retrieve_Y_Coordinate(Monomer) + .5
    Z_Coordinate = Retrieve_Z_Coordinate(Monomer) 
    return X_Coordinate , Y_Coordinate , Z_Coordinate
end

"""
Movement_Eight_Coordinate_Exception(Monomer)

Handles special boundary cases for diagonal movement in the XY plane when near corners or edges.

# Arguments
- `Monomer::Tuple{Float64, Float64, Float64}`: The current coordinates of the monomer.

# Returns
- `Tuple{Float64, Float64, Float64}`: The coordinates adjusted for edge or corner boundary conditions.
"""


function Movement_Eight_Coordinate_Exception(Monomer) 
    X_Coordinate = Retrieve_X_Coordinate(Monomer)
    Y_Coordinate = Retrieve_Y_Coordinate(Monomer)
    Z_Coordinate = Retrieve_Z_Coordinate(Monomer) 
    return X_Coordinate , Y_Coordinate , Z_Coordinate
end

"""
Movement_Eight_Coordinate_Exception_Second(Monomer)

Handles flipped XY diagonal transition during movement with partial edge conditions.

# Arguments
- `Monomer::Tuple{Float64, Float64, Float64}`

# Returns
- `Tuple{Float64, Float64, Float64}`
"""

function Movement_Eight_Coordinate_Exception_Second(Monomer) 
    X_Coordinate = Retrieve_Y_Coordinate(Monomer)
    Y_Coordinate = Retrieve_X_Coordinate(Monomer)
    Z_Coordinate = Retrieve_Z_Coordinate(Monomer) 
    return X_Coordinate , Y_Coordinate, Z_Coordinate
end

"""
Movement_Eight_Coordinate_Exception_Third(Monomer)

Wraps from X = 0 and Y = Lattice_Size across boundary diagonally.

# Arguments
- `Monomer::Tuple{Float64, Float64, Float64}`

# Global variables read
- `Lattice_Size`

# Returns
- `Tuple{Float64, Float64, Float64}`
"""

function Movement_Eight_Coordinate_Exception_Third(Monomer) 
    X_Coordinate = Lattice_Size
    Y_Coordinate = 0
    Z_Coordinate = Retrieve_Z_Coordinate(Monomer) 
    return Y_Coordinate , X_Coordinate, Z_Coordinate
end

"""
Movement_Nine_Coordinate(Monomer)

Calculates the new coordinates after a diagonal movement of +0.5 in Y and -0.5 in Z.

# Arguments
- `Monomer::Tuple{Float64, Float64, Float64}`: The current coordinates of the monomer.

# Returns
- `Tuple{Float64, Float64, Float64}`: The updated coordinates after the diagonal movement in the YZ plane.
"""


function Movement_Nine_Coordinate(Monomer) 
    X_Coordinate = Retrieve_X_Coordinate(Monomer)
    Y_Coordinate = Retrieve_Y_Coordinate(Monomer) + .5
    Z_Coordinate = Retrieve_Z_Coordinate(Monomer) - .5
    return X_Coordinate , Y_Coordinate , Z_Coordinate 
end

"""
Movement_Nine_Coordinate_Exception(Monomer)

Handles special boundary cases for diagonal movement in the YZ plane when near corners or edges.

# Arguments
- `Monomer::Tuple{Float64, Float64, Float64}`: The current coordinates of the monomer.

# Returns
- `Tuple{Float64, Float64, Float64}`: The coordinates adjusted for edge or corner boundary conditions.
"""


function Movement_Nine_Coordinate_Exception(Monomer) 
    X_Coordinate = Retrieve_X_Coordinate(Monomer)
    Y_Coordinate = Retrieve_Y_Coordinate(Monomer)
    Z_Coordinate = Retrieve_Z_Coordinate(Monomer)
    return X_Coordinate , Y_Coordinate , Z_Coordinate 
end

"""
Movement_Nine_Coordinate_Exception_Second(Monomer)

Handles YZ diagonal movement boundary exceptions from sides or corners.

# Arguments
- `Monomer::Tuple{Float64, Float64, Float64}`

# Returns
- `Tuple{Float64, Float64, Float64}`
"""

function Movement_Nine_Coordinate_Exception_Second(Monomer) 
    X_Coordinate = Retrieve_X_Coordinate(Monomer)
    Y_Coordinate = Retrieve_Z_Coordinate(Monomer)
    Z_Coordinate = Retrieve_Y_Coordinate(Monomer)
    return X_Coordinate , Y_Coordinate, Z_Coordinate
end

"""
Movement_Nine_Coordinate_Exception_Third(Monomer)

Edge case where Y = Lattice_Size and Z = 0 for diagonal motion.

# Arguments
- `Monomer::Tuple{Float64, Float64, Float64}`

# Returns
- `Tuple{Float64, Float64, Float64}`
"""

function Movement_Nine_Coordinate_Exception_Third(Monomer) 
    X_Coordinate = Retrieve_X_Coordinate(Monomer)
    Y_Coordinate = 0
    Z_Coordinate = Lattice_Size
    return X_Coordinate , Z_Coordinate , Y_Coordinate
end

"""
Movement_Ten_Coordinate(Monomer)

Calculates the new coordinates after a diagonal movement of -0.5 in Y and +0.5 in Z.

# Arguments
- `Monomer::Tuple{Float64, Float64, Float64}`: The current coordinates of the monomer.

# Returns
- `Tuple{Float64, Float64, Float64}`: The updated coordinates after the diagonal movement in the YZ plane.
"""


function Movement_Ten_Coordinate(Monomer) 
    X_Coordinate = Retrieve_X_Coordinate(Monomer)
    Y_Coordinate = Retrieve_Y_Coordinate(Monomer) - .5
    Z_Coordinate = Retrieve_Z_Coordinate(Monomer) + .5
    return X_Coordinate , Y_Coordinate , Z_Coordinate
end

"""
Movement_Ten_Coordinate_Exception(Monomer)

Handles special boundary cases for diagonal movement in the YZ plane when near corners or edges.

# Arguments
- `Monomer::Tuple{Float64, Float64, Float64}`: The current coordinates of the monomer.

# Returns
- `Tuple{Float64, Float64, Float64}`: The coordinates adjusted for edge or corner boundary conditions.
"""


function Movement_Ten_Coordinate_Exception(Monomer) 
    X_Coordinate = Retrieve_X_Coordinate(Monomer)
    Y_Coordinate = Retrieve_Y_Coordinate(Monomer)
    Z_Coordinate = Retrieve_Z_Coordinate(Monomer) 
    return X_Coordinate , Y_Coordinate , Z_Coordinate
end

"""
Movement_Ten_Coordinate_Exception_Second(Monomer)

Handles motion across the ZY plane with boundary condition logic.

# Arguments
- `Monomer::Tuple{Float64, Float64, Float64}`

# Returns
- `Tuple{Float64, Float64, Float64}`
"""

function Movement_Ten_Coordinate_Exception_Second(Monomer) 
    X_Coordinate = Retrieve_X_Coordinate(Monomer)
    Y_Coordinate = Retrieve_Z_Coordinate(Monomer) 
    Z_Coordinate = Retrieve_Y_Coordinate(Monomer) 
    return X_Coordinate , Y_Coordinate , Z_Coordinate
end

"""
Movement_Ten_Coordinate_Exception_Third(Monomer)

Wraps from Z = Lattice_Size and Y = 0 for half-step motion.

# Arguments
- `Monomer::Tuple{Float64, Float64, Float64}`

# Returns
- `Tuple{Float64, Float64, Float64}`
"""

function Movement_Ten_Coordinate_Exception_Third(Monomer) 
    X_Coordinate = Retrieve_X_Coordinate(Monomer)
    Y_Coordinate = Lattice_Size
    Z_Coordinate = 0
    return X_Coordinate , Y_Coordinate , Z_Coordinate
end

"""
Movement_Eleven_Coordinate(Monomer)

Calculates the new coordinates after a diagonal movement of +0.5 in X and -0.5 in Z.

# Arguments
- `Monomer::Tuple{Float64, Float64, Float64}`: The current coordinates of the monomer.

# Returns
- `Tuple{Float64, Float64, Float64}`: The updated coordinates after the diagonal movement in the XZ plane.
"""


function Movement_Eleven_Coordinate(Monomer) 
    X_Coordinate = Retrieve_X_Coordinate(Monomer) + .5
    Y_Coordinate = Retrieve_Y_Coordinate(Monomer) 
    Z_Coordinate = Retrieve_Z_Coordinate(Monomer) - .5
    return X_Coordinate , Y_Coordinate , Z_Coordinate
end

"""
Movement_Eleven_Coordinate_Exception(Monomer)

Handles special boundary cases for diagonal movement in the XZ plane when near corners or edges.

# Arguments
- `Monomer::Tuple{Float64, Float64, Float64}`: The current coordinates of the monomer.

# Returns
- `Tuple{Float64, Float64, Float64}`: The coordinates adjusted for edge or corner boundary conditions.
"""

function Movement_Eleven_Coordinate_Exception(Monomer) 
    X_Coordinate = Retrieve_X_Coordinate(Monomer)
    Y_Coordinate = Retrieve_Y_Coordinate(Monomer)
    Z_Coordinate = Retrieve_Z_Coordinate(Monomer) 
    return X_Coordinate , Y_Coordinate , Z_Coordinate
end

"""
Movement_Eleven_Coordinate_Exception_Second(Monomer)

Handles diagonal XZ movement with Z-based boundary flipping.

# Arguments
- `Monomer::Tuple{Float64, Float64, Float64}`

# Returns
- `Tuple{Float64, Float64, Float64}`
"""

function Movement_Eleven_Coordinate_Exception_Second(Monomer) 
    X_Coordinate = Retrieve_Z_Coordinate(Monomer) 
    Y_Coordinate = Retrieve_Y_Coordinate(Monomer)
    Z_Coordinate = Retrieve_X_Coordinate(Monomer) 
    return X_Coordinate , Y_Coordinate , Z_Coordinate
end

"""
Movement_Eleven_Coordinate_Exception_Third(Monomer)

Edge wrap where Z = 0 and X = Lattice_Size.

# Arguments
- `Monomer::Tuple{Float64, Float64, Float64}`

# Returns
- `Tuple{Float64, Float64, Float64}`
"""

function Movement_Eleven_Coordinate_Exception_Third(Monomer) 
    X_Coordinate = 0
    Y_Coordinate = Retrieve_Y_Coordinate(Monomer)
    Z_Coordinate = Lattice_Size
    return X_Coordinate , Y_Coordinate , Z_Coordinate
end

"""
Movement_Twelve_Coordinate(Monomer)

Calculates the new coordinates after a diagonal movement of -0.5 in X and +0.5 in Z.

# Arguments
- `Monomer::Tuple{Float64, Float64, Float64}`: The current coordinates of the monomer.

# Returns
- `Tuple{Float64, Float64, Float64}`: The updated coordinates after the diagonal movement in the XZ plane.
"""

function Movement_Twelve_Coordinate(Monomer)
    X_Coordinate = Retrieve_X_Coordinate(Monomer) - .5
    Y_Coordinate = Retrieve_Y_Coordinate(Monomer)
    Z_Coordinate = Retrieve_Z_Coordinate(Monomer) + .5
    return X_Coordinate , Y_Coordinate , Z_Coordinate
end

"""
Movement_Twelve_Coordinate_Exception(Monomer)

Handles special boundary cases for diagonal movement in the XZ plane when near corners or edges.

# Arguments
- `Monomer::Tuple{Float64, Float64, Float64}`: The current coordinates of the monomer.

# Returns
- `Tuple{Float64, Float64, Float64}`: The coordinates adjusted for edge or corner boundary conditions.
"""


function Movement_Twelve_Coordinate_Exception(Monomer) 
    X_Coordinate = Retrieve_X_Coordinate(Monomer)
    Y_Coordinate = Retrieve_Y_Coordinate(Monomer)
    Z_Coordinate = Retrieve_Z_Coordinate(Monomer) 
    return X_Coordinate , Y_Coordinate , Z_Coordinate
end

"""
Movement_Twelve_Coordinate_Exception_Second(Monomer)

Boundary case where movement in the XZ plane requires coordinate swap.

# Arguments
- `Monomer::Tuple{Float64, Float64, Float64}`

# Returns
- `Tuple{Float64, Float64, Float64}`
"""

function Movement_Twelve_Coordinate_Exception_Second(Monomer) 
    X_Coordinate = Retrieve_Z_Coordinate(Monomer)
    Y_Coordinate = Retrieve_Y_Coordinate(Monomer)
    Z_Coordinate = Retrieve_X_Coordinate(Monomer)
    return X_Coordinate , Y_Coordinate , Z_Coordinate
end

"""
Movement_Twelve_Coordinate_Exception_Third(Monomer)

Wraps X = 0 and Z = Lattice_Size for half-step motion.

# Arguments
- `Monomer::Tuple{Float64, Float64, Float64}`

# Global variables read
- `Lattice_Size`

# Returns
- `Tuple{Float64, Float64, Float64}`
"""

function Movement_Twelve_Coordinate_Exception_Third(Monomer) 
    X_Coordinate = Lattice_Size
    Y_Coordinate = Retrieve_Y_Coordinate(Monomer)
    Z_Coordinate = 0
    return X_Coordinate , Y_Coordinate , Z_Coordinate
end

"""
Movement_Thirteen_Coordinate_Exception(Monomer)

Handles boundary conditions for diagonal movement involving the XZ plane when at specific corner or edge cases.

# Arguments
- `Monomer::Tuple{Float64, Float64, Float64}`: The current coordinates of the monomer.

# Returns
- `Tuple{Float64, Float64, Float64}`: The coordinates adjusted for boundary conditions during Movement Thirteen.
"""


function Movement_Thirteen_Coordinate_Exception(Monomer) 
    X_Coordinate = Retrieve_X_Coordinate(Monomer)
    Y_Coordinate = Retrieve_Y_Coordinate(Monomer)
    Z_Coordinate = Retrieve_Z_Coordinate(Monomer)
    return X_Coordinate , Y_Coordinate , Z_Coordinate 
end

"""
Movement_Thirteen_Coordinate_Exception_Second(Monomer)

Wraps when X = 0 and Z is internal for diagonal backward motion.

# Returns
- `Tuple{Float64, Float64, Float64}`
"""

function Movement_Thirteen_Coordinate_Exception_Second(Monomer)
    X_Coordinate = Lattice_Size - Retrieve_Z_Coordinate(Monomer)
    Y_Coordinate = Retrieve_Y_Coordinate(Monomer)
    Z_Coordinate = Lattice_Size
    return X_Coordinate , Y_Coordinate , Z_Coordinate 
end

"""
Movement_Thirteen_Coordinate_Exception_Third(Monomer)

Handles XZ plane mirror edge at Z = 0.

# Returns
- `Tuple{Float64, Float64, Float64}`
"""

function Movement_Thirteen_Coordinate_Exception_Third(Monomer) 
    X_Coordinate = Lattice_Size
    Y_Coordinate = Retrieve_Y_Coordinate(Monomer)
    Z_Coordinate = Lattice_Size - Retrieve_X_Coordinate(Monomer)
    return X_Coordinate , Y_Coordinate , Z_Coordinate 
end

"""
Movement_Thirteen_Coordinate_Exception_Fourth(Monomer)

Directly sets coordinate to Lattice_Size on X and Z axis.

# Returns
- `Tuple{Float64, Float64, Float64}`
"""

function Movement_Thirteen_Coordinate_Exception_Fourth(Monomer) 
    X_Coordinate = Lattice_Size
    Y_Coordinate = Retrieve_Y_Coordinate(Monomer)
    Z_Coordinate = Lattice_Size
    return X_Coordinate , Y_Coordinate , Z_Coordinate 
end

"""
Movement_Thirteen_Coordinate(Monomer)

Calculates the new coordinates after a diagonal movement of -0.5 in both X and Z directions.

# Arguments
- `Monomer::Tuple{Float64, Float64, Float64}`: The current coordinates of the monomer.

# Returns
- `Tuple{Float64, Float64, Float64}`: The updated coordinates after the diagonal movement in the XZ plane.
"""


function Movement_Thirteen_Coordinate(Monomer) 
    X_Coordinate = Retrieve_X_Coordinate(Monomer) - .5
    Y_Coordinate = Retrieve_Y_Coordinate(Monomer)
    Z_Coordinate = Retrieve_Z_Coordinate(Monomer) -.5
    return X_Coordinate , Y_Coordinate , Z_Coordinate 
end

"""
Movement_Fourteen_Coordinate_Exception(Monomer)

Handles boundary conditions for diagonal movement in the XZ plane during Movement Fourteen when near corners or edges.

# Arguments
- `Monomer::Tuple{Float64, Float64, Float64}`: The current coordinates of the monomer.

# Returns
- `Tuple{Float64, Float64, Float64}`: The coordinates adjusted for edge or corner boundary conditions.
"""


function Movement_Fourteen_Coordinate_Exception(Monomer)  
    X_Coordinate = Retrieve_X_Coordinate(Monomer)
    Y_Coordinate = Retrieve_Y_Coordinate(Monomer)
    Z_Coordinate = Retrieve_Z_Coordinate(Monomer)
    return X_Coordinate , Y_Coordinate , Z_Coordinate 
end

"""
Movement_Fourteen_Coordinate_Exception_Second(Monomer)

Flips coordinate from Z to X across the boundary plane.

# Returns
- `Tuple{Float64, Float64, Float64}`
"""

function Movement_Fourteen_Coordinate_Exception_Second(Monomer) 
    X_Coordinate = 0
    Y_Coordinate = Retrieve_Y_Coordinate(Monomer)
    Z_Coordinate = Lattice_Size - (Retrieve_X_Coordinate(Monomer))
    return X_Coordinate , Y_Coordinate , Z_Coordinate 
end

"""
Movement_Fourteen_Coordinate_Exception_Third(Monomer)

Wraps across X-Z boundary when X = Lattice_Size and Z is interior.

# Returns
- `Tuple{Float64, Float64, Float64}`
"""

function Movement_Fourteen_Coordinate_Exception_Third(Monomer) 
    X_Coordinate = Lattice_Size - (Retrieve_Z_Coordinate(Monomer))
    Y_Coordinate = Retrieve_Y_Coordinate(Monomer)
    Z_Coordinate = 0
    return X_Coordinate , Y_Coordinate , Z_Coordinate 
end

"""
Movement_Fourteen_Coordinate_Exception_Fourth(Monomer)

Boundary fix for X = Lattice_Size and Z = Lattice_Size

# Returns
- `Tuple{Float64, Float64, Float64}`
"""

function Movement_Fourteen_Coordinate_Exception_Fourth(Monomer) 
    X_Coordinate = Lattice_Size 
    Y_Coordinate = Retrieve_Y_Coordinate(Monomer)
    Z_Coordinate = Lattice_Size
    return X_Coordinate , Y_Coordinate , Z_Coordinate 
end

"""
Movement_Fourteen_Coordinate(Monomer)

Calculates the new coordinates after a diagonal movement of +0.5 in both X and Z directions.

# Arguments
- `Monomer::Tuple{Float64, Float64, Float64}`: The current coordinates of the monomer.

# Returns
- `Tuple{Float64, Float64, Float64}`: The updated coordinates after the diagonal movement in the XZ plane.
"""


function Movement_Fourteen_Coordinate(Monomer) 
    X_Coordinate = Retrieve_X_Coordinate(Monomer) + .5
    Y_Coordinate = Retrieve_Y_Coordinate(Monomer)
    Z_Coordinate = Retrieve_Z_Coordinate(Monomer) + .5
    return X_Coordinate , Y_Coordinate , Z_Coordinate 
end

"""
Movement_Fifteen_Coordinate_Exception(Monomer)

Handles boundary conditions for diagonal movement in the YZ plane during Movement Fifteen when near corners or edges.

# Arguments
- `Monomer::Tuple{Float64, Float64, Float64}`: The current coordinates of the monomer.

# Returns
- `Tuple{Float64, Float64, Float64}`: The coordinates adjusted for edge or corner boundary conditions.
"""


function Movement_Fifteen_Coordinate_Exception(Monomer)  
    X_Coordinate = Retrieve_X_Coordinate(Monomer)
    Y_Coordinate = Retrieve_Y_Coordinate(Monomer)
    Z_Coordinate = Retrieve_Z_Coordinate(Monomer)
    return X_Coordinate , Y_Coordinate , Z_Coordinate 
end

"""
Movement_Fifteen_Coordinate_Exception_Second(Monomer)

Handles diagonal YZ wrapping from Y = 0 and Z internal.

# Returns
- `Tuple{Float64, Float64, Float64}`
"""

function Movement_Fifteen_Coordinate_Exception_Second(Monomer)  
    X_Coordinate = Retrieve_X_Coordinate(Monomer)
    Y_Coordinate = Lattice_Size - Retrieve_Z_Coordinate(Monomer)
    Z_Coordinate = Lattice_Size
    return X_Coordinate , Y_Coordinate , Z_Coordinate 
end

"""
Movement_Fifteen_Coordinate_Exception_Third(Monomer)

Wraps from Z = 0 and internal Y for diagonal transition.

# Returns
- `Tuple{Float64, Float64, Float64}`
"""

function Movement_Fifteen_Coordinate_Exception_Third(Monomer)  
    X_Coordinate = Retrieve_X_Coordinate(Monomer)
    Y_Coordinate = Lattice_Size
    Z_Coordinate = Lattice_Size - Retrieve_Y_Coordinate(Monomer)
    return X_Coordinate , Y_Coordinate , Z_Coordinate 
end

"""
Movement_Fifteen_Coordinate_Exception_Fourth(Monomer)

Edge fix for Y = 0 and Z = 0

# Returns
- `Tuple{Float64, Float64, Float64}`
"""

function Movement_Fifteen_Coordinate_Exception_Fourth(Monomer)
    X_Coordinate = Retrieve_X_Coordinate(Monomer)
    Y_Coordinate = Lattice_Size
    Z_Coordinate = Lattice_Size 
    return X_Coordinate , Y_Coordinate , Z_Coordinate 
end

"""
Movement_Fifteen_Coordinate(Monomer)

Calculates the new coordinates after a diagonal movement of -0.5 in both Y and Z directions.

# Arguments
- `Monomer::Tuple{Float64, Float64, Float64}`: The current coordinates of the monomer.

# Returns
- `Tuple{Float64, Float64, Float64}`: The updated coordinates after the diagonal movement in the YZ plane.
"""

function Movement_Fifteen_Coordinate(Monomer)  
    X_Coordinate = Retrieve_X_Coordinate(Monomer)
    Y_Coordinate = Retrieve_Y_Coordinate(Monomer) - .5 
    Z_Coordinate = Retrieve_Z_Coordinate(Monomer) - .5
    return X_Coordinate , Y_Coordinate , Z_Coordinate 
end

"""
Movement_Sixteen_Coordinate_Exception(Monomer)

Handles boundary conditions for diagonal movement in the YZ plane during Movement Sixteen when near corners or edges.

# Arguments
- `Monomer::Tuple{Float64, Float64, Float64}`: The current coordinates of the monomer.

# Returns
- `Tuple{Float64, Float64, Float64}`: The coordinates adjusted for edge or corner boundary conditions.
"""


function Movement_Sixteen_Coordinate_Exception(Monomer) 
    X_Coordinate = Retrieve_X_Coordinate(Monomer)
    Y_Coordinate = Retrieve_Y_Coordinate(Monomer)
    Z_Coordinate = Retrieve_Z_Coordinate(Monomer)
    return X_Coordinate , Y_Coordinate , Z_Coordinate 
end

"""
Movement_Sixteen_Coordinate_Exception_Second(Monomer)

Wraps to origin from opposing Y and Z limits.

# Returns
- `Tuple{Float64, Float64, Float64}`
"""

function Movement_Sixteen_Coordinate_Exception_Second(Monomer)  
    X_Coordinate = Retrieve_X_Coordinate(Monomer)
    Y_Coordinate = 0
    Z_Coordinate = 0
    return X_Coordinate , Y_Coordinate , Z_Coordinate 
end

"""
Movement_Sixteen_Coordinate_Exception_Third(Monomer)

Reflects Y into Z axis for diagonal wrapping.

# Returns
- `Tuple{Float64, Float64, Float64}`
"""

function Movement_Sixteen_Coordinate_Exception_Third(Monomer)  
    X_Coordinate = Retrieve_X_Coordinate(Monomer)
    Y_Coordinate = 0
    Z_Coordinate = Lattice_Size - Retrieve_Y_Coordinate(Monomer)
    return X_Coordinate , Y_Coordinate , Z_Coordinate 
end

"""
Movement_Sixteen_Coordinate_Exception_Fourth(Monomer)

Reflects Z into Y axis for diagonal wrapping.

# Returns
- `Tuple{Float64, Float64, Float64}`
"""

function Movement_Sixteen_Coordinate_Exception_Fourth(Monomer)  
    X_Coordinate = Retrieve_X_Coordinate(Monomer)
    Y_Coordinate = Lattice_Size - Retrieve_Z_Coordinate(Monomer)
    Z_Coordinate = 0
    return X_Coordinate , Y_Coordinate , Z_Coordinate 
end

"""
Movement_Sixteen_Coordinate(Monomer)

Calculates the new coordinates after a diagonal movement of +0.5 in both Y and Z directions.

# Arguments
- `Monomer::Tuple{Float64, Float64, Float64}`: The current coordinates of the monomer.

# Returns
- `Tuple{Float64, Float64, Float64}`: The updated coordinates after the diagonal movement in the YZ plane.
"""

function Movement_Sixteen_Coordinate(Monomer)  
    X_Coordinate = Retrieve_X_Coordinate(Monomer)
    Y_Coordinate = Retrieve_Y_Coordinate(Monomer) + .5
    Z_Coordinate = Retrieve_Z_Coordinate(Monomer) + .5
    return X_Coordinate , Y_Coordinate , Z_Coordinate 
end

"""
Movement_Seventeen_Coordinate_Exception(Monomer)

Handles boundary conditions for diagonal movement in the XY plane during Movement Seventeen when near corners or edges.

# Arguments
- `Monomer::Tuple{Float64, Float64, Float64}`: The current coordinates of the monomer.

# Returns
- `Tuple{Float64, Float64, Float64}`: The coordinates adjusted for edge or corner boundary conditions.
"""

function Movement_Seventeen_Coordinate_Exception(Monomer) 
    X_Coordinate = Retrieve_X_Coordinate(Monomer)
    Y_Coordinate = Retrieve_Y_Coordinate(Monomer)
    Z_Coordinate = Retrieve_Z_Coordinate(Monomer)
    return X_Coordinate , Y_Coordinate , Z_Coordinate 
end

"""
Movement_Seventeen_Coordinate_Exception_Second(Monomer)

Fixes edge case where both X and Y are at lattice origin or max.

# Returns
- `Tuple{Float64, Float64, Float64}`
"""

function Movement_Seventeen_Coordinate_Exception_Second(Monomer) 
    X_Coordinate = Lattice_Size
    Y_Coordinate = Lattice_Size
    Z_Coordinate = Retrieve_Z_Coordinate(Monomer)
    return X_Coordinate , Y_Coordinate , Z_Coordinate 
end

"""
Movement_Seventeen_Coordinate_Exception_Third(Monomer)

Reflects Y value across X for periodic logic.

# Returns
- `Tuple{Float64, Float64, Float64}`
"""

function Movement_Seventeen_Coordinate_Exception_Third(Monomer) 
    X_Coordinate = Lattice_Size - Retrieve_Y_Coordinate(Monomer)
    Y_Coordinate = Lattice_Size
    Z_Coordinate = Retrieve_Z_Coordinate(Monomer)
    return X_Coordinate , Y_Coordinate , Z_Coordinate 
end

"""
Movement_Seventeen_Coordinate_Exception_Fourth(Monomer)

Reflects X value across Y for periodic logic.

# Returns
- `Tuple{Float64, Float64, Float64}`
"""

function Movement_Seventeen_Coordinate_Exception_Fourth(Monomer) 
    X_Coordinate = Lattice_Size 
    Y_Coordinate = Lattice_Size - Retrieve_X_Coordinate(Monomer)
    Z_Coordinate = Retrieve_Z_Coordinate(Monomer)
    return X_Coordinate , Y_Coordinate , Z_Coordinate 
end

"""
Movement_Seventeen_Coordinate(Monomer)

Calculates the new coordinates after a diagonal movement of -0.5 in both X and Y directions.

# Arguments
- `Monomer::Tuple{Float64, Float64, Float64}`: The current coordinates of the monomer.

# Returns
- `Tuple{Float64, Float64, Float64}`: The updated coordinates after the diagonal movement in the XY plane.
"""

function Movement_Seventeen_Coordinate(Monomer)  
    X_Coordinate = Retrieve_X_Coordinate(Monomer) - .5
    Y_Coordinate = Retrieve_Y_Coordinate(Monomer) - .5
    Z_Coordinate = Retrieve_Z_Coordinate(Monomer)
    return X_Coordinate , Y_Coordinate , Z_Coordinate 
end

"""
Movement_Eighteen_Coordinate_Exception(Monomer)

Handles boundary conditions for diagonal movement in the XY plane during Movement Eighteen when near corners or edges.

# Arguments
- `Monomer::Tuple{Float64, Float64, Float64}`: The current coordinates of the monomer.

# Returns
- `Tuple{Float64, Float64, Float64}`: The coordinates adjusted for edge or corner boundary conditions.
"""

function Movement_Eighteen_Coordinate_Exception(Monomer) 
    X_Coordinate = Retrieve_X_Coordinate(Monomer)
    Y_Coordinate = Retrieve_Y_Coordinate(Monomer)
    Z_Coordinate = Retrieve_Z_Coordinate(Monomer)
    return X_Coordinate , Y_Coordinate , Z_Coordinate 
end

"""
Movement_Eighteen_Coordinate_Exception_Second(Monomer)

Handles full diagonal wrap from both X and Y boundaries.

# Returns
- `Tuple{Float64, Float64, Float64}`
"""

function Movement_Eighteen_Coordinate_Exception_Second(Monomer) 
    X_Coordinate = 0
    Y_Coordinate = 0
    Z_Coordinate = Retrieve_Z_Coordinate(Monomer)
    return X_Coordinate , Y_Coordinate , Z_Coordinate 
end

"""
Movement_Eighteen_Coordinate_Exception_Third(Monomer)

Reflects Y into X at boundary for diagonal shift.

# Returns
- `Tuple{Float64, Float64, Float64}`
"""

function Movement_Eighteen_Coordinate_Exception_Third(Monomer) 
    X_Coordinate = Lattice_Size - Retrieve_Y_Coordinate(Monomer)
    Y_Coordinate = 0
    Z_Coordinate = Retrieve_Z_Coordinate(Monomer)
    return X_Coordinate , Y_Coordinate , Z_Coordinate 
end

"""
Movement_Eighteen_Coordinate_Exception_Fourth(Monomer)

Reflects X into Y at boundary for diagonal shift.

# Returns
- `Tuple{Float64, Float64, Float64}`
"""
function Movement_Eighteen_Coordinate_Exception_Fourth(Monomer) 
    X_Coordinate = 0
    Y_Coordinate = Lattice_Size - Retrieve_X_Coordinate(Monomer)
    Z_Coordinate = Retrieve_Z_Coordinate(Monomer)
    return X_Coordinate , Y_Coordinate , Z_Coordinate 
end

"""
Movement_Eighteen_Coordinate(Monomer)

Calculates the new coordinates after a diagonal movement of +0.5 in both X and Y directions.

# Arguments
- `Monomer::Tuple{Float64, Float64, Float64}`: The current coordinates of the monomer.

# Returns
- `Tuple{Float64, Float64, Float64}`: The updated coordinates after the diagonal movement in the XY plane.
"""

function Movement_Eighteen_Coordinate(Monomer) 
    X_Coordinate = Retrieve_X_Coordinate(Monomer) + .5
    Y_Coordinate = Retrieve_Y_Coordinate(Monomer) + .5
    Z_Coordinate = Retrieve_Z_Coordinate(Monomer)
    return X_Coordinate , Y_Coordinate , Z_Coordinate 
end

"""
MovementFunctions

A dictionary that maps movement direction labels (as strings) to a list of movement rules. Each rule is represented as a tuple of a condition function and a corresponding movement function.

This structure is used to determine how a monomer should move on an FCC lattice while handling periodic boundary conditions and lattice edge exceptions. The conditions are evaluated in order, and the first condition that returns `true` determines which movement function is applied.

# Structure
- Keys: `String` labels for each movement direction (e.g., `"One"`, `"Two"`, ..., `"Eighteen"`).
- Values: A `Vector` of tuples, each containing:
  - A predicate function `(X, Y, Z) -> Bool` that checks a spatial condition.
  - A movement function `f(Monomer)` that returns new coordinates based on that condition.

# Example
```julia
MovementFunctions["One"] = [
    ((X, Y, Z) -> X == Lattice_Size, Movement_One_Coordinate_Exception),
    ((X, Y, Z) -> X == Lattice_Size - 0.5, Movement_One_Coordinate_Exception_Second),
    ((X, Y, Z) -> true, Movement_One_Coordinate)
]
"""

# Define a dictionary mapping movement types to functions and conditions
const MovementFunctions = Dict(
    "One" => [
        ((X, Y, Z) -> X == Lattice_Size, Movement_One_Coordinate_Exception),
        ((X, Y, Z) -> X == Lattice_Size - 0.5, Movement_One_Coordinate_Exception_Second),
        ((X, Y, Z) -> true, Movement_One_Coordinate)
    ],
    "Two" => [
        ((X, Y, Z) -> Y == Lattice_Size, Movement_Two_Coordinate_Exception),
        ((X, Y, Z) -> Y == Lattice_Size - 0.5, Movement_Two_Coordinate_Exception_Second),
        ((X, Y, Z) -> true, Movement_Two_Coordinate)
    ],
    "Three" => [
        ((X, Y, Z) -> Z == Lattice_Size, Movement_Three_Coordinate_Exception),
        ((X, Y, Z) -> Z == Lattice_Size - 0.5, Movement_Three_Coordinate_Exception_Second),
        ((X, Y, Z) -> true, Movement_Three_Coordinate)
    ],
    "Four" => [
        ((X, Y, Z) -> X == 0, Movement_Four_Coordinate_Exception),
        ((X, Y, Z) -> X == 0.5, Movement_Four_Coordinate_Exception_Second),
        ((X, Y, Z) -> true, Movement_Four_Coordinate)
    ],
    "Five" => [
        ((X, Y, Z) -> Y == 0, Movement_Five_Coordinate_Exception),
        ((X, Y, Z) -> Y == 0.5, Movement_Five_Coordinate_Exception_Second),
        ((X, Y, Z) -> true, Movement_Five_Coordinate)
    ],
    "Six" => [
        ((X, Y, Z) -> Z == 0, Movement_Six_Coordinate_Exception),
        ((X, Y, Z) -> Z == 0.5, Movement_Six_Coordinate_Exception_Second),
        ((X, Y, Z) -> true, Movement_Six_Coordinate)
    ],
    "Seven" => [
        ((X, Y, Z) -> (X == 0 && Y == 0) || (X == Lattice_Size && Y == Lattice_Size), Movement_Seven_Coordinate_Exception),
        ((X, Y, Z) -> (X == Lattice_Size && Y < Lattice_Size && Y > 0) || (Y == 0 && X > 0 && X < Lattice_Size), Movement_Seven_Coordinate_Exception_Second),
        ((X, Y, Z) -> (X == Lattice_Size && Y == 0), Movement_Seven_Coordinate_Exception_Third),
        ((X, Y, Z) -> true, Movement_Seven_Coordinate)
    ],
    "Eight" => [
        ((X, Y, Z) -> (X == 0 && Y == 0) || (X == Lattice_Size && Y == Lattice_Size), Movement_Eight_Coordinate_Exception),
        ((X, Y, Z) -> (X == 0 && Y > 0 && Y < Lattice_Size) || (X < Lattice_Size && X > 0 && Y == Lattice_Size), Movement_Eight_Coordinate_Exception_Second),
        ((X, Y, Z) -> (X == 0 && Y == Lattice_Size), Movement_Eight_Coordinate_Exception_Third),
        ((X, Y, Z) -> true, Movement_Eight_Coordinate)
    ],
    "Nine" => [
        ((X, Y, Z) -> (Y == 0 && Z == 0) || (Y == Lattice_Size && Z == Lattice_Size), Movement_Nine_Coordinate_Exception),
        ((X, Y, Z) -> (Y > 0 && Y < Lattice_Size && Z == 0) || (Y == Lattice_Size && Z < Lattice_Size && Z > 0), Movement_Nine_Coordinate_Exception_Second),
        ((X, Y, Z) -> (Y == Lattice_Size && Z == 0), Movement_Nine_Coordinate_Exception_Third),
        ((X, Y, Z) -> true, Movement_Nine_Coordinate)
    ],
    "Ten" => [
        ((X, Y, Z) -> (Y == 0 && Z == 0) || (Y == Lattice_Size && Z == Lattice_Size), Movement_Ten_Coordinate_Exception),
        ((X, Y, Z) -> (Y == 0 && Z > 0 && Z < Lattice_Size) || (Z == Lattice_Size && Y < Lattice_Size && Y > 0), Movement_Ten_Coordinate_Exception_Second),
        ((X, Y, Z) -> (Z == Lattice_Size && Y == 0), Movement_Ten_Coordinate_Exception_Third),
        ((X, Y, Z) -> true, Movement_Ten_Coordinate)
    ],
    "Eleven" => [
        ((X, Y, Z) -> (X == 0 && Z == 0) || (X == Lattice_Size && Z == Lattice_Size), Movement_Eleven_Coordinate_Exception),
        ((X, Y, Z) -> (Z == 0 && X > 0 && X < Lattice_Size) || (X == Lattice_Size && Z < Lattice_Size && Z > 0), Movement_Eleven_Coordinate_Exception_Second),
        ((X, Y, Z) -> (X == Lattice_Size && Z == 0), Movement_Eleven_Coordinate_Exception_Third),
        ((X, Y, Z) -> true, Movement_Eleven_Coordinate)
    ],
    "Twelve" => [
        ((X, Y, Z) -> (X == 0 && Z == 0) || (X == Lattice_Size && Z == Lattice_Size), Movement_Twelve_Coordinate_Exception),
        ((X, Y, Z) -> (X == 0 && Z > 0 && Z < Lattice_Size) || (Z == Lattice_Size && X < Lattice_Size && X > 0), Movement_Twelve_Coordinate_Exception_Second),
        ((X, Y, Z) -> (X == 0 && Z == Lattice_Size), Movement_Twelve_Coordinate_Exception_Third),
        ((X, Y, Z) -> true, Movement_Twelve_Coordinate)
    ],
    "Thirteen" => [
        ((X, Y, Z) -> (X == Lattice_Size && Z == 0) || (X == 0 && Z == Lattice_Size), Movement_Thirteen_Coordinate_Exception),
        ((X, Y, Z) -> (X == 0 && Z > 0 && Z < Lattice_Size), Movement_Thirteen_Coordinate_Exception_Second),
        ((X, Y, Z) -> (Z == 0 && X > 0 && X < Lattice_Size), Movement_Thirteen_Coordinate_Exception_Third),
        ((X, Y, Z) -> (X == 0 && Z == 0), Movement_Thirteen_Coordinate_Exception_Fourth),
        ((X, Y, Z) -> true, Movement_Thirteen_Coordinate)
    ],
    "Fourteen" => [
        ((X, Y, Z) -> (X == Lattice_Size && Z == 0) || (X == 0 && Z == Lattice_Size), Movement_Fourteen_Coordinate_Exception),
        ((X, Y, Z) -> (Z == Lattice_Size && X > 0 && X < Lattice_Size), Movement_Fourteen_Coordinate_Exception_Second),
        ((X, Y, Z) -> (X == Lattice_Size && Z > 0 && Z < Lattice_Size), Movement_Fourteen_Coordinate_Exception_Third),
        ((X, Y, Z) -> (X == Lattice_Size && Z == Lattice_Size), Movement_Fourteen_Coordinate_Exception_Fourth),
        ((X, Y, Z) -> true, Movement_Fourteen_Coordinate)
    ],
    "Fifteen" => [
        ((X, Y, Z) -> (Y == 0 && Z == Lattice_Size) || (Y == Lattice_Size && Z == 0), Movement_Fifteen_Coordinate_Exception),
        ((X, Y, Z) -> (Y == 0 && Z > 0 && Z < Lattice_Size), Movement_Fifteen_Coordinate_Exception_Second),
        ((X, Y, Z) -> (Z == 0 && Y > 0 && Y < Lattice_Size), Movement_Fifteen_Coordinate_Exception_Third),
        ((X, Y, Z) -> (Y == 0 && Z == 0), Movement_Fifteen_Coordinate_Exception_Fourth),
        ((X, Y, Z) -> true, Movement_Fifteen_Coordinate)
    ],
    "Sixteen" => [
        ((X, Y, Z) -> (Y == Lattice_Size && Z == 0) || (Y == 0 && Z == Lattice_Size), Movement_Sixteen_Coordinate_Exception),
        ((X, Y, Z) -> (Y == Lattice_Size && Z == Lattice_Size), Movement_Sixteen_Coordinate_Exception_Second),
        ((X, Y, Z) -> (Z == Lattice_Size), Movement_Sixteen_Coordinate_Exception_Third),
        ((X, Y, Z) -> (Y == Lattice_Size), Movement_Sixteen_Coordinate_Exception_Fourth),
        ((X, Y, Z) -> true, Movement_Sixteen_Coordinate)
    ],
    "Seventeen" => [
        ((X, Y, Z) -> (X == Lattice_Size && Y == 0) || (X == 0 && Y == Lattice_Size), Movement_Seventeen_Coordinate_Exception),
        ((X, Y, Z) -> (X == 0 && Y == 0), Movement_Seventeen_Coordinate_Exception_Second),
        ((X, Y, Z) -> (X == 0), Movement_Seventeen_Coordinate_Exception_Third),
        ((X, Y, Z) -> (Y == 0), Movement_Seventeen_Coordinate_Exception_Fourth),
        ((X, Y, Z) -> true, Movement_Seventeen_Coordinate)
    ],
    "Eighteen" => [
        ((X, Y, Z) -> (X == Lattice_Size && Y == 0) || (X == 0 && Y == Lattice_Size), Movement_Eighteen_Coordinate_Exception),
        ((X, Y, Z) -> (X == Lattice_Size && Y == Lattice_Size), Movement_Eighteen_Coordinate_Exception_Second),
        ((X, Y, Z) -> (X == Lattice_Size), Movement_Eighteen_Coordinate_Exception_Third),
        ((X, Y, Z) -> (Y == Lattice_Size), Movement_Eighteen_Coordinate_Exception_Fourth),
        ((X, Y, Z) -> true, Movement_Eighteen_Coordinate)
    ]
) 
    let
    valid_sampled_moves = Set(filter(x -> x != "None", Possible_Movement_Options))
    dispatch_moves = Set(keys(MovementFunctions))

    missing_from_dispatch = setdiff(valid_sampled_moves, dispatch_moves)
    extra_in_dispatch = setdiff(dispatch_moves, valid_sampled_moves)

    if !isempty(missing_from_dispatch) || !isempty(extra_in_dispatch)
        error("""
        Movement option / dispatch mismatch detected.

        Missing from MovementFunctions: $(collect(missing_from_dispatch))
        Extra in MovementFunctions: $(collect(extra_in_dispatch))
        """)
    end
end

"""
Distinguishing_Monomers(Monomer, Desired_Location, Type_of_Movement)

Distinguishes between different types of monomers and calls the appropriate movement/interaction function.

# Arguments
- `Monomer::Tuple{Float64, Float64, Float64}`: The current coordinates of the monomer.
- `Desired_Location::Tuple{Float64, Float64, Float64}`: The desired location for the monomer.
- `Type_of_Movement::String`: The type of movement being attempted.

# Calls
- `Retrieve_State_Monomer`
- `Native_Move_or_Conformational_Change`
- `AggregateProne_Aggregation_or_Movement`
- `Oligomer_Move`
- `Fibril_Move`
"""

function Distinguishing_Monomers(Monomer, Desired_Location, Type_of_Movement)
    State= Retrieve_State_Monomer(Monomer)

    if State == 1
        Native_Move_or_Conformational_Change(Monomer, Desired_Location)
    elseif State == 2
        AggregateProne_Aggregation_or_Movement(Monomer, Desired_Location)
    elseif State == 3
        Oligomer_Move(Monomer, Desired_Location, Type_of_Movement)
    elseif State == 4
        Fibril_Move(Monomer, Desired_Location, Type_of_Movement)
    end

end

"""
Key_Array_Locations_and_States()

Returns an array of all keys (coordinates) in the Locations_and_States_Dict.

# Global variables read
- `Locations_and_States_Dict`

# Returns
- `Vector{Tuple{Float64, Float64, Float64}}`: An array of all coordinates in the dictionary.
"""

function Key_Array_Locations_and_States()
    return collect(keys(Locations_and_States_Dict))
end

"""
Monomer_State(State)

Returns the string representation of a monomer state.

# Arguments
- `State::Int`: The integer representation of the monomer state.

# Returns
- `String`: The string representation of the monomer state.

# Throws
- `ErrorException`: If an unknown state is encountered.
"""
    
function Monomer_State(State)
    if State == 1
        return "Native"
    elseif State == 2
        return "AggregateProne"
    elseif State == 3
        return "Oligomer"
    elseif State == 4
        return "Fibril"
    elseif State == 5
        return "Sphere"
    elseif State == 0
        return "Empty"
    else
        error("Unknown state: $State")
    end
end

"""
Monomer_Availability(Monomer)

Checks if a given location is available (empty).

# Arguments
- `Monomer::Tuple{Float64, Float64, Float64}`: The coordinates of the location to check.

# Global variables read
- `Locations_and_States_Dict`

# Returns
- `Bool`: `true` if the location is available (state is 0), `false` otherwise.
"""

function Monomer_Availability(Monomer)
    # Check if the Monomer coordinate exists in the dictionary
    if haskey(Locations_and_States_Dict, Monomer)
        State, _ = Locations_and_States_Dict[Monomer]

        # Check if the state is non-zero (indicating it's occupied)
        return State == 0  # true if available, false if occupied
    else
        # If the coordinate is not found in the dictionary, consider it available (as no state is assigned)
        return true
    end
end

"""
One_Monomer_Movement(Monomer, Desired_Location)

Moves a monomer from one location to another.

# Arguments
- `Monomer::Tuple{Float64, Float64, Float64}`: The current coordinates of the monomer.
- `Desired_Location::Tuple{Float64, Float64, Float64}`: The desired coordinates for the monomer.

# Global variables modified
- `Locations_and_States_Dict`

# Calls
- `Update_Locations_States` (twice)
"""

function One_Monomer_Movement(Monomer, Desired_Location)
   # Retrieve the state and unique number of the monomer once
   monomer_state, monomer_unique_number = Locations_and_States_Dict[Monomer]

   Update_Locations_States(Desired_Location, monomer_state, monomer_unique_number)

   # Set the original location's state to 0 while keeping the unique number unchanged
   Update_Locations_States(Monomer, 0, 0)
end

"""
Empties_Locations_and_States(Monomer)

Sets the state of a location to 0 (empty).

# Arguments
- `Monomer::Tuple{Float64, Float64, Float64}`: The coordinates of the location to empty.

# Global variables modified
- `Locations_and_States_Dict`
"""

function Empties_Locations_and_States(Monomer)
    Locations_and_States_Dict[Monomer] = (0, Locations_and_States_Dict[Monomer][2])
end

"""
Conformational_Change(Monomer)

Changes the state of a monomer between Native (1) and AggregateProne (2).

# Arguments
- `Monomer::Tuple{Float64, Float64, Float64}`: The coordinates of the monomer.

# Global variables modified
- `Locations_and_States_Dict`

# Calls
- `Update_Locations_States`
"""

function Conformational_Change(Monomer)
    # Retrieve the state and unique number from the dictionary using the Monomer coordinates
    State, Unique_Number = Locations_and_States_Dict[Monomer]

    if State == 1
        # Change state from Native (1) to AggregateProne (2)
        Update_Locations_States(Monomer, 2, Unique_Number)
    elseif State == 2
        # Change state from AggregateProne (2) to Native (1)
        Update_Locations_States(Monomer, 1, Unique_Number)
    end


end

"""
Native_Move_or_Conformational_Change(Monomer, Desired_Location)

Handles the movement or conformational change of a Native monomer.

# Arguments
- `Monomer::Tuple{Float64, Float64, Float64}`: The current coordinates of the monomer.
- `Desired_Location::Tuple{Float64, Float64, Float64}`: The desired coordinates for the monomer.

# Global variables read
- `Native_to_AggregateProne`

# Calls
- `Retrieve_State_Monomer`
- `One_Monomer_Movement`
- `Random_Dice`
- `Conformational_Change`
"""

function Native_Move_or_Conformational_Change(Monomer, Desired_Location)

    State = Retrieve_State_Monomer(Desired_Location)

    if State == 0  # Available if state is empty (0)
        One_Monomer_Movement(Monomer, Desired_Location)
        # Only allow conformational change within this block
        if Random_Dice() < Native_to_AggregateProne
            Conformational_Change(Desired_Location)
        end
    else
        # Handle conformational change only if movement doesn't happen
        if Random_Dice() < Native_to_AggregateProne
            Conformational_Change(Monomer)
        end
    end

end

"""
AggregateProne_Aggregation_or_Movement(Monomer, Desired_Location)

Handles the aggregation or movement of an AggregateProne monomer.

# Arguments
- `Monomer::Tuple{Float64, Float64, Float64}`: The current coordinates of the monomer.
- `Desired_Location::Tuple{Float64, Float64, Float64}`: The desired coordinates for the monomer.

# Global variables read
- `AggregateProne_to_Native`

# Calls
- `Retrieve_State_Monomer`
- `Monomer_State`
- `One_Monomer_Movement`
- `Random_Dice`
- `Conformational_Change`
- `AggregateProne_Aggregation`
"""

function AggregateProne_Aggregation_or_Movement(Monomer, Desired_Location)

    Status_Desired_Location = Retrieve_State_Monomer(Desired_Location)
    location_state = Monomer_State(Status_Desired_Location)

    if Status_Desired_Location == 0
        # Case 1: Move the monomer to an empty location
        One_Monomer_Movement(Monomer, Desired_Location)
        if Random_Dice() < AggregateProne_to_Native
            Conformational_Change(Desired_Location)
        end
    elseif Status_Desired_Location == 2 || Status_Desired_Location == 3
        # Case 2: Handle AggregateProne aggregation
        AggregateProne_Aggregation(Monomer, Desired_Location, location_state)
        Status_Monomer = Retrieve_State_Monomer(Monomer)
        if Status_Monomer == 2 && Random_Dice() < AggregateProne_to_Native
            Conformational_Change(Monomer)
        end
    elseif Random_Dice() < AggregateProne_to_Native #If the state is 1, 4, or 5
        # Case 3: Apply conformational change directly if no other conditions are met
        Conformational_Change(Monomer)
    end

end

"""
Unique_Number_Generator!(available_numbers)

Generates a unique number for a new aggregate.

# Arguments
- `available_numbers::Vector{Int64}`: A vector of available unique numbers.

# Returns
- `Int64`: A unique number.

# Throws
- `ErrorException`: If there are no more unique numbers available.
"""

function Unique_Number_Generator!(available_numbers)
    if length(available_numbers) == 0
        error("No more unique numbers available!")
    end

    # Randomly select and remove a number from the available list
    idx = rand(1:length(available_numbers))
    UniqueCode = available_numbers[idx]
    deleteat!(available_numbers, idx)  # Remove the selected number to ensure uniqueness
    
    return UniqueCode
end

"""
No_Repeating_Unique_Number_Locked(UniqueCode)

Checks if a unique number is already in use.

# Arguments
- `UniqueCode::Int64`: The unique number to check.

# Global variables read
- `Locations_and_States_Dict`

# Returns
- `Bool`: `true` if the unique number is already in use, `false` otherwise.
"""

function No_Repeating_Unique_Number_Locked(UniqueCode)
    # Iterate through the values of the dictionary to check for the unique code
    for (_, unique_number) in values(Locations_and_States_Dict)
        if unique_number == UniqueCode
            return true
        end
    end
    return false
end

"""
AggregateProne_AggregateProne_Lock(Monomer, Desired_Location)

Handles the aggregation of two AggregateProne monomers into an Oligomer.

# Arguments
- `Monomer::Tuple{Float64, Float64, Float64}`: The coordinates of the first AggregateProne monomer.
- `Desired_Location::Tuple{Float64, Float64, Float64}`: The coordinates of the second AggregateProne monomer.

# Global variables modified
- `Locations_and_States_Dict`
- `Initial_Locations_and_States_Dict`
- `available_numbers`

# Calls
- `Unique_Number_Generator!`
- `Retrieve_Unique_Number_Monomer` (twice)
- `Update_Locations_States` (twice)
- `Remove_Center_of_Mass_Info` (twice)
- `Calculate_Center_of_Mass`
- `Initialize_New_Center_of_Mass`
"""

function AggregateProne_AggregateProne_Lock(Monomer, Desired_Location)
    # Generate a unique code for locking
    Unique_Code = Unique_Number_Generator!(available_numbers)

    # Retrieve the unique number of both the monomer and the desired location
    Monomer_Unique_Number = Retrieve_Unique_Number_Monomer(Monomer)
    Desired_Location_Unique_Number = Retrieve_Unique_Number_Monomer(Desired_Location)

    # Update the states to Oligomer (3) and set the unique number for both locations
    Update_Locations_States(Monomer, 3, Unique_Code)
    Update_Locations_States(Desired_Location, 3, Unique_Code)

    #Delete the info of monomer and the desired location monomer from Initial_Locations_and_States_Dict
    Remove_Center_of_Mass_Info(Monomer_Unique_Number)
    Remove_Center_of_Mass_Info(Desired_Location_Unique_Number)

    #Calculate the new center of mass for the new oligomer
    X, Y, Z = Calculate_Center_of_Mass(Monomer, Desired_Location)

    #Input that information with its unique code into Initial_Locations_and_States_Dict
    Initialize_New_Center_of_Mass(X, Y, Z, 3, Unique_Code)
end

"""
Calculate_Center_of_Mass(Monomer, Desired_Location)

Calculates the center of mass of two monomers.

# Arguments
- `Monomer::Tuple{Float64, Float64, Float64}`: The coordinates of the first monomer.
- `Desired_Location::Tuple{Float64, Float64, Float64}`: The coordinates of the second monomer.

# Returns
- `Tuple{Float64, Float64, Float64}`: The (X, Y, Z) coordinates of the center of mass.
"""

function Calculate_Center_of_Mass(Monomer, Desired_Location)
    Center_of_Mass_X = (Monomer[1] + Desired_Location[1]) / 2
    Center_of_Mass_Y = (Monomer[2] + Desired_Location[2]) / 2
    Center_of_Mass_Z = (Monomer[3] + Desired_Location[3]) / 2
    return Center_of_Mass_X, Center_of_Mass_Y, Center_of_Mass_Z
end

"""
Initialize_New_Center_of_Mass(X_Coordinate, Y_Coordinate, Z_Coordinate, State_Value, Unique_Code)

Initializes the center of mass information for a new aggregate.

# Arguments
- `X_Coordinate::Float64`: The X coordinate of the center of mass.
- `Y_Coordinate::Float64`: The Y coordinate of the center of mass.
- `Z_Coordinate::Float64`: The Z coordinate of the center of mass.
- `State_Value::Int64`: The state value of the aggregate (3 for Oligomer, 4 for Fibril).
- `Unique_Code::Int64`: The unique code of the aggregate.

# Global variables modified
- `Initial_Locations_and_States_Dict`
"""

function Initialize_New_Center_of_Mass(X_Coordinate, Y_Coordinate, Z_Coordinate, State_Value, Unique_Code)
    Initial_Locations_and_States_Dict[X_Coordinate, Y_Coordinate, Z_Coordinate] = (State_Value, Unique_Code)
end

"""
Delete_Monomer_Information_from_Initial_Locations_and_States(unique_number, state)

Deletes the information of a monomer from the Initial_Locations_and_States_Dict.

# Arguments
- `unique_number::Int64`: The unique number of the monomer to delete.
- `state::Int64`: The state of the monomer to delete.

# Global variables modified
- `Initial_Locations_and_States_Dict`
"""

function Delete_Monomer_Information_from_Initial_Locations_and_States(unique_number, state)
     # Find the coordinate that corresponds to the unique number
     monomer_to_delete = nothing

     for (coordinate, (current_state, uid)) in Initial_Locations_and_States_Dict
         if uid == unique_number && current_state == state
             monomer_to_delete = coordinate
             break  # Stop after finding the first match
         end
     end
 
     # If a matching coordinate is found, delete it
     if monomer_to_delete !== nothing
         delete!(Initial_Locations_and_States_Dict, monomer_to_delete)
     end
end

"""
Retrieve_Unique_Number_Monomer(Monomer)

Retrieves the unique number associated with a monomer.

# Arguments
- `Monomer::Tuple{Float64, Float64, Float64}`: The coordinates of the monomer.

# Global variables read
- `Locations_and_States_Dict`

# Returns
- `Int64`: The unique number of the monomer.
"""

function Retrieve_Unique_Number_Monomer(Monomer)
     # Get the unique number associated with this Monomer (coordinate)
     _, Unique_Number = Locations_and_States_Dict[Monomer]

     return Unique_Number
end

"""
AggregateProne_Aggregation(Monomer, Desired_Location, Status_Desired_Location)

Handles the aggregation of an AggregateProne monomer with another monomer or aggregate.

# Arguments
- `Monomer::Tuple{Float64, Float64, Float64}`: The coordinates of the AggregateProne monomer.
- `Desired_Location::Tuple{Float64, Float64, Float64}`: The coordinates of the other monomer or aggregate.
- `Status_Desired_Location::String`: The state of the other monomer or aggregate.

# Global variables read
- `Oligomer_Formation`
- `Fibril_Formation`

# Calls
- `Random_Dice` (twice)
- `AggregateProne_AggregateProne_Lock`
- `AggregateProne_Oligomer_Lock_Two`
"""

function AggregateProne_Aggregation(Monomer, Desired_Location, Status_Desired_Location)
    # Check if the status of the desired location is "AggregateProne" and if conditions meet for forming an oligomer
    if Status_Desired_Location == "AggregateProne" && Random_Dice() < Oligomer_Formation && Monomer != Desired_Location #If you are next to a 2
        # Perform AggregateProne-AggregateProne locking
        AggregateProne_AggregateProne_Lock(Monomer, Desired_Location)
    elseif Status_Desired_Location == "Oligomer" && Random_Dice() < Fibril_Formation #If you are next to a 3
        # Perform AggregateProne-Oligomer locking
        AggregateProne_Oligomer_Lock_Two(Monomer, Desired_Location)
    end
end

"""
AggregateProne_Oligomer_Lock_Two(Monomer_AggregateProne, Desired_Location_Oligomer)

Handles the aggregation of an AggregateProne monomer with an Oligomer to form a Fibril.

# Arguments
- `Monomer_AggregateProne::Tuple{Float64, Float64, Float64}`: The coordinates of the AggregateProne monomer.
- `Desired_Location_Oligomer::Tuple{Float64, Float64, Float64}`: The coordinates of the Oligomer.

# Global variables modified
- `Locations_and_States_Dict`
- `Initial_Locations_and_States_Dict`

# Calls
- `Retrieve_Unique_Number_Monomer` (twice)
- `Key_Array_Locations_and_States`
- `Retrieve_State_Monomer`
- `Update_Locations_States` (multiple times)
- `Remove_Center_of_Mass_Info` (twice)
- `Calculate_Center_of_Mass_Fibril`
- `Initialize_New_Center_of_Mass`
"""

function AggregateProne_Oligomer_Lock_Two(Monomer_AggregateProne, Desired_Location_Oligomer)
    # Retrieve the unique code associated with the AggregateProne and oligomer
    Unique_Code_Oligomer = Retrieve_Unique_Number_Monomer(Desired_Location_Oligomer)
    Unique_Code_AggregateProne = Retrieve_Unique_Number_Monomer(Monomer_AggregateProne)

    Key_Array = Filter_Oligomer()

    # Iterate over all coordinates in the dictionary to find matching oligomers
    @threads for i in 1:length(Key_Array)
        #for i in 1:length(Key_Array)
        Current_Coordinate =  Key_Array[i]
        Unique_Code = Retrieve_Unique_Number_Monomer(Current_Coordinate)
        State = Retrieve_State_Monomer(Current_Coordinate)

        # Check if the unique code matches and the state is 3 (Oligomer)
        if Unique_Code == Unique_Code_Oligomer && State == 3
            # Update the state to 4 (Fibril)
            Update_Locations_States(Current_Coordinate, 4, Unique_Code_Oligomer)
        end
    end

    # Update the state of the AggregateProne monomer to 4 (Fibril) with the same unique code
    Update_Locations_States(Monomer_AggregateProne, 4, Unique_Code_Oligomer)

    #Delete the info of oligomer and AggregateProne from Initial_Locations_and_States_Dict
    Remove_Center_of_Mass_Info(Unique_Code_Oligomer)
    Remove_Center_of_Mass_Info(Unique_Code_AggregateProne)

    #Calculate the new center of mass for the new fibril
    X, Y, Z = Calculate_Center_of_Mass_Fibril(Unique_Code_Oligomer)

    #Calculate new center of mass for new fibril
    Initialize_New_Center_of_Mass(X, Y, Z, 4, Unique_Code_Oligomer)
end

"""
Oligomer_Move(Monomer, Desired_Location, Type_of_Movement)

Handles the movement of an Oligomer.

# Arguments
- `Monomer::Tuple{Float64, Float64, Float64}`: The coordinates of the Oligomer.
- `Desired_Location::Tuple{Float64, Float64, Float64}`: The desired coordinates for the Oligomer.
- `Type_of_Movement::String`: The type of movement being attempted.

# Global variables read
- `Oligomer_Dissociation_rate`

# Calls
- `Retrieve_State_Monomer` (twice)
- `Retrieve_Unique_Number_Monomer`
- `Gather_all_Aggregate_Monomers_Oligomer`
- `Random_Dice`
- `Oligomer_Dissociation`
- `Oligomer_Aggregate`
"""

function Oligomer_Move(Monomer, Desired_Location, Type_of_Movement)
     # Check the state of the desired location

    State_Desired_Location = Retrieve_State_Monomer(Desired_Location)
    Unique_Number_Monomer = Retrieve_Unique_Number_Monomer(Monomer)

    # If the desired location is empty, move oligomer
    if State_Desired_Location == 0
       Gather_all_Aggregate_Monomers_Oligomer(Monomer, Type_of_Movement)
       #Unique_Number_Desired_Location = Retrieve_Unique_Number_Monomer(Desired_Location)
       if Random_Dice() < Oligomer_Dissociation_rate #&& Unique_Number_Monomer == Unique_Number_Desired_Location #Check and see if it will dissociate
            Oligomer_Dissociation(Unique_Number_Monomer)
       end
    else #If it is next to any monomer it will go in here
        Oligomer_Aggregate(Monomer, Desired_Location)
    end

end

"""
Oligomer_Aggregate(Monomer, Desired_Location)

Handles the aggregation or dissociation of an Oligomer.

# Arguments
- `Monomer::Tuple{Float64, Float64, Float64}`: The coordinates of the Oligomer.
- `Desired_Location::Tuple{Float64, Float64, Float64}`: The coordinates of the other monomer or aggregate.

# Global variables read
- `Fibril_Formation`
- `Oligomer_Dissociation_rate`

# Calls
- `Retrieve_State_Monomer` (twice)
- `Random_Dice` (twice)
- `AggregateProne_Oligomer_Lock`
- `Retrieve_Unique_Number_Monomer`
- `Oligomer_Dissociation`
"""

function Oligomer_Aggregate(Monomer, Desired_Location)
    # Retrieve the states of the monomer and the desired location directly from the dictionary


    State_Desired_Location = Retrieve_State_Monomer(Desired_Location) 
    State_Monomer = Retrieve_State_Monomer(Monomer) 

    # Check for the conditions for aggregation or dissociation
    if State_Desired_Location == 2 && Random_Dice() < Fibril_Formation
        AggregateProne_Oligomer_Lock(Monomer, Desired_Location)
    elseif Random_Dice() < Oligomer_Dissociation_rate && State_Monomer == 3
        Unique_Number_Oligomer = Retrieve_Unique_Number_Monomer(Monomer)
        Oligomer_Dissociation(Unique_Number_Oligomer)
    end

end

"""
Oligomer_Dissociation(Unique_Number_Oligomer)

Handles the dissociation of an Oligomer into AggregateProne monomers.

# Arguments
- `Unique_Number_Oligomer::Int64`: The unique number of the Oligomer.

# Global variables modified
- `Locations_and_States_Dict`
- `available_numbers`
- `Initial_Locations_and_States_Dict`

# Calls
- `Retrieve_Unique_Number_Monomer`
- `Key_Array_Locations_and_States`
- `Randomly_Choosing_Unique_Number_Monomer`
- `Update_Locations_States` (multiple times)
- `Remove_Center_of_Mass_Info`
- `Restore_Dissociated_Monomers`
"""

function Oligomer_Dissociation(Unique_Number_Oligomer)
    # Retrieve the unique number associated with the monomer directly from the dictionary

    # Iterate through all coordinates in the dictionary
    Key_Array = Filter_Oligomer()
    @threads for i in 1:length(Key_Array)
        #for i in 1:length(Key_Array)
        coordinate = Key_Array[i]
        Unique_Code = Retrieve_Unique_Number_Monomer(coordinate)
        # If the unique code matches the monomer's unique code, modify the state and unique number
        if Unique_Code == Unique_Number_Oligomer
            Unique_Number_Monomer = Randomly_Choosing_Unique_Number_Monomer()
            Update_Locations_States(coordinate, 2, Unique_Number_Monomer)
            Restore_Dissociated_Monomers(coordinate, Unique_Number_Monomer)
        end
    end
    Remove_Center_of_Mass_Info(Unique_Number_Oligomer)
end

"""
Remove_Center_of_Mass_Info(Unique_Number)

Removes the center of mass information of an aggregate from Initial_Locations_and_States_Dict.

# Arguments
- `Unique_Number::Int64`: The unique number of the aggregate.

# Global variables modified
- `Initial_Locations_and_States_Dict`
"""

function Remove_Center_of_Mass_Info(Unique_Number)
    for (initial_CoM, (state, uid)) in Initial_Locations_and_States_Dict
        if uid == Unique_Number  # Find the corresponding oligomer CoM
            delete!(Initial_Locations_and_States_Dict, initial_CoM)
            break  # Exit once found to avoid unnecessary checks
        end
    end
end

"""
Restore_Dissociated_Monomers(Monomer, Unique_Number)

Restores the state and unique number of dissociated monomers.

# Arguments
- `Monomer::Tuple{Float64, Float64, Float64}`: The coordinates of the dissociated monomer.
- `Unique_Number::Int64`: The unique number of the dissociated monomer.

# Global variables modified
- `Initial_Locations_and_States_Dict`
"""

function Restore_Dissociated_Monomers(Monomer, Unique_Number)
        Initial_Locations_and_States_Dict[Monomer] = (2, Unique_Number)
end

"""
AggregateProne_Oligomer_Lock(Monomer_Oligomer, Desired_Location_AggregateProne)

Handles the aggregation of an Oligomer with an AggregateProne monomer to form a Fibril.

# Arguments
- `Monomer_Oligomer::Tuple{Float64, Float64, Float64}`: The coordinates of the Oligomer.
- `Desired_Location_AggregateProne::Tuple{Float64, Float64, Float64}`: The coordinates of the AggregateProne monomer.

# Global variables modified
- `Locations_and_States_Dict`
- `Initial_Locations_and_States_Dict`

# Calls
- `Retrieve_Unique_Number_Monomer` (twice)
- `Key_Array_Locations_and_States`
- `Update_Locations_States` (multiple times)
- `Remove_Center_of_Mass_Info` (twice)
- `Calculate_Center_of_Mass_Fibril`
- `Initialize_New_Center_of_Mass`
"""

function AggregateProne_Oligomer_Lock(Monomer_Oligomer, Desired_Location_AggregateProne)
    # Retrieve the unique number associated with the oligomer 
    Unique_Code_Oligomer = Retrieve_Unique_Number_Monomer(Monomer_Oligomer)

    # Retrieve the unique number associated with the AggregateProne 
    Unique_Code_AggregateProne = Retrieve_Unique_Number_Monomer(Desired_Location_AggregateProne)

    # Iterate through all coordinates in the dictionary
    #Key_Array = Key_Array_Locations_and_States()
    Key_Array = Filter_Oligomer()

    @threads for i in 1:length(Key_Array)
        #for i in 1:length(Key_Array)
        coordinate = Key_Array[i]
        Unique_Code = Retrieve_Unique_Number_Monomer(coordinate)

        # If the unique code matches the oligomer's unique code, change the state to 4 (Fibril)
        if Unique_Code == Unique_Code_Oligomer
            Update_Locations_States(coordinate, 4, Unique_Code_Oligomer)
        end
    end

    # Update the state of the desired location (AggregateProne) to 4 (Fibril) and set its unique code to the oligomer's unique code
    Update_Locations_States(Desired_Location_AggregateProne, 4, Unique_Code_Oligomer)

    #Remove Info of Oligomer and AggregateProne from Initial_Locations_and_States_Dict
    Remove_Center_of_Mass_Info(Unique_Code_Oligomer)
    Remove_Center_of_Mass_Info(Unique_Code_AggregateProne)

    #Calculate new center of mass for new fibril 
    X, Y, Z = Calculate_Center_of_Mass_Fibril(Unique_Code_Oligomer)

    #Update that information from Initial_Locations_and_States_Dict
    Initialize_New_Center_of_Mass(X, Y, Z, 4, Unique_Code_Oligomer)


end

"""
Calculate_Center_of_Mass_Fibril(Unique_Code)

Calculates the center of mass of a Fibril.

# Arguments
- `Unique_Code::Int64`: The unique code of the Fibril.

# Global variables read
- `Locations_and_States_Dict`

# Global variables modified
- `dict_lock`

# Calls
- `Key_Array_Locations_and_States`
- `Retrieve_Unique_Number_Monomer`

# Returns
- `Tuple{Float64, Float64, Float64}`: The (CoM_x, CoM_y, CoM_z) coordinates of the center of mass.
"""


function Calculate_Center_of_Mass_Fibril(Unique_Code)

    Key_Array = Key_Array_Locations_and_States()

    # Initialize lists to store coordinates
    x_coords = Float64[]
    y_coords = Float64[]
    z_coords = Float64[]

    @threads for i in 1:length(Key_Array)
        Coordinate = Key_Array[i]
        Unique_Code_Coordinate = Retrieve_Unique_Number_Monomer(Coordinate)

        if Unique_Code == Unique_Code_Coordinate
            lock(dict_lock) do
                push!(x_coords, Coordinate[1])
                push!(y_coords, Coordinate[2])
                push!(z_coords, Coordinate[3])
            end
        end
    end
# Compute the center of mass
CoM_x = sum(x_coords) / length(x_coords)
CoM_y = sum(y_coords) / length(y_coords)
CoM_z = sum(z_coords) / length(z_coords)

return (CoM_x, CoM_y, CoM_z)

end

"""
Fibril_Move(Monomer, Desired_Location, Type_of_Movement)

Handles the movement of a Fibril.

# Arguments
- `Monomer::Tuple{Float64, Float64, Float64}`: The coordinates of the Fibril.
- `Desired_Location::Tuple{Float64, Float64, Float64}`: The desired coordinates for the Fibril.
- `Type_of_Movement::String`: The type of movement being attempted.

# Calls
- `Retrieve_State_Monomer`
- `Gather_all_Aggregate_Monomers_Aggregate`
- `Fibril_Aggregate`
"""

function Fibril_Move(Monomer, Desired_Location, Type_of_Movement)
 # Check the state of the desired location
 State_Desired_Location = Retrieve_State_Monomer(Desired_Location)

 # If the desired location is empty, check if the aggregate can move
 if State_Desired_Location == 0
    Gather_all_Aggregate_Monomers_Aggregate(Monomer, Type_of_Movement)
 else 
    Fibril_Aggregate(Monomer, Desired_Location)
 end


end

"""
Fibril_Aggregate(Monomer, Desired_Location)

Handles the aggregation of a Fibril with an AggregateProne monomer.

# Arguments
- `Monomer::Tuple{Float64, Float64, Float64}`: The coordinates of the Fibril.
- `Desired_Location::Tuple{Float64, Float64, Float64}`: The coordinates of the other monomer or aggregate.

# Global variables read
- `Fibril_Growth`

# Calls
- `Retrieve_State_Monomer`
- `Random_Dice`
- `Fibril_Lock`
"""

function Fibril_Aggregate(Monomer, Desired_Location)

 # Retrieve the state of the desired location directly from the dictionary
 State_Desired_Location = Retrieve_State_Monomer(Desired_Location)

 # Check if the desired location is "AggregateProne" and a random chance allows fibril growth
 if State_Desired_Location == 2 && Random_Dice() < Fibril_Growth
     Fibril_Lock(Monomer, Desired_Location)
 end

end

"""
Fibril_Lock(Fibril, Desired_Location_AggregateProne)

Handles the locking of a Fibril with an AggregateProne monomer.

# Arguments
- `Fibril::Tuple{Float64, Float64, Float64}`: The coordinates of the Fibril.
- `Desired_Location_AggregateProne::Tuple{Float64, Float64, Float64}`: The coordinates of the AggregateProne monomer.

# Global variables modified
- `Locations_and_States_Dict`
- `Initial_Locations_and_States_Dict`

# Calls
- `Retrieve_Unique_Number_Monomer` (twice)
- `Update_Locations_States`
- `Remove_Center_of_Mass_Info` (twice)
- `Calculate_Center_of_Mass_Fibril`
- `Initialize_New_Center_of_Mass`
"""

function Fibril_Lock(Fibril, Desired_Location_AggregateProne)
    # Retrieve the unique number associated with the fibril directly from the dictionary
    Unique_Code_Fibril = Retrieve_Unique_Number_Monomer(Fibril)
    Unique_Code_AggregateProne = Retrieve_Unique_Number_Monomer(Desired_Location_AggregateProne)

    # Update the state of the desired location (AggregateProne) to 4 (Fibril) and set its unique code
    Update_Locations_States(Desired_Location_AggregateProne, 4, Unique_Code_Fibril)
    
    #Remove Info from Initial_Locations_and_States_Dict
    Remove_Center_of_Mass_Info(Unique_Code_Fibril)
    Remove_Center_of_Mass_Info(Unique_Code_AggregateProne)

    #Calculate new center of mass with bigger fibril 
    X, Y, Z = Calculate_Center_of_Mass_Fibril(Unique_Code_Fibril)

    #Input that info into Initial_Locations_and_States_Dict
    Initialize_New_Center_of_Mass(X, Y, Z, 4, Unique_Code_Fibril)
end

"""
Random_Dice()

Simulates a random dice roll (returns a random number between 0 and 1).

# Returns
- `Float64`: A random number between 0 and 1.
"""

function Random_Dice()
    Random_Dice = 1-rand()
    return Random_Dice
end

"""
Gather_all_Aggregate_Monomers_Oligomer(Monomer, Type_of_Movement)

Gathers all monomers that belong to the same Oligomer and checks if they can move.

# Arguments
- `Monomer::Tuple{Float64, Float64, Float64}`: The coordinates of a monomer in the Oligomer.
- `Type_of_Movement::String`: The type of movement being attempted.

# Global variables modified
- `Possible_Coordinate_Movements_Dict`

# Global variables read
- `Locations_and_States_Dict`
- `MovementFunctions`

# Calls
- `Empty_Possible_Coordinates_Movement_Dict`
- `Retrieve_Unique_Number_Monomer`
- `Filter_Oligomer`
- `Retrieve_X_Coordinate`
- `Retrieve_Y_Coordinate`
- `Retrieve_Z_Coordinate`
- Movement functions (from `MovementFunctions`)
- `Retrieve_State_Monomer`
- `lock`
- `Appending_Location_Possible_Coordinate_Dict`
- `Move_Aggregate`
"""

function Gather_all_Aggregate_Monomers_Oligomer(Monomer, Type_of_Movement)
    Empty_Possible_Coordinates_Movement_Dict()
    Unique_Code_Monomer = Retrieve_Unique_Number_Monomer(Monomer)
    Movement_Options = MovementFunctions[Type_of_Movement]
    # Filter the dictionary to get only entries with state == 3
    relevant_keys_array = Filter_Oligomer()
    dictionary_emptied = Threads.Atomic{Bool}(false) # Flag to check if dictionary has been emptied

    # Now iterate over the filtered collection with relevant keys only
    @threads for i in 1:length(relevant_keys_array)
        Current_Coordinate = relevant_keys_array[i]
        Unique_Code_Current_Coordinate = Retrieve_Unique_Number_Monomer(Current_Coordinate)

        if Unique_Code_Monomer == Unique_Code_Current_Coordinate
            X_Coordinate, Y_Coordinate, Z_Coordinate = Retrieve_X_Coordinate(Current_Coordinate), Retrieve_Y_Coordinate(Current_Coordinate), Retrieve_Z_Coordinate(Current_Coordinate)

            for (Condition, Movement_Function) in Movement_Options
                if Condition(X_Coordinate, Y_Coordinate, Z_Coordinate)
                    New_Location = Movement_Function(Current_Coordinate)

                    State_New_Location = Retrieve_State_Monomer(New_Location)

                        if State_New_Location == 0 
                            Appending_Location_Possible_Coordinate_Dict(Current_Coordinate, New_Location)
                        else
                            Threads.atomic_or!(dictionary_emptied, true) 
                 
                        end 
                    break
                end
            end
        end
    end

    # Perform aggregate movement only if the dictionary was not emptied
    if dictionary_emptied[]
        Empty_Possible_Coordinates_Movement_Dict()
    else
        Move_Aggregate(Monomer)
    end
    
end


"""
Gather_all_Aggregate_Monomers_Aggregate(Monomer, Type_of_Movement)

Gathers all monomers that belong to the same Aggregate (Fibril) and checks if they can move.

# Arguments
- `Monomer::Tuple{Float64, Float64, Float64}`: The coordinates of a monomer in the Aggregate.
- `Type_of_Movement::String`: The type of movement being attempted.

# Global variables modified
- `Possible_Coordinate_Movements_Dict`

# Global variables read
- `Locations_and_States_Dict`
- `MovementFunctions`

# Calls
- `Empty_Possible_Coordinates_Movement_Dict`
- `Retrieve_Unique_Number_Monomer`
- `Filter_Aggregate`
- `Retrieve_X_Coordinate`
- `Retrieve_Y_Coordinate`
- `Retrieve_Z_Coordinate`
- Movement functions (from `MovementFunctions`)
- `Retrieve_State_Monomer`
- `lock`
- `Appending_Location_Possible_Coordinate_Dict`
- `Move_Aggregate`
"""

function Gather_all_Aggregate_Monomers_Aggregate(Monomer, Type_of_Movement)
    Empty_Possible_Coordinates_Movement_Dict()
    Unique_Code_Monomer = Retrieve_Unique_Number_Monomer(Monomer)
    Movement_Options = MovementFunctions[Type_of_Movement]
    # Filter the dictionary to get only entries with state == 4
    relevant_keys_array = Filter_Aggregate()
    dictionary_emptied = Threads.Atomic{Bool}(false)  # Flag to check if dictionary has been emptied

    # Now iterate over the filtered collection with relevant keys only
    @threads for i in 1:length(relevant_keys_array)
        Current_Coordinate = relevant_keys_array[i]
        Unique_Code_Current_Coordinate = Retrieve_Unique_Number_Monomer(Current_Coordinate)

        if Unique_Code_Monomer == Unique_Code_Current_Coordinate
            X_Coordinate, Y_Coordinate, Z_Coordinate = Retrieve_X_Coordinate(Current_Coordinate), Retrieve_Y_Coordinate(Current_Coordinate), Retrieve_Z_Coordinate(Current_Coordinate)

            for (Condition, Movement_Function) in Movement_Options
                if Condition(X_Coordinate, Y_Coordinate, Z_Coordinate)
                    New_Location = Movement_Function(Current_Coordinate)

                    State_New_Location = Retrieve_State_Monomer(New_Location)

                        if State_New_Location == 0
                            Appending_Location_Possible_Coordinate_Dict(Current_Coordinate, New_Location)
                        else
                            Threads.atomic_or!(dictionary_emptied, true)  # Set flag to true when dictionary is emptied
                            
                        end
                    break
                end
            end
        end
    end

    if dictionary_emptied[]
        Empty_Possible_Coordinates_Movement_Dict()
    else
        Move_Aggregate(Monomer)
    end
    
end

"""
Filter_Oligomer()

Filters the Locations_and_States_Dict to return an array of coordinates of Oligomers (state == 3).

# Global variables read
- `Locations_and_States_Dict`

# Returns
- `Vector{Tuple{Float64, Float64, Float64}}`: An array of coordinates of Oligomers.
"""

function Filter_Oligomer()
    return [k for (k, (state, _)) in Locations_and_States_Dict if state == 3]
end

"""
Filter_Aggregate()

Filters the Locations_and_States_Dict to return an array of coordinates of Aggregates (state == 4).

# Global variables read
- `Locations_and_States_Dict`

# Returns
- `Vector{Tuple{Float64, Float64, Float64}}`: An array of coordinates of Aggregates.
"""

function Filter_Aggregate()
    return [k for (k, (state, _)) in Locations_and_States_Dict if state == 4]
end

"""
Filter_Native()

Filters the Locations_and_States_Dict to return an array of coordinates of Native monomers (state == 1).

# Global variables read
- `Locations_and_States_Dict`

# Returns
- `Vector{Tuple{Float64, Float64, Float64}}`: An array of coordinates of Native monomers.
"""

function Filter_Native()
    return [k for (k, (state, _)) in Locations_and_States_Dict if state == 1]
end

"""
Filter_AggregateProne()

Filters the Locations_and_States_Dict to return an array of coordinates of AggregateProne monomers (state == 2).

# Global variables read
- `Locations_and_States_Dict`

# Returns
- `Vector{Tuple{Float64, Float64, Float64}}`: An array of coordinates of AggregateProne monomers.
"""

function Filter_AggregateProne()
    return [k for (k, (state, _)) in Locations_and_States_Dict if state == 2]
end


"""
Filter_Monomers()

Filters the Locations_and_States_Dict to return an array of coordinates of Monomers (state == 1 or 2).

# Global variables read
- `Locations_and_States_Dict`

# Returns
- `Vector{Tuple{Float64, Float64, Float64}}`: An array of coordinates of Monomers.
"""

function Filter_Monomers()
    return [k for (k, (state, _)) in Locations_and_States_Dict if state == 1 || state == 2]
end

"""
Empty_Possible_Coordinates_Movement_Dict()

Empties the Possible_Coordinate_Movements_Dict.

# Global variables modified
- `Possible_Coordinate_Movements_Dict`
"""

function Empty_Possible_Coordinates_Movement_Dict()
    empty!(Possible_Coordinate_Movements_Dict)
end

"""
Appending_Location_Possible_Coordinate_Dict(Current_Coordinate, New_Location)

Appends a possible movement to the Possible_Coordinate_Movements_Dict in a thread-safe manner using a lock.

# Arguments
- `Current_Coordinate::Tuple{Float64, Float64, Float64}`: The current coordinates of the monomer.
- `New_Location::Tuple{Float64, Float64, Float64}`: The new target coordinates for the movement.

# Global variables modified
- `Possible_Coordinate_Movements_Dict`

# Global variables read
- `dict_lock`: A ReentrantLock used to ensure safe concurrent writes.
"""

function Appending_Location_Possible_Coordinate_Dict(Current_Coordinate, New_Location)
    lock(dict_lock) do
    Possible_Coordinate_Movements_Dict[Current_Coordinate] = New_Location
    end
end

"""
Move_Aggregate(Monomer)

Moves an aggregate (Oligomer or Fibril) to a new location.

# Arguments
- `Monomer::Tuple{Float64, Float64, Float64}`: The coordinates of a monomer in the aggregate.

# Global variables modified
- `Locations_and_States_Dict`
- `Possible_Coordinate_Movements_Dict`

# Calls
- `Retrieve_Unique_Number_Monomer`
- `Retrieve_State_Monomer`
- `Update_Locations_States` (multiple times)
- `Empty_Possible_Coordinates_Movement_Dict`
"""

function Move_Aggregate(Monomer)
    Unique_Code_Monomer = Retrieve_Unique_Number_Monomer(Monomer)
    State_Monomer = Retrieve_State_Monomer(Monomer)
    possible_keys = collect(keys(Possible_Coordinate_Movements_Dict))


    @threads for i in 1:length(possible_keys)
        #for i in 1:length(possible_keys)
        Old_Position = possible_keys[i]
        New_Position = Possible_Coordinate_Movements_Dict[Old_Position]


        # Set old position to empty
        Update_Locations_States(Old_Position, 0, 0)


        # Set new position to the monomer's state
        Update_Locations_States(New_Position, State_Monomer, Unique_Code_Monomer)
    end
    Empty_Possible_Coordinates_Movement_Dict()
end

"""
Update_Locations_States(Position, State, Unique_Code)

Updates the state and unique code of a location in the lattice.

# Arguments
- `Position::Tuple{Float64, Float64, Float64}`: The coordinates of the location.
- `State::Int64`: The new state of the location.
- `Unique_Code::Int64`: The new unique code of the location.

# Global variables modified
- `Locations_and_States_Dict`
- `dict_lock` (for thread safety)

# Calls
- `lock`
"""

function Update_Locations_States(Position, State, Unique_Code)
    # Lock to ensure only one thread writes to the dictionary at a time
    lock(dict_lock) do
        Locations_and_States_Dict[Position] = (State, Unique_Code)
    end
end

"""
Retrieve_State_Monomer(Monomer)

Retrieves the state of a monomer.

# Arguments
- `Monomer::Tuple{Float64, Float64, Float64}`: The coordinates of the monomer.

# Global variables read
- `Locations_and_States_Dict`

# Returns
- `Int64`: The state of the monomer.
"""

function Retrieve_State_Monomer(Monomer)
    State_Monomer, _ = Locations_and_States_Dict[Monomer]
    return State_Monomer
end

"""
Count_Total_Monomers()

Counts the total number of monomers in the lattice.

# Global variables read
- `Locations_and_States_Dict`

# Returns
- `Int64`: The total number of monomers.
"""

function Count_Total_Monomers()
    total_monomers = 0

    # Iterate through all entries in Locations_and_States_Dict
    for (_, (state, _)) in Locations_and_States_Dict
        # Count only the entries where the state is non-zero (indicating a monomer exists)
        if state != 0
            total_monomers += 1
        end
    end

    return total_monomers
end

"""
Count_Monomers_With_Unique_Number(unique_number)

Counts the number of monomers with a given unique number.

# Arguments
- `unique_number::Int64`: The unique number to count.

# Global variables read
- `Locations_and_States_Dict`

# Returns
- `Int64`: The number of monomers with the given unique number.
"""

function Count_Monomers_With_Unique_Number(unique_number)
    count = 0
    
    # Iterate through all entries in Locations_and_States_Dict
    for (_, (_, unique_code)) in Locations_and_States_Dict
        # Count only entries where the unique code matches the specified unique_number
        if unique_code == unique_number
            count += 1
        end
    end

    return count
end

"""
compute_MSD()

Computes the Mean Squared Displacement (MSD) of monomers.

# Global variables read
- `Initial_Locations_and_States_Dict`
- `Locations_and_States_Dict`
- `Lattice_Size`

# Returns
- `Float64`: The Mean Squared Displacement of monomers.
"""

function compute_MSD()
    global Initial_Locations_and_States_Dict, Locations_and_States_Dict

    # Create lookup dictionaries for unique_id → position mapping
    initial_positions = Dict(uid => pos for (pos, (state, uid)) in Initial_Locations_and_States_Dict if state == 1 || state == 2)
    current_positions = Dict(uid => pos for (pos, (state, uid)) in Locations_and_States_Dict if state == 1 || state == 2)

    num_monomers = Threads.Atomic{Int}(0)
    total_displacement = Threads.Atomic{Float64}(0.0)

    # Extract all unique monomer IDs
    monomer_ids = collect(keys(initial_positions))

    # Parallel loop over unique monomer IDs
    @threads for i in 1:length(monomer_ids)
        unique_id = monomer_ids[i]

        # Retrieve initial and current positions
         # Get the current position from the lookup dictionary
         current_position = get(current_positions, unique_id, nothing)

         if current_position === nothing
             continue
         end
 
         initial_position = get(initial_positions, unique_id, nothing)
         if initial_position === nothing
             continue
         end

        # Compute displacement considering periodic boundaries
        dx = min(abs(current_position[1] - initial_position[1]), Lattice_Size - abs(current_position[1] - initial_position[1]))
        dy = min(abs(current_position[2] - initial_position[2]), Lattice_Size - abs(current_position[2] - initial_position[2]))
        dz = min(abs(current_position[3] - initial_position[3]), Lattice_Size - abs(current_position[3] - initial_position[3]))

        squared_displacement = dx^2 + dy^2 + dz^2

        # Atomically update displacement and count
        Threads.atomic_add!(total_displacement, squared_displacement)
        Threads.atomic_add!(num_monomers, 1)
    end

    return total_displacement[] / max(num_monomers[], 1)  # Avoid division by zero
end

"""
Filter_Initial_Aggregates()

Filters the Initial_Locations_and_States_Dict to return an array of coordinates of initial Aggregates (Oligomers and Fibrils).

# Global variables read
- `Initial_Locations_and_States_Dict`

# Returns
- `Vector{Tuple{Float64, Float64, Float64}}`: An array of coordinates of initial Aggregates.
"""

function Filter_Initial_Aggregates()
    return [k for (k, (state, _)) in Initial_Locations_and_States_Dict if state == 3 || state == 4]
end

"""
Filter_Current_Aggregates()

Filters the Locations_and_States_Dict to return an array of coordinates of current Aggregates (Oligomers and Fibrils).

# Calls
- `Filter_Oligomer`
- `Filter_Aggregate`

# Returns
- `Vector{Tuple{Float64, Float64, Float64}}`: An array of coordinates of current Aggregates.
"""

function Filter_Current_Aggregates()
    return vcat(Filter_Oligomer(), Filter_Aggregate())  # Merges both lists
end

"""
Group_Aggregates_By_Unique_Number()

Groups the coordinates of current Aggregates by their unique number.

# Global variables read
- `Locations_and_States_Dict`

# Calls
- `Filter_Current_Aggregates`

# Returns
- `Dict{Int64, Vector{Tuple{Float64, Float64, Float64}}}`: A dictionary mapping unique numbers to lists of coordinates.
"""

function Group_Aggregates_By_Unique_Number()
    aggregate_groups = Dict{Int, Vector{Tuple{Float64, Float64, Float64}}}()
    
    for coord in Filter_Current_Aggregates()
        state, unique_id = Locations_and_States_Dict[coord]
        
        if haskey(aggregate_groups, unique_id)
            push!(aggregate_groups[unique_id], coord)
        else
            aggregate_groups[unique_id] = [coord]
        end
    end

    return aggregate_groups  # Returns a dictionary of {unique_id => list of coordinates}
end

"""
compute_MSD_aggregates()

Computes the Mean Squared Displacement (MSD) of aggregates.

# Global variables read
- `Initial_Locations_and_States_Dict`
- `Locations_and_States_Dict`
- `Lattice_Size`

# Calls
- `Filter_Initial_Aggregates`
- `Group_Aggregates_By_Unique_Number`

# Returns
- `Float64`: The Mean Squared Displacement of aggregates.
"""

function compute_MSD_aggregates()
    # Step 1: Get initial centers of mass for oligomers and aggregates (by unique number)
    initial_coms = Dict{Int, Tuple{Float64, Float64, Float64}}()

    for coord in Filter_Initial_Aggregates()
        if haskey(Initial_Locations_and_States_Dict, coord)  # Ensure key exists
            _, unique_id = Initial_Locations_and_States_Dict[coord]
            initial_coms[unique_id] = coord  # Store initial CoM
        end
    end

    # Step 2: Compute current centers of mass grouped by unique number
    current_coms = Dict{Int, Tuple{Float64, Float64, Float64}}()
    aggregate_groups = Group_Aggregates_By_Unique_Number()  # Get {unique_id => [coords]}

    for (unique_id, coords) in aggregate_groups
        num_coords = length(coords)
        center_x = sum(c[1] for c in coords) / num_coords
        center_y = sum(c[2] for c in coords) / num_coords
        center_z = sum(c[3] for c in coords) / num_coords
        current_coms[unique_id] = (center_x, center_y, center_z)
    end

    # Debugging Step: Compare Initial and Current Aggregates
    initial_ids = Set(keys(initial_coms))
    current_ids = Set(keys(current_coms))

    missing_ids = setdiff(initial_ids, current_ids)  # Aggregates that disappeared
    new_ids = setdiff(current_ids, initial_ids)      # New aggregates that were not in initial

    # Step 3: Compute MSD based on unique ID matching
    total_displacement = 0.0
    num_aggregates = 0

    for (unique_id, current_CoM) in current_coms
        if !haskey(initial_coms, unique_id)  # Skip if initial CoM is missing
            continue
        end
        
        initial_CoM = initial_coms[unique_id]

        # Step 4: Compute displacement considering periodic boundaries
        dx = min(abs(current_CoM[1] - initial_CoM[1]), Lattice_Size - abs(current_CoM[1] - initial_CoM[1]))
        dy = min(abs(current_CoM[2] - initial_CoM[2]), Lattice_Size - abs(current_CoM[2] - initial_CoM[2]))
        dz = min(abs(current_CoM[3] - initial_CoM[3]), Lattice_Size - abs(current_CoM[3] - initial_CoM[3]))

        squared_displacement = dx^2 + dy^2 + dz^2
        total_displacement += squared_displacement
        num_aggregates += 1
    end

    return total_displacement / max(num_aggregates, 1)  # Avoid division by zero
end






"""
reset_simulation!()

Resets counters and in-memory data structures so a new run can be executed in the same Julia session.
"""
function reset_simulation!()
    # Agent/lattice state
    reset_agents!()

    # Counters
    global timesteps = 0
    global CurrentTimeStep = 0
    Total_Cleared_Monomers[] = 0

    # Data collection tables
    global results_df       = DataFrame(Timestep = Int[], Oligomers = Int[], Aggregates = Int[])
    global results_df_two   = DataFrame(Timestep = Int[], Native = Int[], AggregateProne = Int[])
    global msd_data         = DataFrame(Timestep = Int[], MSD_Monomer = Float64[], MSD_Aggregate = Float64[])
    global Number_Monomers_Cleared = DataFrame(Timestep = Int[], Number_Monomers_Cleared = Int[])

    # Movement bookkeeping
    empty!(Possible_Coordinate_Movements_Dict)
    global available_numbers = collect(1:200000)

    return nothing
end

"""
run_simulation(; output_root::Union{Nothing,AbstractString}=nothing,
               run_id::AbstractString=safe_timestamp())

Resets the simulation state, loads parameters, applies them, initializes
the lattice, creates the output directory, and runs the full simulation.

# Keyword Arguments
- `output_root::Union{Nothing,AbstractString}=nothing`: Optional output root directory.
- `run_id::AbstractString=safe_timestamp()`: Identifier used for the simulation folder.

# Returns
- `String`: Absolute path to the simulation output directory.

# Calls
- `reset_simulation!`
- `load_csv_parameters`
- `apply_parameters!`
- `initialize_simulation!`
- `Make_Directory`
- `Movement`
"""
function run_simulation(; output_root::Union{Nothing,AbstractString}=nothing,
                        run_id::AbstractString=safe_timestamp())

    reset_simulation!()

    global Parameters = load_csv_parameters(PARAMETER_FILE)

    apply_parameters!()

    initialize_simulation!()
    Make_Directory(output_root=output_root, run_id=run_id, parameter_file=PARAMETER_FILE)
    Movement()

    return directory
end