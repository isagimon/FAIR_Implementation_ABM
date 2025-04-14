using Random
using Plots
using DataFrames
using CSV
using Dates
using XLSX
using Profile
using Base.Threads


include("Agents.jl")

#INPUT VALUES
MAX_NumberMovements = 5000
Native_to_Amyloid = 0.2
Amyloid_to_Native = 0.2
Oligomer_Formation = 0.05
Oligomer_Dissociation_rate = 0.005
Fibril_Formation = 0.01
Fibril_Growth = 0.9
Fibril_No_Growth = 0.0

#################

global timesteps = 0
global CurrentTimeStep = 0
global Fibril_Length_Count
Possible_Movement_Options = ["One", "Two", "Three", "Four", "Five", "Six", "Seven", "Eight", "Nine", "Ten", "Eleven", "Twelve" , "Thirteen", "Fourteen", "Fifteen" , "Sixteen", "Seventeen", "Eighteen" , "None"] 
Possible_Coordinate_Movements_Dict = Dict{Tuple{Float64, Float64, Float64}, Tuple{Float64, Float64, Float64}}()
results_df = DataFrame(Timestep = Int[], Oligomers = Int[], Aggregates = Int[])
results_df_two = DataFrame(Timestep = Int[], Native = Int[], Amyloid = Int[])
msd_data = DataFrame(Timestep = Int[], MSD_Monomer = Float64[], MSD_Aggregate = Float64[])
available_numbers = collect(1:200000)
global max_fibril_size = Max_NumberMonomers_Amyloid + Max_NumberMonomers_Native


# Create a global lock
const dict_lock = ReentrantLock()



current_time = now()
hour_now = hour(current_time)
minute_now = minute(current_time)
second_now = second(current_time)
millisecond_now = Dates.millisecond(current_time)
am_pm = ifelse(hour_now < 12, "AM", "PM")
adjusted_hour = ifelse(hour_now > 12, hour_now - 12, ifelse(hour_now == 0, 12, hour_now))

#global timestamp = string(month(current_time), ":", day(current_time), ":", year(current_time), " Time ", adjusted_hour, "-", minute_now, "-", second_now, ".", millisecond_now, " ", am_pm)
global timestamp = string(simulation)

global use_windows_dir = false #true  # switch for conner to use his (windows) file system (set to true)

function Make_Directory()
    if use_windows_dir 
        global directory = "C:/Users/sande/OneDrive/Protein Aggregation/FCC_AmorphousAggr_V1_RESULTS/Simulation_Test"
    else
        # Use mkpath to create the entire directory path, including parent directories
        global directory = "/afs/crc.nd.edu/user/i/igimon/FCC_AmorphousAggr_V1_CRC/Data_Collection/Simulation_$timestamp"
        #println("The number of threads used: ",Threads.nthreads())
        
        # Create subdirectories
        mkpath("$directory/Bar_Graphs")
        mkpath("$directory/Heat_Maps") 
        mkpath("$directory/Line_Graphs")
        mkpath("$directory/Analysis")
    end
    # Call any additional functions if needed
    Input_Parameters()
end


function Input_Parameters()
    if use_windows_dir 
        File_Path = "C:/Users/sande/OneDrive/Protein Aggregation/FCC_AmorphousAggr_V1_RESULTS/Simulation_Information.csv" 
    else
        File_Path = "/afs/crc.nd.edu/user/i/igimon/FCC_AmorphousAggr_V1_CRC/Data_Collection/Simulation_Information.csv"  
    end
    
    Data = DataFrame(
    File_Name = [timestamp], Lattice_Size = [Lattice_Size::Int], Number_Native_Monomers = [Max_NumberMonomers_Native::Int],
    Number_Amyloid_Monomers = [Max_NumberMonomers_Amyloid::Int], Timesteps = [MAX_NumberMovements::Int], Probability_Native_to_Amyloid = [Native_to_Amyloid::Float64],
    Probability_Amyloid_to_Native = [Amyloid_to_Native::Float64], Probability_Oligomer_Formation = [Oligomer_Formation::Float64], Probability_Oligomer_Dissociation = [Oligomer_Dissociation_rate::Float64],
    Probability_Fibril_Formation = [Fibril_Formation::Float64], Probability_Fibril_Growth = [Fibril_Growth::Float64], Probability_Fibril_No_Growth = [Fibril_No_Growth::Float64]
)
    Append_Input_Parameters(File_Path, Data)
end

function Append_Input_Parameters(File_Path, Data)

    existing_data = CSV.File(File_Path) |> DataFrame
    
    appended_data = vcat(existing_data, Data)
    
    CSV.write(File_Path, appended_data)
end

function Counting_Timesteps()
    global timesteps = timesteps + 1
    println("----------------------------------")
    println("This is what timesteps looks like: ",timesteps)
    return timesteps
end
function Randomly_Chooses_Monomer()

    while true
        # Choose a random key (coordinate) from the dictionary
        RandomLocation = rand(keys(Locations_and_States_Dict))
        State, _ = Locations_and_States_Dict[RandomLocation]
        
        # Check if the state is non-zero (indicating it's a monomer)
        if State != 0
            return RandomLocation
        else
            # println("inside randomly_chooses_monomer else statement")
        end
    end
end


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

function Intial_Conditions()
    global results_df_two
    global results_df
    push!(results_df_two, (0,Max_NumberMonomers_Native, Max_NumberMonomers_Amyloid ))
    push!(results_df, (0, 0, 0))
end

function Save_MSD_Data(timesteps)
    global msd_data
    msd_calculated_monomers = compute_MSD() #This only meant for native or amyloid monomers
    msd_calculated_aggregates = compute_MSD_aggregates()
    push!(msd_data, (timesteps, msd_calculated_monomers, msd_calculated_aggregates))
end

function Export_MSD_Data()
    global directory = "/afs/crc.nd.edu/user/i/igimon/FCC_AmorphousAggr_V1_CRC/Data_Collection/Simulation_$timestamp"
    
    # Write the CSV file
    file_path = "$directory/MSD_Data.csv"
    CSV.write(file_path, msd_data)
end


function Randomly_Chooses_Movement() 
    Move = Possible_Movement_Options[rand(1:13)]
    #println("This is Move: ",Move)
    return Move
end


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


function Current_TimeStep()
    global CurrentTimeStep = CurrentTimeStep + 1
    return CurrentTimeStep
end


function Export_DataFrame(df, CurrentTimestep)

    if use_windows_dir 
        File_Path = "C:/Users/sande/OneDrive/Protein Aggregation/FCC_AmorphousAggr_V1_RESULTS/Simulation_Information.csv" 
        if CurrentTimeStep <10
            CSV.write("C:/Users/sande/OneDrive/Protein Aggregation/FCC_AmorphousAggr_V1_RESULTS/Simulation_Test/Timestep00$CurrentTimeStep.csv", df)
        elseif CurrentTimeStep >=10 && CurrentTimeStep < 100
            CSV.write("C:/Users/sande/OneDrive/Protein Aggregation/FCC_AmorphousAggr_V1_RESULTS/Simulation_Test/Timestep0$CurrentTimeStep.csv", df)
        elseif CurrentTimeStep >=100 && CurrentTimeStep < 100000
           CSV.write("C:/Users/sande/OneDrive/Protein Aggregation/FCC_AmorphousAggr_V1_RESULTS/Simulation_Test/Timestep$CurrentTimeStep.csv", df)
        end
    else
        if CurrentTimeStep <10
            CSV.write("/afs/crc.nd.edu/user/i/igimon/FCC_AmorphousAggr_V1_CRC/Data_Collection/Simulation_$timestamp/Timestep00$CurrentTimeStep.csv", df)
        elseif CurrentTimeStep >=10 && CurrentTimeStep < 100
            CSV.write("/afs/crc.nd.edu/user/i/igimon/FCC_AmorphousAggr_V1_CRC/Data_Collection/Simulation_$timestamp/Timestep0$CurrentTimeStep.csv", df)
        elseif CurrentTimeStep >=100 && CurrentTimeStep < 100000
           CSV.write("/afs/crc.nd.edu/user/i/igimon/FCC_AmorphousAggr_V1_CRC/Data_Collection/Simulation_$timestamp/Timestep$CurrentTimeStep.csv", df)
        end
    end

    
end


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

function Count_Native_Amyloid()
    native_count = 0
    amyloid_count = 0

    # Loop through the dictionary to count monomers in state 3 and 4
    for (_, (state, _)) in Locations_and_States_Dict
        if state == 1
            native_count += 1
        elseif state == 2
            amyloid_count += 1
        end
    end

    return native_count, amyloid_count

end

function Save_Data(timestep, oligomers, aggregates)
    # Append the results to the DataFrame in memory
    global results_df
    push!(results_df, (timestep, oligomers, aggregates))
end

function Save_Data_Two(timesteps, native, amyloid)
    global results_df_two
    push!(results_df_two, (timesteps, native, amyloid))
end

function Create_Fibril_Length_DataFrame()
    # Get the data and column names from the formatter
    data, column_names = Format_Fibril_Length_DataFrame()

    # Initialize the global Fibril_Length_Count DataFrame
    global Fibril_Length_Count = DataFrame(data, Symbol.(column_names))
end

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

function Export_Final_Results()
    global directory = "/afs/crc.nd.edu/user/i/igimon/FCC_AmorphousAggr_V1_CRC/Data_Collection/Simulation_$timestamp"
    
    # Write the CSV file
    file_path = "$directory/Simulation_Results.csv"
    CSV.write(file_path, results_df)
end

function Export_Final_Results_Two()
    global directory = "/afs/crc.nd.edu/user/i/igimon/FCC_AmorphousAggr_V1_CRC/Data_Collection/Simulation_$timestamp"
    
    # Write the CSV file
    file_path = "$directory/Native_and_Amyloid_Count_Results.csv"
    CSV.write(file_path, results_df_two)
end

function Export_Fibril_Length_Count()
    #global Directory = directory * "Simulation_$timestamp"
    file_path = "$directory/Fibril_Length_Count_Results.csv"
    CSV.write(file_path, Fibril_Length_Count)
end

function Total_Monomers()
    Total_Monomers = Max_NumberMonomers_Amyloid + Max_NumberMonomers_Native
    return Total_Monomers
end

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

function Update_Fibril_Length_DataFrame(timestep, fibril_lengths, max_fibril_size)
    global Fibril_Length_Count

    # Iterate over all possible fibril sizes (1 to max_fibril_size)
    for size in 1:max_fibril_size
        column_name = Symbol(string(size))
        # Set the count for this size at the current timestep
        Fibril_Length_Count[timestep, column_name] = get(fibril_lengths, size, 0)
    end
end

function Movement()
    Create_Fibril_Length_DataFrame()
    Intial_Conditions()
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

            #println("Processing Monomer: $Monomer with State: $state")
            
            # Randomly choose a movement type (e.g., "One" to "Eighteen")
            Movement = Randomly_Chooses_Movement()

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
          oligomers, aggregates = Count_Oligomers_Aggregates()
          native, amyloid = Count_Native_Amyloid()
          
          # Save data for this timestep
          Save_Data(timesteps, oligomers, aggregates)
          Save_Data_Two(timesteps, native, amyloid)
          Save_MSD_Data(timesteps)
          Count_Fibril_Length(timesteps)
          #println("This is Locations_and_States_Dict: ",Locations_and_States_Dict)

        # Timing for the timestep
        Past_Time_Raw = Current_Time_Raw
        Current_Time_Raw = now()
        Time_Taken_For_This_Timestep = Current_Time_Raw - Past_Time_Raw
        println("Time taken for timestep $timesteps: $Time_Taken_For_This_Timestep")
        
    end
    
    # Export results after all timesteps
    Export_Final_Results()
    Export_Final_Results_Two()
    Export_MSD_Data()
    Export_Fibril_Length_Count()
end

function Retrieve_X_Coordinate(Monomer)
    return Monomer[1]  # Return the X coordinate
end

function Retrieve_Y_Coordinate(Monomer)
    return Monomer[2]  # Return the Y coordinate
end

function Retrieve_Z_Coordinate(Monomer)
    return Monomer[3]  # Return the Z coordinate
end

function Movement_One_Coordinate(Monomer) 
    X_Coordinate = Retrieve_X_Coordinate(Monomer) + 1
    Y_Coordinate = Retrieve_Y_Coordinate(Monomer)
    Z_Coordinate = Retrieve_Z_Coordinate(Monomer) 
    return X_Coordinate , Y_Coordinate , Z_Coordinate
end

function Movement_One_Coordinate_Exception(Monomer) 
    X_Coordinate = 0
    Y_Coordinate = Retrieve_Y_Coordinate(Monomer)
    Z_Coordinate = Retrieve_Z_Coordinate(Monomer) 
    return X_Coordinate , Y_Coordinate , Z_Coordinate
end

function Movement_One_Coordinate_Exception_Second(Monomer) 
    X_Coordinate = (Retrieve_X_Coordinate(Monomer) + 1 ) - Lattice_Size
    Y_Coordinate = Retrieve_Y_Coordinate(Monomer)
    Z_Coordinate = Retrieve_Z_Coordinate(Monomer) 
    return X_Coordinate , Y_Coordinate , Z_Coordinate
end


function Movement_Two_Coordinate(Monomer) 
    X_Coordinate = Retrieve_X_Coordinate(Monomer)
    Y_Coordinate = Retrieve_Y_Coordinate(Monomer) + 1
    Z_Coordinate = Retrieve_Z_Coordinate(Monomer) 
    return X_Coordinate , Y_Coordinate , Z_Coordinate
end

function Movement_Two_Coordinate_Exception(Monomer) 
    X_Coordinate = Retrieve_X_Coordinate(Monomer)
    Y_Coordinate = 0
    Z_Coordinate = Retrieve_Z_Coordinate(Monomer) 
    return X_Coordinate , Y_Coordinate , Z_Coordinate
end


function Movement_Two_Coordinate_Exception_Second(Monomer) 
    X_Coordinate = Retrieve_X_Coordinate(Monomer)
    Y_Coordinate = (Retrieve_Y_Coordinate(Monomer) + 1 ) - Lattice_Size
    Z_Coordinate = Retrieve_Z_Coordinate(Monomer) 
    return X_Coordinate , Y_Coordinate , Z_Coordinate
end


function Movement_Three_Coordinate(Monomer) 
    X_Coordinate = Retrieve_X_Coordinate(Monomer)
    Y_Coordinate = Retrieve_Y_Coordinate(Monomer)
    Z_Coordinate = Retrieve_Z_Coordinate(Monomer) + 1
    return X_Coordinate , Y_Coordinate , Z_Coordinate
end

function Movement_Three_Coordinate_Exception(Monomer) 
    X_Coordinate = Retrieve_X_Coordinate(Monomer)
    Y_Coordinate = Retrieve_Y_Coordinate(Monomer)
    Z_Coordinate = 0
    return X_Coordinate , Y_Coordinate , Z_Coordinate
end

function Movement_Three_Coordinate_Exception_Second(Monomer) 
    X_Coordinate = Retrieve_X_Coordinate(Monomer)
    Y_Coordinate = Retrieve_Y_Coordinate(Monomer)
    Z_Coordinate = (Retrieve_Z_Coordinate(Monomer) + 1 ) - Lattice_Size
    return X_Coordinate , Y_Coordinate , Z_Coordinate
end


function Movement_Four_Coordinate(Monomer) 
    X_Coordinate = Retrieve_X_Coordinate(Monomer) - 1
    Y_Coordinate = Retrieve_Y_Coordinate(Monomer)
    Z_Coordinate = Retrieve_Z_Coordinate(Monomer) 
    return X_Coordinate , Y_Coordinate , Z_Coordinate
end

function Movement_Four_Coordinate_Exception(Monomer) 
    X_Coordinate = Lattice_Size
    Y_Coordinate = Retrieve_Y_Coordinate(Monomer)
    Z_Coordinate = Retrieve_Z_Coordinate(Monomer) 
    return X_Coordinate , Y_Coordinate , Z_Coordinate
end

function Movement_Four_Coordinate_Exception_Second(Monomer) 
    X_Coordinate = Lattice_Size - .5
    Y_Coordinate = Retrieve_Y_Coordinate(Monomer)
    Z_Coordinate = Retrieve_Z_Coordinate(Monomer) 
    return X_Coordinate , Y_Coordinate , Z_Coordinate
end

function Movement_Five_Coordinate(Monomer) 
    X_Coordinate = Retrieve_X_Coordinate(Monomer)
    Y_Coordinate = Retrieve_Y_Coordinate(Monomer) - 1
    Z_Coordinate = Retrieve_Z_Coordinate(Monomer) 
    return X_Coordinate , Y_Coordinate , Z_Coordinate
end

function Movement_Five_Coordinate_Exception(Monomer) 
    X_Coordinate = Retrieve_X_Coordinate(Monomer)
    Y_Coordinate = Lattice_Size
    Z_Coordinate = Retrieve_Z_Coordinate(Monomer) 
    return X_Coordinate , Y_Coordinate , Z_Coordinate
end

function Movement_Five_Coordinate_Exception_Second(Monomer) 
    X_Coordinate = Retrieve_X_Coordinate(Monomer)
    Y_Coordinate = Lattice_Size - .5
    Z_Coordinate = Retrieve_Z_Coordinate(Monomer) 
    return X_Coordinate , Y_Coordinate , Z_Coordinate
end

function Movement_Six_Coordinate(Monomer) 
    X_Coordinate = Retrieve_X_Coordinate(Monomer)
    Y_Coordinate = Retrieve_Y_Coordinate(Monomer)
    Z_Coordinate = Retrieve_Z_Coordinate(Monomer) - 1
    return X_Coordinate , Y_Coordinate , Z_Coordinate
end

function Movement_Six_Coordinate_Exception(Monomer) 
    X_Coordinate = Retrieve_X_Coordinate(Monomer)
    Y_Coordinate = Retrieve_Y_Coordinate(Monomer)
    Z_Coordinate = Lattice_Size
    return X_Coordinate , Y_Coordinate , Z_Coordinate
end

function Movement_Six_Coordinate_Exception_Second(Monomer) 
    X_Coordinate = Retrieve_X_Coordinate(Monomer)
    Y_Coordinate = Retrieve_Y_Coordinate(Monomer)
    Z_Coordinate = Lattice_Size - .5
    return X_Coordinate , Y_Coordinate , Z_Coordinate
end


function Movement_Seven_Coordinate(Monomer) 
    X_Coordinate = Retrieve_X_Coordinate(Monomer) + .5
    Y_Coordinate = Retrieve_Y_Coordinate(Monomer) - .5
    Z_Coordinate = Retrieve_Z_Coordinate(Monomer)
    return X_Coordinate , Y_Coordinate , Z_Coordinate 
end
 
function Movement_Seven_Coordinate_Exception(Monomer) 
    X_Coordinate = Retrieve_X_Coordinate(Monomer)
    Y_Coordinate = Retrieve_Y_Coordinate(Monomer)
    Z_Coordinate = Retrieve_Z_Coordinate(Monomer)
    return X_Coordinate , Y_Coordinate , Z_Coordinate 
end

function Movement_Seven_Coordinate_Exception_Second(Monomer) 
    X_Coordinate = Retrieve_Y_Coordinate(Monomer)
    Y_Coordinate = Retrieve_X_Coordinate(Monomer)
    Z_Coordinate = Retrieve_Z_Coordinate(Monomer)
    return  X_Coordinate , Y_Coordinate, Z_Coordinate 
end

function Movement_Seven_Coordinate_Exception_Third(Monomer) 
    X_Coordinate = 0
    Y_Coordinate = Lattice_Size
    Z_Coordinate = Retrieve_Z_Coordinate(Monomer)
    return X_Coordinate , Y_Coordinate , Z_Coordinate 
end

function Movement_Eight_Coordinate(Monomer) 
    X_Coordinate = Retrieve_X_Coordinate(Monomer) - .5
    Y_Coordinate = Retrieve_Y_Coordinate(Monomer) + .5
    Z_Coordinate = Retrieve_Z_Coordinate(Monomer) 
    return X_Coordinate , Y_Coordinate , Z_Coordinate
end

function Movement_Eight_Coordinate_Exception(Monomer) 
    X_Coordinate = Retrieve_X_Coordinate(Monomer)
    Y_Coordinate = Retrieve_Y_Coordinate(Monomer)
    Z_Coordinate = Retrieve_Z_Coordinate(Monomer) 
    return X_Coordinate , Y_Coordinate , Z_Coordinate
end

function Movement_Eight_Coordinate_Exception_Second(Monomer) 
    X_Coordinate = Retrieve_Y_Coordinate(Monomer)
    Y_Coordinate = Retrieve_X_Coordinate(Monomer)
    Z_Coordinate = Retrieve_Z_Coordinate(Monomer) 
    return X_Coordinate , Y_Coordinate, Z_Coordinate
end

function Movement_Eight_Coordinate_Exception_Third(Monomer) 
    X_Coordinate = Lattice_Size
    Y_Coordinate = 0
    Z_Coordinate = Retrieve_Z_Coordinate(Monomer) 
    return Y_Coordinate , X_Coordinate, Z_Coordinate
end


function Movement_Nine_Coordinate(Monomer) 
    X_Coordinate = Retrieve_X_Coordinate(Monomer)
    Y_Coordinate = Retrieve_Y_Coordinate(Monomer) + .5
    Z_Coordinate = Retrieve_Z_Coordinate(Monomer) - .5
    return X_Coordinate , Y_Coordinate , Z_Coordinate 
end

function Movement_Nine_Coordinate_Exception(Monomer) 
    X_Coordinate = Retrieve_X_Coordinate(Monomer)
    Y_Coordinate = Retrieve_Y_Coordinate(Monomer)
    Z_Coordinate = Retrieve_Z_Coordinate(Monomer)
    return X_Coordinate , Y_Coordinate , Z_Coordinate 
end

function Movement_Nine_Coordinate_Exception_Second(Monomer) 
    X_Coordinate = Retrieve_X_Coordinate(Monomer)
    Y_Coordinate = Retrieve_Z_Coordinate(Monomer)
    Z_Coordinate = Retrieve_Y_Coordinate(Monomer)
    return X_Coordinate , Y_Coordinate, Z_Coordinate
end

function Movement_Nine_Coordinate_Exception_Third(Monomer) 
    X_Coordinate = Retrieve_X_Coordinate(Monomer)
    Y_Coordinate = 0
    Z_Coordinate = Lattice_Size
    return X_Coordinate , Z_Coordinate , Y_Coordinate
end


function Movement_Ten_Coordinate(Monomer) 
    X_Coordinate = Retrieve_X_Coordinate(Monomer)
    Y_Coordinate = Retrieve_Y_Coordinate(Monomer) - .5
    Z_Coordinate = Retrieve_Z_Coordinate(Monomer) + .5
    return X_Coordinate , Y_Coordinate , Z_Coordinate
end

function Movement_Ten_Coordinate_Exception(Monomer) 
    X_Coordinate = Retrieve_X_Coordinate(Monomer)
    Y_Coordinate = Retrieve_Y_Coordinate(Monomer)
    Z_Coordinate = Retrieve_Z_Coordinate(Monomer) 
    return X_Coordinate , Y_Coordinate , Z_Coordinate
end

function Movement_Ten_Coordinate_Exception_Second(Monomer) 
    X_Coordinate = Retrieve_X_Coordinate(Monomer)
    Y_Coordinate = Retrieve_Z_Coordinate(Monomer) 
    Z_Coordinate = Retrieve_Y_Coordinate(Monomer) 
    return X_Coordinate , Y_Coordinate , Z_Coordinate
end

function Movement_Ten_Coordinate_Exception_Third(Monomer) 
    X_Coordinate = Retrieve_X_Coordinate(Monomer)
    Y_Coordinate = Lattice_Size
    Z_Coordinate = 0
    return X_Coordinate , Y_Coordinate , Z_Coordinate
end

function Movement_Eleven_Coordinate(Monomer) 
    X_Coordinate = Retrieve_X_Coordinate(Monomer) + .5
    Y_Coordinate = Retrieve_Y_Coordinate(Monomer) 
    Z_Coordinate = Retrieve_Z_Coordinate(Monomer) - .5
    return X_Coordinate , Y_Coordinate , Z_Coordinate
end

function Movement_Eleven_Coordinate_Exception(Monomer) 
    X_Coordinate = Retrieve_X_Coordinate(Monomer)
    Y_Coordinate = Retrieve_Y_Coordinate(Monomer)
    Z_Coordinate = Retrieve_Z_Coordinate(Monomer) 
    return X_Coordinate , Y_Coordinate , Z_Coordinate
end

function Movement_Eleven_Coordinate_Exception_Second(Monomer) 
    X_Coordinate = Retrieve_Z_Coordinate(Monomer) 
    Y_Coordinate = Retrieve_Y_Coordinate(Monomer)
    Z_Coordinate = Retrieve_X_Coordinate(Monomer) 
    return X_Coordinate , Y_Coordinate , Z_Coordinate
end

function Movement_Eleven_Coordinate_Exception_Third(Monomer) 
    X_Coordinate = 0
    Y_Coordinate = Retrieve_Y_Coordinate(Monomer)
    Z_Coordinate = Lattice_Size
    return X_Coordinate , Y_Coordinate , Z_Coordinate
end

function Movement_Twelve_Coordinate(Monomer)
    X_Coordinate = Retrieve_X_Coordinate(Monomer) - .5
    Y_Coordinate = Retrieve_Y_Coordinate(Monomer)
    Z_Coordinate = Retrieve_Z_Coordinate(Monomer) + .5
    return X_Coordinate , Y_Coordinate , Z_Coordinate
end

function Movement_Twelve_Coordinate_Exception(Monomer) 
    X_Coordinate = Retrieve_X_Coordinate(Monomer)
    Y_Coordinate = Retrieve_Y_Coordinate(Monomer)
    Z_Coordinate = Retrieve_Z_Coordinate(Monomer) 
    return X_Coordinate , Y_Coordinate , Z_Coordinate
end

function Movement_Twelve_Coordinate_Exception_Second(Monomer) 
    X_Coordinate = Retrieve_Z_Coordinate(Monomer)
    Y_Coordinate = Retrieve_Y_Coordinate(Monomer)
    Z_Coordinate = Retrieve_X_Coordinate(Monomer)
    return X_Coordinate , Y_Coordinate , Z_Coordinate
end

function Movement_Twelve_Coordinate_Exception_Third(Monomer) 
    X_Coordinate = Lattice_Size
    Y_Coordinate = Retrieve_Y_Coordinate(Monomer)
    Z_Coordinate = 0
    return X_Coordinate , Y_Coordinate , Z_Coordinate
end

function Movement_Thirteen_Coordinate_Exception(Monomer) 
    X_Coordinate = Retrieve_X_Coordinate(Monomer)
    Y_Coordinate = Retrieve_Y_Coordinate(Monomer)
    Z_Coordinate = Retrieve_Z_Coordinate(Monomer)
    return X_Coordinate , Y_Coordinate , Z_Coordinate 
end

function Movement_Thirteen_Coordinate_Exception_Second(Monomer)
    X_Coordinate = Lattice_Size - Retrieve_Z_Coordinate(Monomer)
    Y_Coordinate = Retrieve_Y_Coordinate(Monomer)
    Z_Coordinate = Lattice_Size
    return X_Coordinate , Y_Coordinate , Z_Coordinate 
end

function Movement_Thirteen_Coordinate_Exception_Third(Monomer) 
    X_Coordinate = Lattice_Size
    Y_Coordinate = Retrieve_Y_Coordinate(Monomer)
    Z_Coordinate = Lattice_Size - Retrieve_X_Coordinate(Monomer)
    return X_Coordinate , Y_Coordinate , Z_Coordinate 
end

function Movement_Thirteen_Coordinate_Exception_Fourth(Monomer) 
    X_Coordinate = Lattice_Size
    Y_Coordinate = Retrieve_Y_Coordinate(Monomer)
    Z_Coordinate = Lattice_Size
    return X_Coordinate , Y_Coordinate , Z_Coordinate 
end

function Movement_Thirteen_Coordinate(Monomer) 
    X_Coordinate = Retrieve_X_Coordinate(Monomer) - .5
    Y_Coordinate = Retrieve_Y_Coordinate(Monomer)
    Z_Coordinate = Retrieve_Z_Coordinate(Monomer) -.5
    return X_Coordinate , Y_Coordinate , Z_Coordinate 
end


function Movement_Fourteen_Coordinate_Exception(Monomer)  
    X_Coordinate = Retrieve_X_Coordinate(Monomer)
    Y_Coordinate = Retrieve_Y_Coordinate(Monomer)
    Z_Coordinate = Retrieve_Z_Coordinate(Monomer)
    return X_Coordinate , Y_Coordinate , Z_Coordinate 
end

function Movement_Fourteen_Coordinate_Exception_Second(Monomer) 
    X_Coordinate = 0
    Y_Coordinate = Retrieve_Y_Coordinate(Monomer)
    Z_Coordinate = Lattice_Size - (Retrieve_X_Coordinate(Monomer))
    return X_Coordinate , Y_Coordinate , Z_Coordinate 
end

function Movement_Fourteen_Coordinate_Exception_Third(Monomer) 
    X_Coordinate = Lattice_Size - (Retrieve_Z_Coordinate(Monomer))
    Y_Coordinate = Retrieve_Y_Coordinate(Monomer)
    Z_Coordinate = 0
    return X_Coordinate , Y_Coordinate , Z_Coordinate 
end

function Movement_Fourteen_Coordinate_Exception_Fourth(Monomer) 
    X_Coordinate = Lattice_Size 
    Y_Coordinate = Retrieve_Y_Coordinate(Monomer)
    Z_Coordinate = Lattice_Size
    return X_Coordinate , Y_Coordinate , Z_Coordinate 
end


function Movement_Fourteen_Coordinate(Monomer) 
    X_Coordinate = Retrieve_X_Coordinate(Monomer) + .5
    Y_Coordinate = Retrieve_Y_Coordinate(Monomer)
    Z_Coordinate = Retrieve_Z_Coordinate(Monomer) + .5
    return X_Coordinate , Y_Coordinate , Z_Coordinate 
end


function Movement_Fifteen_Coordinate_Exception(Monomer)  
    X_Coordinate = Retrieve_X_Coordinate(Monomer)
    Y_Coordinate = Retrieve_Y_Coordinate(Monomer)
    Z_Coordinate = Retrieve_Z_Coordinate(Monomer)
    return X_Coordinate , Y_Coordinate , Z_Coordinate 
end

function Movement_Fifteen_Coordinate_Exception_Second(Monomer)  
    X_Coordinate = Retrieve_X_Coordinate(Monomer)
    Y_Coordinate = Lattice_Size - Retrieve_Z_Coordinate(Monomer)
    Z_Coordinate = Lattice_Size
    return X_Coordinate , Y_Coordinate , Z_Coordinate 
end

function Movement_Fifteen_Coordinate_Exception_Third(Monomer)  
    X_Coordinate = Retrieve_X_Coordinate(Monomer)
    Y_Coordinate = Lattice_Size
    Z_Coordinate = Lattice_Size - Retrieve_Y_Coordinate(Monomer)
    return X_Coordinate , Y_Coordinate , Z_Coordinate 
end

function Movement_Fifteen_Coordinate_Exception_Fourth(Monomer)
    X_Coordinate = Retrieve_X_Coordinate(Monomer)
    Y_Coordinate = Lattice_Size
    Z_Coordinate = Lattice_Size 
    return X_Coordinate , Y_Coordinate , Z_Coordinate 
end

function Movement_Fifteen_Coordinate(Monomer)  
    X_Coordinate = Retrieve_X_Coordinate(Monomer)
    Y_Coordinate = Retrieve_Y_Coordinate(Monomer) - .5 
    Z_Coordinate = Retrieve_Z_Coordinate(Monomer) - .5
    return X_Coordinate , Y_Coordinate , Z_Coordinate 
end

function Movement_Sixteen_Coordinate_Exception(Monomer) 
    X_Coordinate = Retrieve_X_Coordinate(Monomer)
    Y_Coordinate = Retrieve_Y_Coordinate(Monomer)
    Z_Coordinate = Retrieve_Z_Coordinate(Monomer)
    return X_Coordinate , Y_Coordinate , Z_Coordinate 
end

function Movement_Sixteen_Coordinate_Exception_Second(Monomer)  
    X_Coordinate = Retrieve_X_Coordinate(Monomer)
    Y_Coordinate = 0
    Z_Coordinate = 0
    return X_Coordinate , Y_Coordinate , Z_Coordinate 
end

function Movement_Sixteen_Coordinate_Exception_Third(Monomer)  
    X_Coordinate = Retrieve_X_Coordinate(Monomer)
    Y_Coordinate = 0
    Z_Coordinate = Lattice_Size - Retrieve_Y_Coordinate(Monomer)
    return X_Coordinate , Y_Coordinate , Z_Coordinate 
end

function Movement_Sixteen_Coordinate_Exception_Fourth(Monomer)  
    X_Coordinate = Retrieve_X_Coordinate(Monomer)
    Y_Coordinate = Lattice_Size - Retrieve_Z_Coordinate(Monomer)
    Z_Coordinate = 0
    return X_Coordinate , Y_Coordinate , Z_Coordinate 
end

function Movement_Sixteen_Coordinate(Monomer)  
    X_Coordinate = Retrieve_X_Coordinate(Monomer)
    Y_Coordinate = Retrieve_Y_Coordinate(Monomer) + .5
    Z_Coordinate = Retrieve_Z_Coordinate(Monomer) + .5
    return X_Coordinate , Y_Coordinate , Z_Coordinate 
end

function Movement_Seventeen_Coordinate_Exception(Monomer) 
    X_Coordinate = Retrieve_X_Coordinate(Monomer)
    Y_Coordinate = Retrieve_Y_Coordinate(Monomer)
    Z_Coordinate = Retrieve_Z_Coordinate(Monomer)
    return X_Coordinate , Y_Coordinate , Z_Coordinate 
end

function Movement_Seventeen_Coordinate_Exception_Second(Monomer) 
    X_Coordinate = Lattice_Size
    Y_Coordinate = Lattice_Size
    Z_Coordinate = Retrieve_Z_Coordinate(Monomer)
    return X_Coordinate , Y_Coordinate , Z_Coordinate 
end

function Movement_Seventeen_Coordinate_Exception_Third(Monomer) 
    X_Coordinate = Lattice_Size - Retrieve_Y_Coordinate(Monomer)
    Y_Coordinate = Lattice_Size
    Z_Coordinate = Retrieve_Z_Coordinate(Monomer)
    return X_Coordinate , Y_Coordinate , Z_Coordinate 
end

function Movement_Seventeen_Coordinate_Exception_Fourth(Monomer) 
    X_Coordinate = Lattice_Size 
    Y_Coordinate = Lattice_Size - Retrieve_X_Coordinate(Monomer)
    Z_Coordinate = Retrieve_Z_Coordinate(Monomer)
    return X_Coordinate , Y_Coordinate , Z_Coordinate 
end

function Movement_Seventeen_Coordinate(Monomer)  
    X_Coordinate = Retrieve_X_Coordinate(Monomer) - .5
    Y_Coordinate = Retrieve_Y_Coordinate(Monomer) - .5
    Z_Coordinate = Retrieve_Z_Coordinate(Monomer)
    return X_Coordinate , Y_Coordinate , Z_Coordinate 
end

function Movement_Eighteen_Coordinate_Exception(Monomer) 
    X_Coordinate = Retrieve_X_Coordinate(Monomer)
    Y_Coordinate = Retrieve_Y_Coordinate(Monomer)
    Z_Coordinate = Retrieve_Z_Coordinate(Monomer)
    return X_Coordinate , Y_Coordinate , Z_Coordinate 
end

function Movement_Eighteen_Coordinate_Exception_Second(Monomer) 
    X_Coordinate = 0
    Y_Coordinate = 0
    Z_Coordinate = Retrieve_Z_Coordinate(Monomer)
    return X_Coordinate , Y_Coordinate , Z_Coordinate 
end

function Movement_Eighteen_Coordinate_Exception_Third(Monomer) 
    X_Coordinate = Lattice_Size - Retrieve_Y_Coordinate(Monomer)
    Y_Coordinate = 0
    Z_Coordinate = Retrieve_Z_Coordinate(Monomer)
    return X_Coordinate , Y_Coordinate , Z_Coordinate 
end

function Movement_Eighteen_Coordinate_Exception_Fourth(Monomer) 
    X_Coordinate = 0
    Y_Coordinate = Lattice_Size - Retrieve_X_Coordinate(Monomer)
    Z_Coordinate = Retrieve_Z_Coordinate(Monomer)
    return X_Coordinate , Y_Coordinate , Z_Coordinate 
end

function Movement_Eighteen_Coordinate(Monomer) 
    X_Coordinate = Retrieve_X_Coordinate(Monomer) + .5
    Y_Coordinate = Retrieve_Y_Coordinate(Monomer) + .5
    Z_Coordinate = Retrieve_Z_Coordinate(Monomer)
    return X_Coordinate , Y_Coordinate , Z_Coordinate 
end



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

function Distinguishing_Monomers(Monomer, Desired_Location, Type_of_Movement)
    State= Retrieve_State_Monomer(Monomer)

    if State == 1
        Native_Move_or_Conformational_Change(Monomer, Desired_Location)
    elseif State == 2
        Amyloid_Aggregation_or_Movement(Monomer, Desired_Location)
    elseif State == 3
        Oligomer_Move(Monomer, Desired_Location, Type_of_Movement)
    elseif State == 4
        Fibril_Move(Monomer, Desired_Location, Type_of_Movement)
    end

end

function Key_Array_Locations_and_States()
    return collect(keys(Locations_and_States_Dict))
end
    
function Monomer_State(State)
    if State == 1
        return "Native"
    elseif State == 2
        return "Amyloid"
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

function One_Monomer_Movement(Monomer, Desired_Location)
   # Retrieve the state and unique number of the monomer once
   monomer_state, monomer_unique_number = Locations_and_States_Dict[Monomer]

   Update_Locations_States(Desired_Location, monomer_state, monomer_unique_number)

   # Set the original location's state to 0 while keeping the unique number unchanged
   Update_Locations_States(Monomer, 0, 0)
end

function Empties_Locations_and_States(Monomer)
    Locations_and_States_Dict[Monomer] = (0, Locations_and_States_Dict[Monomer][2])
end

function Conformational_Change(Monomer)
    # Retrieve the state and unique number from the dictionary using the Monomer coordinates
    State, Unique_Number = Locations_and_States_Dict[Monomer]

    if State == 1
        # Change state from Native (1) to Amyloid (2)
        Update_Locations_States(Monomer, 2, Unique_Number)
    elseif State == 2
        # Change state from Amyloid (2) to Native (1)
        Update_Locations_States(Monomer, 1, Unique_Number)
    end


end


function Native_Move_or_Conformational_Change(Monomer, Desired_Location)

    State = Retrieve_State_Monomer(Desired_Location)

    if State == 0  # Available if state is empty (0)
        One_Monomer_Movement(Monomer, Desired_Location)
        # Only allow conformational change within this block
        if Random_Dice() < Native_to_Amyloid
            Conformational_Change(Desired_Location)
        end
    else
        # Handle conformational change only if movement doesn't happen
        if Random_Dice() < Native_to_Amyloid
            #println("I am going to do a conformational change")
            Conformational_Change(Monomer)
        end
    end

end

function Amyloid_Aggregation_or_Movement(Monomer, Desired_Location)

    Status_Desired_Location = Retrieve_State_Monomer(Desired_Location)
    location_state = Monomer_State(Status_Desired_Location)

    if Status_Desired_Location == 0
        # Case 1: Move the monomer to an empty location
        One_Monomer_Movement(Monomer, Desired_Location)
        if Random_Dice() < Amyloid_to_Native
            Conformational_Change(Desired_Location)
        end
    elseif Status_Desired_Location == 2 || Status_Desired_Location == 3
        # Case 2: Handle amyloid aggregation
        Amyloid_Aggregation(Monomer, Desired_Location, location_state)
        Status_Monomer = Retrieve_State_Monomer(Monomer)
        if Status_Monomer == 2 && Random_Dice() < Amyloid_to_Native
            Conformational_Change(Monomer)
        end
    elseif Random_Dice() < Amyloid_to_Native #If the state is 1, 4, or 5
        # Case 3: Apply conformational change directly if no other conditions are met
        Conformational_Change(Monomer)
    end

end



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


function No_Repeating_Unique_Number_Locked(UniqueCode)
    # Iterate through the values of the dictionary to check for the unique code
    for (_, unique_number) in values(Locations_and_States_Dict)
        if unique_number == UniqueCode
            return true
        end
    end
    return false
end
function Amyloid_Amyloid_Lock(Monomer, Desired_Location)
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

function Calculate_Center_of_Mass(Monomer, Desired_Location)
    Center_of_Mass_X = (Monomer[1] + Desired_Location[1]) / 2
    Center_of_Mass_Y = (Monomer[2] + Desired_Location[2]) / 2
    Center_of_Mass_Z = (Monomer[3] + Desired_Location[3]) / 2
    return Center_of_Mass_X, Center_of_Mass_Y, Center_of_Mass_Z
end

function Initialize_New_Center_of_Mass(X_Coordinate, Y_Coordinate, Z_Coordinate, State_Value, Unique_Code)
    Initial_Locations_and_States_Dict[X_Coordinate, Y_Coordinate, Z_Coordinate] = (State_Value, Unique_Code)
end
function Delete_Monomer_Information_from_Initial_Locations_and_States(unique_number, state)
     # Find the coordinate that corresponds to the unique number
     monomer_to_delete = nothing

     for (coordinate, (current_state, uid)) in Initial_Locations_and_States_Dict
         if uid == unique_number && current_state == state
             monomer_to_delete = coordinate
             println("Coordinate to be deleted in Initial_Locations_and_States_Dict: ",coordinate)
             break  # Stop after finding the first match
         end
     end
 
     # If a matching coordinate is found, delete it
     if monomer_to_delete !== nothing
         delete!(Initial_Locations_and_States_Dict, monomer_to_delete)
     end
end


function Retrieve_Unique_Number_Monomer(Monomer)
     # Get the unique number associated with this Monomer (coordinate)
     _, Unique_Number = Locations_and_States_Dict[Monomer]

     return Unique_Number
end

function Amyloid_Aggregation(Monomer, Desired_Location, Status_Desired_Location)
    # Check if the status of the desired location is "Amyloid" and if conditions meet for forming an oligomer
    if Status_Desired_Location == "Amyloid" && Random_Dice() < Oligomer_Formation && Monomer != Desired_Location #If you are next to a 2
        # Perform Amyloid-Amyloid locking
        Amyloid_Amyloid_Lock(Monomer, Desired_Location)
    elseif Status_Desired_Location == "Oligomer" && Random_Dice() < Fibril_Formation #If you are next to a 3
        # Perform Amyloid-Oligomer locking
        Amyloid_Oligomer_Lock_Two(Monomer, Desired_Location)
    end
end

function Amyloid_Oligomer_Lock_Two(Monomer_Amyloid, Desired_Location_Oligomer)
    # Retrieve the unique code associated with the amyloid and oligomer
    Unique_Code_Oligomer = Retrieve_Unique_Number_Monomer(Desired_Location_Oligomer)
    Unique_Code_Amyloid = Retrieve_Unique_Number_Monomer(Monomer_Amyloid)

    Key_Array = Key_Array_Locations_and_States()

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

    # Update the state of the amyloid monomer to 4 (Fibril) with the same unique code
    Update_Locations_States(Monomer_Amyloid, 4, Unique_Code_Oligomer)

    #Delete the info of oligomer and amyloid from Initial_Locations_and_States_Dict
    Remove_Center_of_Mass_Info(Unique_Code_Oligomer)
    Remove_Center_of_Mass_Info(Unique_Code_Amyloid)

    #Calculate the new center of mass for the new fibril
    X, Y, Z = Calculate_Center_of_Mass_Fibril(Unique_Code_Oligomer)

    #Calculate new center of mass for new fibril
    Initialize_New_Center_of_Mass(X, Y, Z, 4, Unique_Code_Oligomer)
end


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


function Oligomer_Aggregate(Monomer, Desired_Location)
    # Retrieve the states of the monomer and the desired location directly from the dictionary


    State_Desired_Location = Retrieve_State_Monomer(Desired_Location) 
    State_Monomer = Retrieve_State_Monomer(Monomer) 

    # Check for the conditions for aggregation or dissociation
    if State_Desired_Location == 2 && Random_Dice() < Fibril_Formation
        Amyloid_Oligomer_Lock(Monomer, Desired_Location)
    elseif Random_Dice() < Oligomer_Dissociation_rate && State_Monomer == 3
        Unique_Number_Oligomer = Retrieve_Unique_Number_Monomer(Monomer)
        Oligomer_Dissociation(Unique_Number_Oligomer)
    end

end

function Oligomer_Dissociation(Unique_Number_Oligomer)
    # Retrieve the unique number associated with the monomer directly from the dictionary

    # Iterate through all coordinates in the dictionary
    Key_Array = Key_Array_Locations_and_States()
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

function Remove_Center_of_Mass_Info(Unique_Number)
    for (initial_CoM, (state, uid)) in Initial_Locations_and_States_Dict
        if uid == Unique_Number  # Find the corresponding oligomer CoM
            delete!(Initial_Locations_and_States_Dict, initial_CoM)
            break  # Exit once found to avoid unnecessary checks
        end
    end
end

function Restore_Dissociated_Monomers(Monomer, Unique_Number)
        Initial_Locations_and_States_Dict[Monomer] = (2, Unique_Number)
end


function Amyloid_Oligomer_Lock(Monomer_Oligomer, Desired_Location_Amyloid)
    # Retrieve the unique number associated with the oligomer 
    Unique_Code_Oligomer = Retrieve_Unique_Number_Monomer(Monomer_Oligomer)

    # Retrieve the unique number associated with the amyloid 
    Unique_Code_Amyloid = Retrieve_Unique_Number_Monomer(Desired_Location_Amyloid)

    # Iterate through all coordinates in the dictionary
    Key_Array = Key_Array_Locations_and_States()

    @threads for i in 1:length(Key_Array)
        #for i in 1:length(Key_Array)
        coordinate = Key_Array[i]
        Unique_Code = Retrieve_Unique_Number_Monomer(coordinate)

        # If the unique code matches the oligomer's unique code, change the state to 4 (Fibril)
        if Unique_Code == Unique_Code_Oligomer
            Update_Locations_States(coordinate, 4, Unique_Code_Oligomer)
        end
    end

    # Update the state of the desired location (Amyloid) to 4 (Fibril) and set its unique code to the oligomer's unique code
    Update_Locations_States(Desired_Location_Amyloid, 4, Unique_Code_Oligomer)

    #Remove Info of Oligomer and Amyloid from Initial_Locations_and_States_Dict
    Remove_Center_of_Mass_Info(Unique_Code_Oligomer)
    Remove_Center_of_Mass_Info(Unique_Code_Amyloid)

    #Calculate new center of mass for new fibril 
    X, Y, Z = Calculate_Center_of_Mass_Fibril(Unique_Code_Oligomer)

    #Update that information from Initial_Locations_and_States_Dict
    Initialize_New_Center_of_Mass(X, Y, Z, 4, Unique_Code_Oligomer)


end

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


function Fibril_Aggregate(Monomer, Desired_Location)

 # Retrieve the state of the desired location directly from the dictionary
 State_Desired_Location = Retrieve_State_Monomer(Desired_Location)

 # Check if the desired location is "Amyloid" and a random chance allows fibril growth
 if State_Desired_Location == 2 && Random_Dice() < Fibril_Growth
     Fibril_Lock(Monomer, Desired_Location)
 end

end

function Fibril_Lock(Fibril, Desired_Location_Amyloid)
    # Retrieve the unique number associated with the fibril directly from the dictionary
    Unique_Code_Fibril = Retrieve_Unique_Number_Monomer(Fibril)
    Unique_Code_Amyloid = Retrieve_Unique_Number_Monomer(Desired_Location_Amyloid)

    # Update the state of the desired location (Amyloid) to 4 (Fibril) and set its unique code
    Update_Locations_States(Desired_Location_Amyloid, 4, Unique_Code_Fibril)
    
    #Remove Info from Initial_Locations_and_States_Dict
    Remove_Center_of_Mass_Info(Unique_Code_Fibril)
    Remove_Center_of_Mass_Info(Unique_Code_Amyloid)

    #Calculate new center of mass with bigger fibril 
    X, Y, Z = Calculate_Center_of_Mass_Fibril(Unique_Code_Fibril)

    #Input that info into Initial_Locations_and_States_Dict
    Initialize_New_Center_of_Mass(X, Y, Z, 4, Unique_Code_Fibril)
end

function Random_Dice()
    Random_Dice = 1-rand()
    #println("The Random_Dice is: $Random_Dice")
    return Random_Dice
end


function Gather_all_Aggregate_Monomers_Oligomer(Monomer, Type_of_Movement)
    Empty_Possible_Coordinates_Movement_Dict()
    Unique_Code_Monomer = Retrieve_Unique_Number_Monomer(Monomer)
    Movement_Options = MovementFunctions[Type_of_Movement]
    # Filter the dictionary to get only entries with state == 3
    relevant_keys_array = Filter_Oligomer()
    dictionary_emptied = false  # Flag to check if dictionary has been emptied

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

                    lock(dict_lock) do
                        if State_New_Location == 0 
                            Appending_Location_Possible_Coordinate_Dict(Current_Coordinate, New_Location)
                        else
                            Empty_Possible_Coordinates_Movement_Dict()
                            dictionary_emptied = true  # Set flag to true when dictionary is emptied
                            return  # Exit the thread when the dictionary is emptied
                        end 
                    end
                    break
                end
            end
        end
    end

    # Perform aggregate movement only if the dictionary was not emptied
    if !dictionary_emptied
        println("This is Possible_Coordinate_Movements_Dict: $Possible_Coordinate_Movements_Dict")
        Move_Aggregate(Monomer)
    end
end


function Gather_all_Aggregate_Monomers_Aggregate(Monomer, Type_of_Movement)
    Empty_Possible_Coordinates_Movement_Dict()
    Unique_Code_Monomer = Retrieve_Unique_Number_Monomer(Monomer)
    Movement_Options = MovementFunctions[Type_of_Movement]
    # Filter the dictionary to get only entries with state == 4
    relevant_keys_array = Filter_Aggregate()
    dictionary_emptied = false  # Flag to check if dictionary has been emptied

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

                    lock(dict_lock) do
                        if State_New_Location == 0
                            Appending_Location_Possible_Coordinate_Dict(Current_Coordinate, New_Location)
                        else
                            Empty_Possible_Coordinates_Movement_Dict()
                            dictionary_emptied = true  # Set flag to true when dictionary is emptied
                            return  # Exit the thread when the dictionary is emptied
                        end
                    end
                    break
                end
            end
        end
    end

    # Perform aggregate movement only if the dictionary was not emptied
    if !dictionary_emptied
        Move_Aggregate(Monomer)
    end
end

function Filter_Oligomer()
    return [k for (k, (state, _)) in Locations_and_States_Dict if state == 3]
end

function Filter_Aggregate()
    return [k for (k, (state, _)) in Locations_and_States_Dict if state == 4]
end

function Filter_Monomers()
    return [k for (k, (state, _)) in Locations_and_States_Dict if state == 1 || state == 2]
end

function Empty_Possible_Coordinates_Movement_Dict()
    empty!(Possible_Coordinate_Movements_Dict)
end

function Appending_Location_Possible_Coordinate_Dict(Current_Coordinate, New_Location)
    Possible_Coordinate_Movements_Dict[Current_Coordinate] = New_Location
end

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
    #println("The number of monomers that make up this oligomer or aggregate before Move_Aggregate: ",Count_Monomers_With_Unique_Number(Unique_Code_Monomer))

    Empty_Possible_Coordinates_Movement_Dict()
end


function Update_Locations_States(Position, State, Unique_Code)
    # Lock to ensure only one thread writes to the dictionary at a time
    lock(dict_lock) do
        #println("A position was changed in Update_Locations_States")
        Locations_and_States_Dict[Position] = (State, Unique_Code)
    end
end



function Retrieve_State_Monomer(Monomer)
    State_Monomer, _ = Locations_and_States_Dict[Monomer]
    return State_Monomer
end

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
             println("⚠️ Warning: unique_id $unique_id not found in current_positions")
             continue
         end
 
         initial_position = get(initial_positions, unique_id, nothing)
         if initial_position === nothing
             println("⚠️ Warning: unique_id $unique_id not found in initial_positions")
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


function Filter_Initial_Aggregates()
    return [k for (k, (state, _)) in Initial_Locations_and_States_Dict if state == 3 || state == 4]
end

function Filter_Current_Aggregates()
    return vcat(Filter_Oligomer(), Filter_Aggregate())  # Merges both lists
end

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

    if !isempty(missing_ids)
        println("Warning: Missing aggregates -> ", missing_ids)
    end
    if !isempty(new_ids)
        println("Note: New aggregates detected -> ", new_ids)
    end

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




#Profile.clear()
Make_Directory()
#@profile Movement()
#Profile.print()
Movement()