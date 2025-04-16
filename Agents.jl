using Random
using Plots
using CSV
using DataFrames

###############
# Constants
###############

# State values for agents
const NativeState_Value         = 1  # Native monomer
const AmyloidProne_Value        = 2  # Amyloid-prone monomer
const OligomerState_Value       = 3  # Oligomer
const FibrilState_Value         = 4  # Fibril
const SphereState_Value         = 5  # Crowder/Sphere

# Unique number ranges
const MONOMER_ID_RANGE     = 200_001:499_999
const SPHERE_ID_RANGE      = 500_000:2_000_000

###############
# Global Simulation State
###############

# Stores the current state of all lattice locations
global Locations_and_States_Dict = Dict{Tuple{Float64, Float64, Float64}, Tuple{Int, Int}}()

# Stores possible coordinates used to assemble a spherical crowder
global Possible_Sphere_Coordinates_Set = Dict{Tuple{Float64, Float64, Float64}, Nothing}()

# Records the initial state of all lattice locations for later reference
global Initial_Locations_and_States_Dict = Dict{Tuple{Float64, Float64, Float64}, Tuple{Int, Int}}()

# Tracks centers that have already been used for sphere generation
global used_centers = Set{Tuple{Float64, Float64, Float64}}()

# Unique ID pools
global Sphere_Unique_Numbers = collect(SPHERE_ID_RANGE)
global Monomer_Unique_Numbers = collect(MONOMER_ID_RANGE)

"""
    load_csv_parameters(file_path::String) -> Dict

Loads simulation parameters from a CSV file into a dictionary.
The CSV should contain at least two columns: parameter names and values.

# Arguments
- `file_path`: Path to the CSV file.

# Returns
- A dictionary mapping parameter names (as strings) to values (parsed as Float64, Bool, or String).
"""

function load_csv_parameters(file_path::String)
    # Read the CSV file into a DataFrame
    df = CSV.read(file_path, DataFrame)

    # Standardize column names to remove hidden characters
    rename!(df, strip.(names(df)))
    println("Updated column names: ", names(df))

    # Define column names
    param_col = names(df)[1]  # First column for parameter names
    value_col = names(df)[2]  # Second column for parameter values

    # Create a dictionary from the first two columns
    params = Dict(
        strip(row[param_col]) => 
        try
            parse(Float64, strip(string(row[value_col])))
        catch
            if strip(string(row[value_col])) == "TRUE"
                true
            elseif strip(string(row[value_col])) == "FALSE"
                false
            else
                strip(string(row[value_col]))
            end
        end
        for row in eachrow(df)
    )

    # Ensure integer values for specific parameters
    for key in ["Lattice_Size", "Max_NumberMonomers_Native", "Max_NumberMonomers_Amyloid", "Obstacle_Radius", "MAX_NumberMovements"]
        if haskey(params, key)
            params[key] = Int(params[key])
        end
    end

    println("Loaded parameters: ", keys(params))
    return params
end



# Path to the input parameters CSV file
file_path = "/Users/isabellagimon/Desktop/FAIR_Implementation_ABM/Input_Parameters.csv"

# Load parameters
Parameters = load_csv_parameters(file_path)

# Assign individual variables
Lattice_Size = Parameters["Lattice_Size"]
Max_NumberMonomers_Native = Parameters["Max_NumberMonomers_Native"]
Max_NumberMonomers_Amyloid = Parameters["Max_NumberMonomers_Amyloid"]
Obstacle_Radius = Parameters["Obstacle_Radius"]
Crowder_Concentration_Spheres = Parameters["Crowder_Concentration_Spheres"]
Obstacle = Parameters["Spheres?"]
Sphere_Volume = Parameters["Sphere_Volume"]

# Debugging: Print all parameters
println("Loaded Parameters: ", Parameters)

"""
    Generate_Coordinates(Lattice_Size::Int)

Generates the 3D FCC lattice of coordinates and assigns initial states to monomers and spheres based on provided parameters.

# Arguments
- `Lattice_Size`: Integer value for the size of the lattice in each dimension.
"""


function Generate_Coordinates(Lattice_Size)
    global Locations_and_States_Dict

    for X in 0:(Lattice_Size - 1)
        for Y in 0:(Lattice_Size - 1)
            for Z in 0:(Lattice_Size - 1)
                State = 0
                Unique_Number = 0

                # Add corner positions for each unit cell
                Add_Position(X, Y, Z, State, Unique_Number)
            end
        end
    end
    if Obstacle == true
         Differentiate_Sphere_Crowder_Radius()
    end
    Randomly_Assigns_Location_Monomers_Native()
    Randomly_Assigns_Location_Monomers_Amyloid()
    Copy_Original_Location()
 
end

"""
    Copy_Original_Location()

Stores a snapshot of the initial state and positions of native and amyloid monomers in the global `Initial_Locations_and_States_Dict`.
"""

function Copy_Original_Location()
    global Initial_Locations_and_States_Dict  # Ensure it modifies the global variable

    for (location, (state, unique_number)) in Locations_and_States_Dict
        if state == NativeState_Value || state == AmyloidProne_Value
            Initial_Locations_and_States_Dict[location] = (state, unique_number)
        end
    end
end

"""
    Add_Position(X, Y, Z, State, Unique_Number)

Adds all face- and corner-centered positions for a unit cell at coordinate (X, Y, Z) into the global `Locations_and_States_Dict`.

# Arguments
- `X`, `Y`, `Z`: Integer coordinates of the unit cell origin.
- `State`: The initial state (typically 0).
- `Unique_Number`: The initial unique identifier.
"""


function Add_Position(X, Y, Z, State, Unique_Number)
    global Locations_and_States_Dict

    # Define all corner and face positions to be added to the dictionary
    positions = [
        First_Corner_Position(X, Y, Z),
        Second_Corner_Position(X, Y, Z),
        Third_Corner_Position(X, Y, Z),
        Fourth_Corner_Position(X, Y, Z),
        Fifth_Corner_Position(X, Y, Z),
        Sixth_Corner_Position(X, Y, Z),
        Seventh_Corner_Position(X, Y, Z),
        Eighth_Corner_Position(X, Y, Z),
        First_Face_Position(X, Y, Z),
        Second_Face_Position(X, Y, Z),
        Third_Face_Position(X, Y, Z),
        Fourth_Face_Position(X, Y, Z),
        Fifth_Face_Position(X, Y, Z),
        Sixth_Face_Position(X, Y, Z)
    ]

    for position in positions
        Locations_and_States_Dict[position] = (State, Unique_Number)
    end
end

"""
    First_Corner_Position(X, Y, Z) -> Tuple

Returns the coordinate of the first corner of the unit cell, which corresponds to the origin (X, Y, Z).
"""

function First_Corner_Position(X, Y, Z)
    return X, Y, Z
end

"""
    Second_Corner_Position(X, Y, Z) -> Tuple

Returns the coordinate of the second corner of the unit cell, offset by +1 in the X direction.
"""

function Second_Corner_Position(X, Y, Z)
    return X + 1, Y, Z
end

"""
    Third_Corner_Position(X, Y, Z) -> Tuple

Returns the coordinate of the third corner of the unit cell, offset by +1 in the Y direction.
"""

function Third_Corner_Position(X, Y, Z)
    return X, Y + 1, Z
end

"""
    Fourth_Corner_Position(X, Y, Z) -> Tuple

Returns the coordinate of the fourth corner of the unit cell, offset by +1 in the Z direction.
"""

function Fourth_Corner_Position(X, Y, Z)
    return X, Y, Z + 1
end

"""
    Fifth_Corner_Position(X, Y, Z) -> Tuple

Returns the coordinate of the fifth corner, offset by +1 in both X and Y directions.
"""

function Fifth_Corner_Position(X, Y, Z)
    return X + 1, Y + 1, Z
end

"""
    Sixth_Corner_Position(X, Y, Z) -> Tuple

Returns the coordinate of the sixth corner, offset by +1 in X and Z directions.
"""

function Sixth_Corner_Position(X, Y, Z)
    return X + 1, Y, Z + 1
end

"""
    Seventh_Corner_Position(X, Y, Z) -> Tuple

Returns the coordinate of the seventh corner, offset by +1 in Y and Z directions.
"""

function Seventh_Corner_Position(X, Y, Z)
    return X, Y + 1, Z + 1
end

"""
    Eighth_Corner_Position(X, Y, Z) -> Tuple

Returns the coordinate of the eighth corner, offset by +1 in X, Y, and Z directions.
"""

function Eighth_Corner_Position(X, Y, Z)
    return X + 1, Y + 1, Z + 1
end

"""
    First_Face_Position(X, Y, Z) -> Tuple

Returns the coordinate of the first face-centered position, centered between X and Y directions.
"""

function First_Face_Position(X, Y, Z) #Face-centered along each axis
    return X + 0.5, Y + .5, Z
end

"""
    Second_Face_Position(X, Y, Z) -> Tuple

Returns the coordinate of the second face-centered position, centered between Y and Z directions.
"""

function Second_Face_Position(X, Y, Z) #Face-centered along each axis
    return X, Y + 0.5, Z + .5
end

"""
    Third_Face_Position(X, Y, Z) -> Tuple

Returns the coordinate of the third face-centered position, centered between X and Z directions.
"""
 
function Third_Face_Position(X, Y, Z) #Face-centered along each axis
    return X + .5, Y, Z + 0.5
end

"""
    Fourth_Face_Position(X, Y, Z) -> Tuple

Returns the coordinate of the fourth face-centered position, located at the top face in Z.
"""

function Fourth_Face_Position(X, Y, Z) #Opposite Faces
    return X + .5, Y + 0.5, Z + 1
end

"""
    Fifth_Face_Position(X, Y, Z) -> Tuple

Returns the coordinate of the fifth face-centered position, located at the front face in Y.
"""

function Fifth_Face_Position(X, Y, Z)  #Opposite Faces
    return X + 0.5, Y + 1, Z + .5
end

"""
    Sixth_Face_Position(X, Y, Z) -> Tuple

Returns the coordinate of the sixth face-centered position, located at the right face in X.
"""

function Sixth_Face_Position(X, Y, Z)  #Opposite Faces
    return X + 1, Y + .5, Z + 0.5
end

"""
    Randomly_Assigns_Location_Monomers_Native()

Randomly assigns native monomer state (1) to available lattice positions.
The number of assignments is determined by `Max_NumberMonomers_Native`.
"""

# Function to assign states randomly to native monomers
function Randomly_Assigns_Location_Monomers_Native()
    global Locations_and_States_Dict
    println("We are in Randomly_Assigns_Location_Monomers_Native")
    Monomers_Made_Native = 0
    keys_list = collect(keys(Locations_and_States_Dict))

    while Monomers_Made_Native < Max_NumberMonomers_Native
        Random_Index = rand(1:length(keys_list))
        Random_Location = keys_list[Random_Index]
        State, _ = Locations_and_States_Dict[Random_Location]

        if State == 0
            Assigns_State_Monomer_Native(Random_Location)
            Monomers_Made_Native += 1
        end
    end
end

"""
    Randomly_Assigns_Location_Monomers_Amyloid()

Randomly assigns amyloid-prone monomer state (2) to available lattice positions.
The number of assignments is determined by `Max_NumberMonomers_Amyloid`.
"""

# Function to assign states randomly to amyloid monomers
function Randomly_Assigns_Location_Monomers_Amyloid()
    global Locations_and_States_Dict
    println("We are in Randomly_Assigns_Location_Monomers_Amyloid")
    Monomers_Made_Amyloid = 0
    keys_list = collect(keys(Locations_and_States_Dict))

    while Monomers_Made_Amyloid < Max_NumberMonomers_Amyloid
        Random_Index = rand(1:length(keys_list))
        Random_Location = keys_list[Random_Index]
        State, _ = Locations_and_States_Dict[Random_Location]

        if State == 0
            Assigns_State_Monomer_Amyloid(Random_Location)
            Monomers_Made_Amyloid += 1
        end
    end
end

"""
    Assigns_State_Monomer_Native(Location::Tuple)

Assigns a native state (1) and a unique identifier to the specified lattice coordinate.

# Arguments
- `Location`: Tuple representing (x, y, z) coordinate.
"""

# Functions to assign states to specific coordinates
function Assigns_State_Monomer_Native(Location)
    global Locations_and_States_Dict
    Unique_Number = Randomly_Choosing_Unique_Number_Monomer()
    
    # Assign state as Native and keep the unique number
    Locations_and_States_Dict[Location] = (NativeState_Value, Unique_Number)
end

"""
    Assigns_State_Monomer_Amyloid(Location::Tuple)

Assigns an amyloid-prone state (2) and a unique identifier to the specified lattice coordinate.

# Arguments
- `Location`: Tuple representing (x, y, z) coordinate.
"""

function Assigns_State_Monomer_Amyloid(Location)
    global Locations_and_States_Dict
    Unique_Number = Randomly_Choosing_Unique_Number_Monomer()
    # Assign state as Native and keep the unique number
    Locations_and_States_Dict[Location] = (AmyloidProne_Value, Unique_Number)
end

"""
    Randomly_Choosing_Unique_Number_Monomer() -> Int

Returns a randomly selected, unused monomer unique identifier.
Removes the identifier from the pool to avoid duplication.

# Returns
- A unique integer identifier for a monomer.
"""

function Randomly_Choosing_Unique_Number_Monomer()
    if length(Monomer_Unique_Numbers) == 0
        error("No more unique numbers available!")
    end

    # Randomly select and remove a number from the available list
    idx = rand(1:length(Monomer_Unique_Numbers))
    UniqueCode = Monomer_Unique_Numbers[idx]
    deleteat!(Monomer_Unique_Numbers, idx)  # Remove the selected number to ensure uniqueness
    
    return UniqueCode
end

#----------Sphere Crowder----------------------
"""
    Differentiate_Sphere_Crowder_Radius()

Calls the appropriate sphere generation routine based on the obstacle radius parameter.
"""

function Differentiate_Sphere_Crowder_Radius()
    if Obstacle_Radius == 0
        Generate_Spherical_Crowders_Radius_0()
    else
        Generate_Spherical_Crowders()
    end

end


"""
    Generate_Spherical_Crowders_Radius_0()

Generates spherical crowders of radius 0 by placing non-overlapping spheres in the lattice.
The number of spheres is determined by the target occupied fraction.
"""

function Generate_Spherical_Crowders_Radius_0()
    println("We are in the funciton Generate_Spherical_Crowders_Radius_1")

    # Determine the target number of spheres to create
    target_spheres = Calculate_Target_Number_of_Spheres()
    generated_spheres = 0

    # Generate spheres until we reach the target count
    while generated_spheres < target_spheres
        # Attempt to generate a new sphere
        success = Making_Spheres_Radius_0()

        if success
            generated_spheres += 1
            println("Sphere #$generated_spheres created successfully.")
        else
            #println("Failed to create a sphere; retrying.")
        end
    end
    println("Finished generating spheres. Total spheres created: $generated_spheres.")

end

"""
    Making_Spheres_Radius_0() -> Bool

Attempts to place a single sphere crowder at a randomly selected valid lattice coordinate.

# Returns
- `true` if the sphere was successfully placed, `false` otherwise.
"""

function Making_Spheres_Radius_0()
    global Locations_and_States_Dict
    keys_list = collect(keys(Locations_and_States_Dict))

        Random_Index = rand(1:length(keys_list))
        Random_Location = keys_list[Random_Index]
        State, _ = Locations_and_States_Dict[Random_Location]

        if State == 0
           Sphere_Unique_Number = Randomly_Choose_Unique_Number_Sphere()
           Assigns_State_Monomer_Sphere(Random_Location,  Sphere_Unique_Number)
           return true
        else 
            return false
        end
end

"""
    Generate_Spherical_Crowders()

Places spherical crowders of arbitrary radius into the lattice using a composite sphere construction routine.
"""

# Function to create a spherical shape of crowders
function Generate_Spherical_Crowders()

    # Determine the target number of spheres to create
    target_spheres = Calculate_Target_Number_of_Spheres()
    generated_spheres = 0

    # Generate spheres until we reach the target count
    while generated_spheres < target_spheres
        # Attempt to generate a new sphere
        success = Calling_Sphere_Coordinate_Functions()

        if success
            generated_spheres += 1
            println("Sphere #$generated_spheres created successfully.")
        else
            #println("Failed to create a sphere; retrying.")
        end
    end
    println("Finished generating spheres. Total spheres created: $generated_spheres.")
    
    #Calculate_Target_Number_of_Spheres()
    #Calling_Sphere_Coordinate_Functions()
end

"""
    Randomly_Decide_Point() -> Tuple{Float64, Float64, Float64}

Randomly selects a valid, unused lattice coordinate to serve as the center of a new sphere.

# Returns
- A tuple containing (x, y, z) coordinates.
"""

function Randomly_Decide_Point()
    # Collect all coordinates (keys) from Locations_and_States_Dict
    keys_list = collect(keys(Locations_and_States_Dict))
    selected_coordinate = nothing
    
    # Continue until we find a coordinate with a state value of 0 and not used before
    while selected_coordinate === nothing
        # Select a random coordinate from keys_list
        random_index = rand(1:length(keys_list))
        candidate_coordinate = keys_list[random_index]
        
        # Check if the coordinate has already been used as a center
        if candidate_coordinate in used_centers
            continue  # Skip to the next iteration
        end
        
        # Check the state of the selected coordinate
        state, _ = Locations_and_States_Dict[candidate_coordinate]
        
        # If the state is 0, select it and add it to the used_centers set
        if state == 0
            selected_coordinate = candidate_coordinate
            push!(used_centers, selected_coordinate)  # Mark this coordinate as used
        end
    end
    
    # Extract the x, y, and z values from the selected coordinate
    centerX, centerY, centerZ = selected_coordinate
    
    # Return the selected coordinate as the center point
    return centerX, centerY, centerZ
end

"""
    Calling_Sphere_Coordinate_Functions() -> Bool

Attempts to create a spherical crowder centered at a selected point by computing surrounding coordinates.

# Returns
- `true` if all directional components succeed, otherwise `false`.
"""

function Calling_Sphere_Coordinate_Functions()
    success = false
    while !success
        # Choose a new center point
        center_x, center_y, center_z = Randomly_Decide_Point()
        
        # Try to create the sphere in each direction
        if calculate_sphere_coordinates_left_center_x(center_x, center_y, center_z) &&
           calculate_sphere_coordinates_right_center_x(center_x, center_y, center_z) &&
           calculate_sphere_coordinates_forward_center_y(center_x, center_y, center_z) &&
           calculate_sphere_coordinates_backward_center_y(center_x, center_y, center_z) &&
           calculate_sphere_coordinates_upward_center_z(center_x, center_y, center_z) &&
           calculate_sphere_coordinates_downward_center_z(center_x, center_y, center_z)
           
            # If all functions succeed without clearing the set, we have successfully created a sphere
            success = true
            println("Successfully created a sphere at center ($center_x, $center_y, $center_z).")
            determine_growth_direction()
            Change_State_of_Sphere_Coordinates()  # Finalize state changes
            Empty_Possible_Sphere_Coordinates()
            return true  # Return true indicating success
        else
            println("Failed to create a valid sphere. Retrying with a new center.")
            return false  # Return false indicating failure to create the sphere
        end
    end
end

"""
    calculate_sphere_coordinates_left_center_x(center_x, center_y, center_z) -> Bool

Attempts to identify and validate lattice coordinates to the left of `center_x` within a specified radius for sphere placement.

Returns `true` if all selected coordinates are valid for sphere assignment.
"""


function calculate_sphere_coordinates_left_center_x(center_x, center_y, center_z)
    global Possible_Sphere_Coordinates_Set
    Push_Coordinate_Sphere(center_x, center_y, center_z)

    # Filter coordinates to those left of center_x
    #left_coordinates = filter_coordinates_left_of_center(center_x)
    #println("These are the left_coordinates: $left_coordinates")

    #for coordinate in left_coordinates
    for coordinate in keys(Locations_and_States_Dict)
        state, _ = Locations_and_States_Dict[coordinate]
        coordinate_x, coordinate_y, coordinate_z = coordinate
        distance = X_Coordinate_Left(center_x, center_y, center_z, coordinate_x, coordinate_y, coordinate_z)

        if distance <= Obstacle_Radius
            if state == 0 
                Push_Coordinate_Sphere(coordinate_x, coordinate_y, coordinate_z)
            else
                Empty_Possible_Sphere_Coordinates()
                return false  # Failure due to an invalid state
            end
        end
    end
    return true  # Success if loop completes without clearing the set
end

"""
    calculate_sphere_coordinates_right_center_x(center_x, center_y, center_z) -> Bool

Attempts to identify and validate lattice coordinates to the right of `center_x` within a specified radius for sphere placement.

Returns `true` if all selected coordinates are valid for sphere assignment.
"""

function calculate_sphere_coordinates_right_center_x(center_x, center_y, center_z)
    global Possible_Sphere_Coordinates_Set

    #right_coordinates = filter_coordinates_right_of_center(center_x)
    #println("These are the right_coordinates: $right_coordinates")

    #for coordinate in right_coordinates
    for coordinate in keys(Locations_and_States_Dict)
        state, _ = Locations_and_States_Dict[coordinate]
        coordinate_x, coordinate_y, coordinate_z = coordinate
        distance = X_Coordinate_Right(center_x, center_y, center_z, coordinate_x, coordinate_y, coordinate_z)

        if distance <= Obstacle_Radius
            if state == 0 
                Push_Coordinate_Sphere(coordinate_x, coordinate_y, coordinate_z)
            else
                Empty_Possible_Sphere_Coordinates()
                return false
            end
        end
    end
    return true
end

"""
    calculate_sphere_coordinates_forward_center_y(center_x, center_y, center_z) -> Bool

Attempts to identify and validate lattice coordinates in the forward (positive Y) direction within a specified radius.

Returns `true` if valid sphere coordinates are found.
"""

function calculate_sphere_coordinates_forward_center_y(center_x, center_y, center_z)
    global Possible_Sphere_Coordinates_Set

    #forward_coordinates = filter_coordinates_forward_of_center(center_y)
    #println("These are the forward_coordinates: $forward_coordinates")

    #for coordinate in forward_coordinates
    for coordinate in keys(Locations_and_States_Dict)
        state, _ = Locations_and_States_Dict[coordinate]
        coordinate_x, coordinate_y, coordinate_z = coordinate
        distance = Y_Coordinate_Forward(center_x, center_y, center_z, coordinate_x, coordinate_y, coordinate_z)

        if distance <= Obstacle_Radius
            if state == 0 
                Push_Coordinate_Sphere(coordinate_x, coordinate_y, coordinate_z)
            else
                Empty_Possible_Sphere_Coordinates()
                return false
            end
        end
    end
    return true
end

"""
    calculate_sphere_coordinates_backward_center_y(center_x, center_y, center_z) -> Bool

Attempts to identify and validate lattice coordinates in the backward (negative Y) direction within a specified radius.

Returns `true` if valid sphere coordinates are found.
"""

function calculate_sphere_coordinates_backward_center_y(center_x, center_y, center_z)
    global Possible_Sphere_Coordinates_Set

    #backward_coordinates = filter_coordinates_backward_of_center(center_y)
    #println("These are the backward_coordinates: $backward_coordinates")

    #for coordinate in backward_coordinates
    for coordinate in keys(Locations_and_States_Dict)
        state, _ = Locations_and_States_Dict[coordinate]
        coordinate_x, coordinate_y, coordinate_z = coordinate
        distance = Y_Coordinate_Backward(center_x, center_y, center_z, coordinate_x, coordinate_y, coordinate_z)

        if distance <= Obstacle_Radius
            if state == 0 
                Push_Coordinate_Sphere(coordinate_x, coordinate_y, coordinate_z)
            else
                Empty_Possible_Sphere_Coordinates()
                return false
            end
        end
    end
    return true
end

"""
    calculate_sphere_coordinates_upward_center_z(center_x, center_y, center_z) -> Bool

Identifies valid sphere coordinates upward (positive Z) from a central point.

Returns `true` if all candidates are unoccupied.
"""

function calculate_sphere_coordinates_upward_center_z(center_x, center_y, center_z)
    global Possible_Sphere_Coordinates_Set

    #upward_coordinates = filter_coordinates_upward_of_center(center_z)
    #println("These are the upward_coordinates: $upward_coordinates")

    #for coordinate in upward_coordinates
    for coordinate in keys(Locations_and_States_Dict)
        state, _ = Locations_and_States_Dict[coordinate]
        coordinate_x, coordinate_y, coordinate_z = coordinate
        distance = Z_Coordiante_Upward(center_x, center_y, center_z, coordinate_x, coordinate_y, coordinate_z)
        if distance <= Obstacle_Radius
            if state == 0 
                Push_Coordinate_Sphere(coordinate_x, coordinate_y, coordinate_z)
            else
                Empty_Possible_Sphere_Coordinates()
                return false
            end
        end
    end
    return true
end

"""
    calculate_sphere_coordinates_downward_center_z(center_x, center_y, center_z) -> Bool

Identifies valid sphere coordinates downward (negative Z) from a central point.

Returns `true` if all selected lattice points are available.
"""

function calculate_sphere_coordinates_downward_center_z(center_x, center_y, center_z)
    global Possible_Sphere_Coordinates_Set

    #downward_coordinates = filter_coordinates_downward_of_center(center_z)
    #println("These are the downward_coordinates: $downward_coordinates")

    #for coordinate in downward_coordinates
    for coordinate in keys(Locations_and_States_Dict)
        state, _ = Locations_and_States_Dict[coordinate]
        coordinate_x, coordinate_y, coordinate_z = coordinate
        distance = Z_Coordinate_Downward(center_x, center_y, center_z, coordinate_x, coordinate_y, coordinate_z)

        if distance <= Obstacle_Radius
            if state == 0 
                Push_Coordinate_Sphere(coordinate_x, coordinate_y, coordinate_z)
            else
                Empty_Possible_Sphere_Coordinates()
                return false
            end
        end
    end
    return true
end


"""
    filter_coordinates_left_of_center(center_x) -> Vector{Tuple}

Returns all coordinates from `Locations_and_States_Dict` where the X-value is less than or equal to `center_x`.
"""

function filter_coordinates_left_of_center(center_x)
    # Return only coordinates with x-value less than center_x
    return filter(coord -> coord[1] <= center_x, keys(Locations_and_States_Dict))
end


"""
    filter_coordinates_right_of_center(center_x) -> Vector{Tuple}

Returns all coordinates where the X-value is greater than `center_x`.
"""

function filter_coordinates_right_of_center(center_x)
    return filter(coord -> coord[1] > center_x, keys(Locations_and_States_Dict))
end

"""
    filter_coordinates_forward_of_center(center_y) -> Vector{Tuple}

Returns all coordinates where the Y-value is greater than `center_y`.
"""

function filter_coordinates_forward_of_center(center_y)
    return filter(coord -> coord[2] > center_y, keys(Locations_and_States_Dict))
end

"""
    filter_coordinates_backward_of_center(center_y) -> Vector{Tuple}

Returns all coordinates where the Y-value is less than or equal to `center_y`.
"""

function filter_coordinates_backward_of_center(center_y)
    return filter(coord -> coord[2] <= center_y, keys(Locations_and_States_Dict))
end


"""
    filter_coordinates_upward_of_center(center_z) -> Vector{Tuple}

Returns all coordinates where the Z-value is greater than `center_z`.
"""

function filter_coordinates_upward_of_center(center_z)
    return filter(coord -> coord[3] > center_z, keys(Locations_and_States_Dict))
end

"""
    filter_coordinates_downward_of_center(center_z) -> Vector{Tuple}

Returns all coordinates where the Z-value is less than or equal to `center_z`.
"""

function filter_coordinates_downward_of_center(center_z)
    return filter(coord -> coord[3] <= center_z, keys(Locations_and_States_Dict))
end

"""
    filter_coordinates_downward_of_center(center_z) -> Vector{Tuple}

Returns all coordinates where the Z-value is less than or equal to `center_z`.
"""

function distance_from_center(center_x, center_y, center_z, coordinate_x, coordinate_y, coordinate_z)
    distance = sqrt((coordinate_x - center_x)^2 + (coordinate_y - center_y)^2 + (coordinate_z - center_z)^2)
    return distance
end

"""
    distance_from_center_X_Coordinate_Exception(...) -> Float64

Calculates Euclidean distance when X coordinate wraps around due to periodic boundary conditions.
"""

function distance_from_center_X_Coordinate_Exception(center_x, center_y, center_z, distance_x, coordinate_y, coordinate_z)
    distance = sqrt((distance_x)^2 + (coordinate_y - center_y)^2 + (coordinate_z - center_z)^2)
    return distance
end

"""
    distance_from_center_Y_Coordinate_Exception(...) -> Float64

Calculates distance when Y coordinate wraps due to boundary crossing.
"""

function distance_from_center_Y_Coordinate_Exception(center_x, center_y, center_z, coordinate_x, distance_y, coordinate_z)
    distance = sqrt((coordinate_x - center_x)^2 + (distance_y)^2 + (coordinate_z - center_z)^2)
    return distance
end

"""
    distance_from_center_Z_Coordinate_Exception(...) -> Float64

Calculates distance when Z coordinate wraps due to boundary crossing.
"""

function distance_from_center_Z_Coordinate_Exception(center_x, center_y, center_z, coordinate_x, coordinate_y, distance_z)
    distance = sqrt((coordinate_x - center_x)^2 + (coordinate_y - center_y)^2 + (distance_z)^2)
    return distance
end

"""
    X_Coordinate_Right(...) -> Float64

Calculates distance from center in the +X direction, accounting for wrapping if needed.
"""

function X_Coordinate_Right(center_x, center_y, center_z, coordinate_x, coordinate_y, coordinate_z)
    if coordinate_x < center_x #MAKE X COORDINATE OUT OF BOUNDS
        distance = (Lattice_Size - center_x) + coordinate_x + 1
        distance_from_center_X_Coordinate_Exception(center_x, center_y, center_z, distance, coordinate_y, coordinate_z)
    else
        distance_from_center(center_x, center_y, center_z, coordinate_x, coordinate_y, coordinate_z)
    end
end

"""
    X_Coordinate_Left(...) -> Float64

Calculates distance from center in the -X direction, including wraparound if applicable.
"""

function X_Coordinate_Left(center_x, center_y, center_z, coordinate_x, coordinate_y, coordinate_z)
    if coordinate_x > center_x #MAKE X COORDINATE OUT OF BOUNDS
        distance = center_x + 1 + (Lattice_Size - coordinate_x)
        #println("For coordinate_x: $coordinate_x the new distance is: $distance")
        distance_from_center_X_Coordinate_Exception(center_x, center_y, center_z, distance, coordinate_y, coordinate_z)
    else
        distance_from_center(center_x, center_y, center_z, coordinate_x, coordinate_y, coordinate_z)
    end 
end

"""
    Y_Coordinate_Forward(...) -> Float64

Calculates distance from center in the +Y direction, with wraparound.
"""

function Y_Coordinate_Forward(center_x, center_y, center_z, coordinate_x, coordinate_y, coordinate_z)
    if coordinate_y < center_y #MAKE Y COORDINATE OUT OF BOUNDS
        distance = (Lattice_Size - center_y) + coordinate_y + 1
        distance_from_center_Y_Coordinate_Exception(center_x, center_y, center_z, coordinate_x, distance, coordinate_z)
    else
        distance_from_center(center_x, center_y, center_z, coordinate_x, coordinate_y, coordinate_z)
    end

end

"""
    Y_Coordinate_Backward(...) -> Float64

Calculates distance from center in the -Y direction, accounting for boundary wrap.
"""

function Y_Coordinate_Backward(center_x, center_y, center_z, coordinate_x, coordinate_y, coordinate_z)
    if coordinate_y > center_y #MAKE Y COORDINATE OUT OF BOUNDS
        distance = center_y + 1 + (Lattice_Size - coordinate_y)
        distance_from_center_Y_Coordinate_Exception(center_x, center_y, center_z, coordinate_x, distance, coordinate_z)
    else
        distance_from_center(center_x, center_y, center_z, coordinate_x, coordinate_y, coordinate_z)
    end
    
end

"""
    Z_Coordiante_Upward(...) -> Float64

Calculates distance from center in the +Z direction, including periodic wrapping.
"""

function Z_Coordiante_Upward(center_x, center_y, center_z, coordinate_x, coordinate_y, coordinate_z)
    if coordinate_z < center_z #MAKE Z COORDINATE OUT OF BOUNDS
        distance = (Lattice_Size - center_z) + coordinate_z + 1
        distance_from_center_Z_Coordinate_Exception(center_x, center_y, center_z, coordinate_x, coordinate_y, distance)
    else
        distance_from_center(center_x, center_y, center_z, coordinate_x, coordinate_y, coordinate_z)
    end
end

"""
    Z_Coordinate_Downward(...) -> Float64

Calculates distance from center in the -Z direction, using exception logic for wraparound.
"""

function Z_Coordinate_Downward(center_x, center_y, center_z, coordinate_x, coordinate_y, coordinate_z)
    if coordinate_z > center_z #MAKE Z COORDINATE OUT OF BOUNDS
        distance =  center_z + 1 + (Lattice_Size - coordinate_z)
        distance_from_center_Z_Coordinate_Exception(center_x, center_y, center_z, coordinate_x, coordinate_y, distance)
    else
        distance_from_center(center_x, center_y, center_z, coordinate_x, coordinate_y, coordinate_z)
    end

end

"""
    Change_State_of_Sphere_Coordinates()

Assigns the Sphere state and a unique identifier to all coordinates in `Possible_Sphere_Coordinates_Set`.
"""

function Change_State_of_Sphere_Coordinates()
    # Iterate through each coordinate in Possible_Sphere_Coordinates_Set
    #println("This is Possible_Sphere_Coordinates_Set: $Possible_Sphere_Coordinates_Set")
    println("Total number of coordinates in Possible_Sphere_Coordinates_Set: ", length(Possible_Sphere_Coordinates_Set))
    
    # Get a unique number for this sphere
    Sphere_Unique_Number = Randomly_Choose_Unique_Number_Sphere()
    println("Sphere_Unique_Number: $Sphere_Unique_Number")
    
    # Use `keys(Possible_Sphere_Coordinates_Set)` to get only the coordinates
    for coordinate in keys(Possible_Sphere_Coordinates_Set)
        # Assign state and unique number to each coordinate
        Assigns_State_Monomer_Sphere(coordinate, Sphere_Unique_Number)
    end
end

"""
    Randomly_Choose_Unique_Number_Sphere() -> Int

Selects a unique identifier for a new sphere and removes it from the available pool.

# Returns
- A unique integer ID for a sphere.
"""

function Randomly_Choose_Unique_Number_Sphere()
    if length(Sphere_Unique_Numbers) == 0
        error("No more unique numbers available!")
    end

    # Randomly select and remove a number from the available list
    idx = rand(1:length(Sphere_Unique_Numbers))
    UniqueCode = Sphere_Unique_Numbers[idx]
    deleteat!(Sphere_Unique_Numbers, idx)  # Remove the selected number to ensure uniqueness
    
    return UniqueCode
end

"""
    Assigns_State_Monomer_Sphere(Location::Tuple, Sphere_Unique_Number::Int)

Assigns sphere state (5) and a unique ID to the given coordinate.
"""

function Assigns_State_Monomer_Sphere(Location, Sphere_Unique_Number)
    global Locations_and_States_Dict
    _, Unique_Number = Locations_and_States_Dict[Location]
    Locations_and_States_Dict[Location] = (SphereState_Value, Sphere_Unique_Number)
end

"""
    Empty_Possible_Sphere_Coordinates()

Clears the global dictionary used for temporary storage of sphere coordinates.
"""

function Empty_Possible_Sphere_Coordinates()
    empty!(Possible_Sphere_Coordinates_Set)  
end

"""
    Push_Coordinate_Sphere(X::Float64, Y::Float64, Z::Float64)

Adds the specified coordinate to `Possible_Sphere_Coordinates_Set`.
"""

function Push_Coordinate_Sphere(X, Y, Z)
    Possible_Sphere_Coordinates_Set[(X, Y, Z)] = nothing
end

"""
    determine_growth_direction()

Analyzes the spread of the current sphere and prints the primary growth axis (X, Y, or Z).
"""


function determine_growth_direction()
    x_values = [coord[1] for coord in keys(Possible_Sphere_Coordinates_Set)]
    y_values = [coord[2] for coord in keys(Possible_Sphere_Coordinates_Set)]
    z_values = [coord[3] for coord in keys(Possible_Sphere_Coordinates_Set)]

    x_range = maximum(x_values) - minimum(x_values)
    y_range = maximum(y_values) - minimum(y_values)
    z_range = maximum(z_values) - minimum(z_values)

    if x_range > y_range && x_range > z_range
        println("The sphere is growing primarily in the X direction.")
    elseif y_range > x_range && y_range > z_range
        println("The sphere is growing primarily in the Y direction.")
    elseif z_range > x_range && z_range > y_range
        println("The sphere is growing primarily in the Z direction.")
    else
        println("The sphere is growing equally in multiple directions.")
    end

    println("X range: $x_range, Y range: $y_range, Z range: $z_range")
end

"""
    Calculate_Target_Number_of_Spheres() -> Int

Estimates how many spheres should be placed to reach the target crowding concentration.

# Returns
- Integer number of spheres to generate.
"""

function Calculate_Target_Number_of_Spheres()
  # Calculate the total number of lattice locations
  total_locations = length(keys(Locations_and_States_Dict))
    
  # Calculate the target number of occupied spaces based on crowder concentration
  occupied_spaces = Crowder_Concentration_Spheres * total_locations
  target_occupied_spaces = round(Int, occupied_spaces)

  # Approximate volume of one sphere in lattice coordinates (from observation)
  #Sphere_Volume = 1 # For a sphere of radius 1 in FCC
  
  # Estimate number of spheres needed
  target_number_of_spheres = round(Int, target_occupied_spaces / Sphere_Volume)

  println("Total lattice locations: $total_locations")
  println("Target occupied spaces based on concentration: $target_occupied_spaces")
  println("Estimated number of spheres needed: $target_number_of_spheres")

  return target_number_of_spheres
end

"""
    Count_Number_Coordinates_Spheres() -> Int

Counts how many lattice coordinates are currently marked as part of a sphere.

# Returns
- Integer count of sphere-assigned coordinates.
"""

function Count_Number_Coordinates_Spheres()
    sphere_count = 0

    # Loop through the dictionary to count monomers in state 5
    for (_, (state, _)) in Locations_and_States_Dict
        if state == 5
            sphere_count += 1
        end
    end

    return  sphere_count
end

#-------------------------------------------------

Generate_Coordinates(Lattice_Size)
