using Random
using Plots
using CSV
using DataFrames

"""
Module: Agents.jl

Handles the initialization and management of agents (monomers and crowders) within a 3D FCC lattice environment.
Defines the structure of the simulation space, state variables, and random assignment of monomers and crowders.

Implements FAIR principles:

- **Findable**: Clear structure for agent generation and unique identifiers.
- **Accessible**: Simulation parameters and initial conditions are externally loadable via CSV files.
- **Interoperable**: Built using Julia-native types and standard libraries.
- **Reusable**: Modular functions with transparent global variables, enabling reproducibility across simulations.

Key Features:
- Supports initialization of Native (N), Amyloid-prone (A) monomers, Oligomers (O), and Fibrils (F).
- Crowders are modeled as static spheres generated at random locations.
- Periodic boundary conditions incorporated for crowder creation.
- All agents are tracked with unique numerical identifiers.

Authors: Santiago Schnell; Conner Sandefur; Isabella Gimon
Dependencies: 
- Random
- Plots
- CSV
- DataFrames

License: http://www.apache.org/licenses/LICENSE-2.0
"""


###############
# Constants
###############

# State values for agents
const NativeState_Value         = 1  # Native monomer
const AggregateProne_Value        = 2  # Amyloid-prone monomer
const OligomerState_Value       = 3  # Oligomer
const FibrilState_Value         = 4  # Fibril
const SphereState_Value         = 5  # Crowder/Sphere

# Unique number ranges
const MONOMER_ID_RANGE     = 200_001:499_999
const SPHERE_ID_RANGE      = 500_000:2_000_000

###############
# Global Simulation State
###############

# Dictionary mapping 3D coordinates (Tuple{Float64, Float64, Float64}) to a Tuple containing (State, UniqueID)
# Stores the current state of the 3D lattice.
global Locations_and_States_Dict = Dict{Tuple{Float64, Float64, Float64}, Tuple{Int, Int}}()

# Temporary storage for potential sphere crowder coordinates.
global Possible_Sphere_Coordinates_Set = Dict{Tuple{Float64, Float64, Float64}, Nothing}()

# Snapshot of initial monomer locations (before simulation dynamics).
global Initial_Locations_and_States_Dict = Dict{Tuple{Float64, Float64, Float64}, Tuple{Int, Int}}()

# Used centers for crowder spheres to avoid overlapping.
global used_centers = Set{Tuple{Float64, Float64, Float64}}()

# Pools of unique IDs
global Sphere_Unique_Numbers = collect(SPHERE_ID_RANGE)
global Monomer_Unique_Numbers = collect(MONOMER_ID_RANGE)

"""
    load_csv_parameters(file_path::String) -> Dict

Loads simulation parameters from a CSV file into a dictionary.

The input CSV should have two columns: one for parameter names and one for their values.  
Automatically detects Float64, Bool (`TRUE`/`FALSE`), or String values.

# Arguments
- `file_path::String`: Path to the CSV file.

# Returns
- `Dict{String, Any}`: Dictionary mapping parameter names to their corresponding parsed values.

# Behavior
- Standardizes column names to remove hidden characters.
- Ensures selected fields (e.g., "Lattice_Size") are stored as integers.
"""
function load_csv_parameters(file_path::String)
    df = CSV.read(file_path, DataFrame)

    rename!(df, strip.(names(df)))
    println("Updated column names: ", names(df))

    param_col = names(df)[1]
    value_col = names(df)[2]

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

    for key in ["Lattice_Size", "Max_NumberMonomers_Native", "Max_NumberMonomers_AggregateProne", "Obstacle_Radius", "MAX_NumberMovements"]
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

# Assign simulation parameters
Lattice_Size = Parameters["Lattice_Size"]                                 # Size of the cubic lattice (Lattice_Size x Lattice_Size x Lattice_Size)
Max_NumberMonomers_Native = Parameters["Max_NumberMonomers_Native"]       # Maximum number of native monomers
Max_NumberMonomers_AggregateProne = Parameters["Max_NumberMonomers_AggregateProne"]     # Maximum number of amyloid-prone monomers
Obstacle_Radius = Parameters["Obstacle_Radius"]                           # Radius of spherical crowders
Crowder_Concentration_Spheres = Parameters["Crowder_Concentration_Spheres"] # Crowder concentration
Obstacle = Parameters["Spheres?"]                                         # Boolean: Enable/Disable crowders
#Sphere_Volume = Parameters["Sphere_Volume"]                               # Volume of a sphere

# Debugging: Print loaded parameters
println("Loaded Parameters: ", Parameters)

"""
    Generate_Coordinates(Lattice_Size::Int)

Generates a complete 3D face-centered cubic lattice.

Populates the global `Locations_and_States_Dict` with all corner and face-centered points
at each unit cell position.  
Also randomly assigns monomer states (Native and Amyloid-prone) and optionally generates spherical crowders.

# Arguments
- `Lattice_Size::Int`: The number of unit cells per axis (X, Y, Z).

# Global variables modified
- `Locations_and_States_Dict`
- `Initial_Locations_and_States_Dict`
- `Possible_Sphere_Coordinates_Set` (if Obstacle enabled)

# Calls
- `Add_Position`
- `Differentiate_Sphere_Crowder_Radius`
- `Randomly_Assigns_Location_Monomers_Native`
- `Randomly_Assigns_Location_Monomers_Amyloid`
- `Copy_Original_Location`
"""
function Generate_Coordinates(Lattice_Size)
    global Locations_and_States_Dict

    for X in 0:(Lattice_Size - 1)
        for Y in 0:(Lattice_Size - 1)
            for Z in 0:(Lattice_Size - 1)
                State = 0
                Unique_Number = 0
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

Copies all Native and Amyloid-prone monomer coordinates from `Locations_and_States_Dict`
into `Initial_Locations_and_States_Dict`.

This function preserves a copy of the initial state of the simulation,
before any dynamics (movement, aggregation) occur.

# Global variables modified
- `Initial_Locations_and_States_Dict`

# Notes
- Only state values 1 (Native) and 2 (Amyloid-prone) are copied.
"""
function Copy_Original_Location()
    global Initial_Locations_and_States_Dict

    for (location, (state, unique_number)) in Locations_and_States_Dict
        if state == NativeState_Value || state == AggregateProne_Value
            Initial_Locations_and_States_Dict[location] = (state, unique_number)
        end
    end
end



"""
    Add_Position(X, Y, Z, State, Unique_Number)

Adds all face-centered and corner-centered coordinates for a unit cell at (X, Y, Z)
into the global `Locations_and_States_Dict`.

# Arguments
- `X, Y, Z`: Integer coordinates defining the origin of the unit cell.
- `State`: Initial state assigned to all created positions (typically 0).
- `Unique_Number`: Initial unique number assigned (typically 0).

# Behavior
- Populates 14 points per unit cell: 8 corners + 6 face centers.
"""
function Add_Position(X, Y, Z, State, Unique_Number)
    global Locations_and_States_Dict

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

######################
# Corner Coordinates
######################

"""
    First_Corner_Position(X, Y, Z) -> Tuple

Returns the coordinate for the first corner of the FCC unit cell: the base origin (X, Y, Z).
"""
function First_Corner_Position(X, Y, Z)
    return X, Y, Z
end

"""
    Second_Corner_Position(X, Y, Z) -> Tuple

Returns the coordinate of the second corner, displaced by +1 along the X-axis.
"""
function Second_Corner_Position(X, Y, Z)
    return X + 1, Y, Z
end

"""
    Third_Corner_Position(X, Y, Z) -> Tuple

Returns the coordinate of the third corner, displaced by +1 along the Y-axis.
"""
function Third_Corner_Position(X, Y, Z)
    return X, Y + 1, Z
end

"""
    Fourth_Corner_Position(X, Y, Z) -> Tuple

Returns the coordinate of the fourth corner, displaced by +1 along the Z-axis.
"""
function Fourth_Corner_Position(X, Y, Z)
    return X, Y, Z + 1
end

"""
    Fifth_Corner_Position(X, Y, Z) -> Tuple

Returns the coordinate at (X+1, Y+1, Z), a corner displaced along X and Y axes.
"""
function Fifth_Corner_Position(X, Y, Z)
    return X + 1, Y + 1, Z
end

"""
    Sixth_Corner_Position(X, Y, Z) -> Tuple

Returns the coordinate at (X+1, Y, Z+1), a corner displaced along X and Z axes.
"""
function Sixth_Corner_Position(X, Y, Z)
    return X + 1, Y, Z + 1
end

"""
    Seventh_Corner_Position(X, Y, Z) -> Tuple

Returns the coordinate at (X, Y+1, Z+1), a corner displaced along Y and Z axes.
"""
function Seventh_Corner_Position(X, Y, Z)
    return X, Y + 1, Z + 1
end

"""
    Eighth_Corner_Position(X, Y, Z) -> Tuple

Returns the coordinate at (X+1, Y+1, Z+1), displaced along all three axes.
"""
function Eighth_Corner_Position(X, Y, Z)
    return X + 1, Y + 1, Z + 1
end

######################
# Face Coordinates
######################

"""
    First_Face_Position(X, Y, Z) -> Tuple

Returns the coordinate at the center between (X, Y) plane: mid-point face center.
"""
function First_Face_Position(X, Y, Z)
    return X + 0.5, Y + 0.5, Z
end

"""
    Second_Face_Position(X, Y, Z) -> Tuple

Returns the coordinate at the center between (Y, Z) plane: mid-point face center.
"""
function Second_Face_Position(X, Y, Z)
    return X, Y + 0.5, Z + 0.5
end

"""
    Third_Face_Position(X, Y, Z) -> Tuple

Returns the coordinate at the center between (X, Z) plane: mid-point face center.
"""
function Third_Face_Position(X, Y, Z)
    return X + 0.5, Y, Z + 0.5
end

"""
    Fourth_Face_Position(X, Y, Z) -> Tuple

Returns the coordinate centered on the top face (+Z direction).
"""
function Fourth_Face_Position(X, Y, Z)
    return X + 0.5, Y + 0.5, Z + 1
end

"""
    Fifth_Face_Position(X, Y, Z) -> Tuple

Returns the coordinate centered on the front face (+Y direction).
"""
function Fifth_Face_Position(X, Y, Z)
    return X + 0.5, Y + 1, Z + 0.5
end

"""
    Sixth_Face_Position(X, Y, Z) -> Tuple

Returns the coordinate centered on the right face (+X direction).
"""
function Sixth_Face_Position(X, Y, Z)
    return X + 1, Y + 0.5, Z + 0.5
end

"""
    Randomly_Assigns_Location_Monomers_Native()

Randomly assigns the Native monomer state (1) to unoccupied lattice sites.

The number of monomers assigned is controlled by `Max_NumberMonomers_Native`.

# Global variables modified
- `Locations_and_States_Dict`

# Behavior
- Picks random empty coordinates until desired number is reached.
"""
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

Randomly assigns the Amyloid-prone monomer state (2) to unoccupied lattice sites.

The number of monomers assigned is controlled by `Max_NumberMonomers_AggregateProne`.

# Global variables modified
- `Locations_and_States_Dict`
"""
function Randomly_Assigns_Location_Monomers_Amyloid()
    global Locations_and_States_Dict
    println("We are in Randomly_Assigns_Location_Monomers_Amyloid")
    Monomers_Made_Amyloid = 0
    keys_list = collect(keys(Locations_and_States_Dict))

    while Monomers_Made_Amyloid < Max_NumberMonomers_AggregateProne
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

Assigns the Native monomer state (1) and a random unique ID to a given coordinate.

# Arguments
- `Location::Tuple`: (X, Y, Z) coordinate.

# Global variables modified
- `Locations_and_States_Dict`
"""
function Assigns_State_Monomer_Native(Location)
    global Locations_and_States_Dict
    Unique_Number = Randomly_Choosing_Unique_Number_Monomer()
    Locations_and_States_Dict[Location] = (NativeState_Value, Unique_Number)
end

"""
    Assigns_State_Monomer_Amyloid(Location::Tuple)

Assigns the Amyloid-prone monomer state (2) and a random unique ID to a given coordinate.

# Arguments
- `Location::Tuple`: (X, Y, Z) coordinate.
"""
function Assigns_State_Monomer_Amyloid(Location)
    global Locations_and_States_Dict
    Unique_Number = Randomly_Choosing_Unique_Number_Monomer()
    Locations_and_States_Dict[Location] = (AggregateProne_Value, Unique_Number)
end

"""
    Randomly_Choosing_Unique_Number_Monomer() -> Int

Selects and removes a unique ID from the monomer pool.

# Returns
- An integer unique identifier.

# Behavior
- Ensures IDs are not reused by deleting used ones.
"""
function Randomly_Choosing_Unique_Number_Monomer()
    if isempty(Monomer_Unique_Numbers)
        error("No more unique monomer IDs available!")
    end

    idx = rand(1:length(Monomer_Unique_Numbers))
    UniqueCode = Monomer_Unique_Numbers[idx]
    deleteat!(Monomer_Unique_Numbers, idx)
    
    return UniqueCode
end

# --- Sphere Crowder Functions (Summary) ---

"""
    Differentiate_Sphere_Crowder_Radius()

Chooses the sphere generation method based on radius parameter.

# Behavior
- Calls `Generate_Spherical_Crowders_Radius_0()` if radius = 0.
- Otherwise, calls `Generate_Spherical_Crowders()`.
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

Generates spherical crowders of radius zero (single occupied lattice points).
"""
function Generate_Spherical_Crowders_Radius_0()
    println("We are in the function Generate_Spherical_Crowders_Radius_0")

    target_spheres = Calculate_Target_Number_of_Spheres()
    generated_spheres = 0

    while generated_spheres < target_spheres
        success = Making_Spheres_Radius_0()
        if success
            generated_spheres += 1
            println("Sphere #$generated_spheres created successfully.")
        end
    end

    println("Finished generating spheres. Total spheres created: $generated_spheres.")
end

"""
    Making_Spheres_Radius_0() -> Bool

Attempts to assign a sphere crowder at a random empty coordinate.

# Returns
- `true` if successful, `false` if the site was already occupied.
"""
function Making_Spheres_Radius_0()
    global Locations_and_States_Dict
    keys_list = collect(keys(Locations_and_States_Dict))

    Random_Index = rand(1:length(keys_list))
    Random_Location = keys_list[Random_Index]
    State, _ = Locations_and_States_Dict[Random_Location]

    if State == 0
        Sphere_Unique_Number = Randomly_Choose_Unique_Number_Sphere()
        Assigns_State_Monomer_Sphere(Random_Location, Sphere_Unique_Number)
        return true
    else 
        return false
    end
end



"""
    Generate_Spherical_Crowders()

Places spherical crowders of arbitrary radius into the lattice using a composite sphere construction routine.

# Behavior
- Calls `Calling_Sphere_Coordinate_Functions()` iteratively until the desired number of spheres is generated.
- Each sphere occupies unoccupied lattice sites within a defined radius.
"""
function Generate_Spherical_Crowders()
    target_spheres = Calculate_Target_Number_of_Spheres()
    generated_spheres = 0

    while generated_spheres < target_spheres
        success = Calling_Sphere_Coordinate_Functions()

        if success
            generated_spheres += 1
            println("Sphere #$generated_spheres created successfully.")
        end
    end

    println("Finished generating spheres. Total spheres created: $generated_spheres.")
end

"""
    Randomly_Decide_Point() -> Tuple{Float64, Float64, Float64}

Randomly selects a valid, unused lattice coordinate to serve as the center of a new sphere.

# Returns
- A tuple `(x, y, z)` representing the selected center coordinate.
"""
function Randomly_Decide_Point()
    keys_list = collect(keys(Locations_and_States_Dict))
    selected_coordinate = nothing
    
    while selected_coordinate === nothing
        candidate_coordinate = keys_list[rand(1:length(keys_list))]
        
        if candidate_coordinate ∉ used_centers
            state, _ = Locations_and_States_Dict[candidate_coordinate]
            if state == 0
                selected_coordinate = candidate_coordinate
                push!(used_centers, selected_coordinate)
            end
        end
    end

    return selected_coordinate
end

"""
    Calling_Sphere_Coordinate_Functions() -> Bool

Attempts to construct a spherical crowder by validating lattice points 
in all six cardinal directions from a random center.

# Returns
- `true` if the sphere is successfully created.
- `false` if placement fails and must retry.
"""
function Calling_Sphere_Coordinate_Functions()
    center_x, center_y, center_z = Randomly_Decide_Point()

    if calculate_sphere_coordinates_left_center_x(center_x, center_y, center_z) &&
       calculate_sphere_coordinates_right_center_x(center_x, center_y, center_z) &&
       calculate_sphere_coordinates_forward_center_y(center_x, center_y, center_z) &&
       calculate_sphere_coordinates_backward_center_y(center_x, center_y, center_z) &&
       calculate_sphere_coordinates_upward_center_z(center_x, center_y, center_z) &&
       calculate_sphere_coordinates_downward_center_z(center_x, center_y, center_z)

        println("✅ Successfully created a sphere at center ($center_x, $center_y, $center_z).")
        determine_growth_direction()
        Change_State_of_Sphere_Coordinates()
        Empty_Possible_Sphere_Coordinates()
        return true
    else
        #println("⚠️ Failed to create a valid sphere. Retrying...")
        return false
    end
end

"""
    calculate_sphere_coordinates_left_center_x(center_x, center_y, center_z) -> Bool

Attempts to identify and validate lattice coordinates to the left of `center_x`
within the given radius for sphere construction.

Returns `true` if all selected coordinates are valid (unoccupied), otherwise `false`.
"""
function calculate_sphere_coordinates_left_center_x(center_x, center_y, center_z)
    global Possible_Sphere_Coordinates_Set
    Push_Coordinate_Sphere(center_x, center_y, center_z)

    for coordinate in keys(Locations_and_States_Dict)
        state, _ = Locations_and_States_Dict[coordinate]
        coordinate_x, coordinate_y, coordinate_z = coordinate
        distance = X_Coordinate_Left(center_x, center_y, center_z, coordinate_x, coordinate_y, coordinate_z)

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
    calculate_sphere_coordinates_right_center_x(center_x, center_y, center_z) -> Bool

Attempts to identify and validate lattice coordinates to the right of `center_x`
within the given radius for sphere construction.

Returns `true` if all selected coordinates are valid (unoccupied), otherwise `false`.
"""
function calculate_sphere_coordinates_right_center_x(center_x, center_y, center_z)
    global Possible_Sphere_Coordinates_Set

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

Attempts to identify and validate lattice coordinates in the forward (positive Y)
direction within the given radius for sphere construction.

Returns `true` if all selected coordinates are valid (unoccupied), otherwise `false`.
"""
function calculate_sphere_coordinates_forward_center_y(center_x, center_y, center_z)
    global Possible_Sphere_Coordinates_Set

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

Attempts to identify and validate lattice coordinates in the backward (negative Y)
direction within the given radius for sphere construction.

Returns `true` if all selected coordinates are valid (unoccupied), otherwise `false`.
"""
function calculate_sphere_coordinates_backward_center_y(center_x, center_y, center_z)
    global Possible_Sphere_Coordinates_Set

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

Attempts to identify and validate lattice coordinates upward (positive Z)
within the given radius for sphere construction.

Returns `true` if all selected coordinates are valid (unoccupied), otherwise `false`.
"""
function calculate_sphere_coordinates_upward_center_z(center_x, center_y, center_z)
    global Possible_Sphere_Coordinates_Set

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

Attempts to identify and validate lattice coordinates downward (negative Z)
within the given radius for sphere construction.

Returns `true` if all selected coordinates are valid (unoccupied), otherwise `false`.
"""
function calculate_sphere_coordinates_downward_center_z(center_x, center_y, center_z)
    global Possible_Sphere_Coordinates_Set

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
    distance_from_center(center_x, center_y, center_z, coordinate_x, coordinate_y, coordinate_z) -> Float64

Calculates the standard Euclidean distance between a center point and a coordinate in the lattice,
without applying any boundary conditions.

# Returns
- Distance as a Float64.
"""
function distance_from_center(center_x, center_y, center_z, coordinate_x, coordinate_y, coordinate_z)
    return sqrt((coordinate_x - center_x)^2 + (coordinate_y - center_y)^2 + (coordinate_z - center_z)^2)
end

"""
    distance_from_center_X_Coordinate_Exception(center_x, center_y, center_z, distance_x, coordinate_y, coordinate_z) -> Float64

Calculates Euclidean distance when an X-coordinate wraps around the periodic boundary.

# Arguments
- `distance_x`: Precomputed wrapped distance in the X direction.

# Returns
- Distance as a Float64.
"""
function distance_from_center_X_Coordinate_Exception(center_x, center_y, center_z, distance_x, coordinate_y, coordinate_z)
    return sqrt(distance_x^2 + (coordinate_y - center_y)^2 + (coordinate_z - center_z)^2)
end

"""
    distance_from_center_Y_Coordinate_Exception(center_x, center_y, center_z, coordinate_x, distance_y, coordinate_z) -> Float64

Calculates Euclidean distance when a Y-coordinate wraps around the periodic boundary.

# Arguments
- `distance_y`: Precomputed wrapped distance in the Y direction.

# Returns
- Distance as a Float64.
"""
function distance_from_center_Y_Coordinate_Exception(center_x, center_y, center_z, coordinate_x, distance_y, coordinate_z)
    return sqrt((coordinate_x - center_x)^2 + distance_y^2 + (coordinate_z - center_z)^2)
end
"""
    distance_from_center_Z_Coordinate_Exception(center_x, center_y, center_z, coordinate_x, coordinate_y, distance_z) -> Float64

Calculates Euclidean distance when a Z-coordinate wraps around the periodic boundary.

# Arguments
- `distance_z`: Precomputed wrapped distance in the Z direction.

# Returns
- Distance as a Float64.
"""
function distance_from_center_Z_Coordinate_Exception(center_x, center_y, center_z, coordinate_x, coordinate_y, distance_z)
    return sqrt((coordinate_x - center_x)^2 + (coordinate_y - center_y)^2 + distance_z^2)
end

"""
    X_Coordinate_Right(center_x, center_y, center_z, coordinate_x, coordinate_y, coordinate_z) -> Float64

Calculates the distance from the center to a coordinate in the positive X direction.
Handles periodic boundary crossing if the coordinate wraps around.

# Returns
- Distance as a Float64.
"""
function X_Coordinate_Right(center_x, center_y, center_z, coordinate_x, coordinate_y, coordinate_z)
    if coordinate_x < center_x  # Crossed periodic boundary
        distance_x = (Lattice_Size - center_x) + coordinate_x + 1
        return distance_from_center_X_Coordinate_Exception(center_x, center_y, center_z, distance_x, coordinate_y, coordinate_z)
    else
        return distance_from_center(center_x, center_y, center_z, coordinate_x, coordinate_y, coordinate_z)
    end
end

"""
    X_Coordinate_Left(center_x, center_y, center_z, coordinate_x, coordinate_y, coordinate_z) -> Float64

Calculates the distance from the center to a coordinate in the negative X direction.
Handles periodic boundary crossing if necessary.

# Returns
- Distance as a Float64.
"""
function X_Coordinate_Left(center_x, center_y, center_z, coordinate_x, coordinate_y, coordinate_z)
    if coordinate_x > center_x  # Crossed periodic boundary
        distance_x = center_x + 1 + (Lattice_Size - coordinate_x)
        return distance_from_center_X_Coordinate_Exception(center_x, center_y, center_z, distance_x, coordinate_y, coordinate_z)
    else
        return distance_from_center(center_x, center_y, center_z, coordinate_x, coordinate_y, coordinate_z)
    end
end

"""
    Y_Coordinate_Forward(center_x, center_y, center_z, coordinate_x, coordinate_y, coordinate_z) -> Float64

Calculates the distance from the center to a coordinate in the positive Y direction.
Handles wrapping if the coordinate crosses the periodic boundary.

# Returns
- Distance as a Float64.
"""
function Y_Coordinate_Forward(center_x, center_y, center_z, coordinate_x, coordinate_y, coordinate_z)
    if coordinate_y < center_y  # Crossed periodic boundary
        distance_y = (Lattice_Size - center_y) + coordinate_y + 1
        return distance_from_center_Y_Coordinate_Exception(center_x, center_y, center_z, coordinate_x, distance_y, coordinate_z)
    else
        return distance_from_center(center_x, center_y, center_z, coordinate_x, coordinate_y, coordinate_z)
    end
end


"""
    Y_Coordinate_Backward(center_x, center_y, center_z, coordinate_x, coordinate_y, coordinate_z) -> Float64

Calculates the distance from the center to a coordinate in the negative Y direction.
Handles wrapping if the coordinate crosses the boundary.

# Returns
- Distance as a Float64.
"""
function Y_Coordinate_Backward(center_x, center_y, center_z, coordinate_x, coordinate_y, coordinate_z)
    if coordinate_y > center_y  # Crossed periodic boundary
        distance_y = center_y + 1 + (Lattice_Size - coordinate_y)
        return distance_from_center_Y_Coordinate_Exception(center_x, center_y, center_z, coordinate_x, distance_y, coordinate_z)
    else
        return distance_from_center(center_x, center_y, center_z, coordinate_x, coordinate_y, coordinate_z)
    end
end


"""
    Z_Coordiante_Upward(center_x, center_y, center_z, coordinate_x, coordinate_y, coordinate_z) -> Float64

Calculates the distance from the center to a coordinate in the positive Z direction.
Handles boundary crossing using periodic wrapping.

# Returns
- Distance as a Float64.
"""
function Z_Coordiante_Upward(center_x, center_y, center_z, coordinate_x, coordinate_y, coordinate_z)
    if coordinate_z < center_z  # Crossed periodic boundary
        distance_z = (Lattice_Size - center_z) + coordinate_z + 1
        return distance_from_center_Z_Coordinate_Exception(center_x, center_y, center_z, coordinate_x, coordinate_y, distance_z)
    else
        return distance_from_center(center_x, center_y, center_z, coordinate_x, coordinate_y, coordinate_z)
    end
end


"""
    Z_Coordinate_Downward(center_x, center_y, center_z, coordinate_x, coordinate_y, coordinate_z) -> Float64

Calculates the distance from the center to a coordinate in the negative Z direction.
Handles boundary crossing with periodic wraparound logic.

# Returns
- Distance as a Float64.
"""
function Z_Coordinate_Downward(center_x, center_y, center_z, coordinate_x, coordinate_y, coordinate_z)
    if coordinate_z > center_z  # Crossed periodic boundary
        distance_z = center_z + 1 + (Lattice_Size - coordinate_z)
        return distance_from_center_Z_Coordinate_Exception(center_x, center_y, center_z, coordinate_x, coordinate_y, distance_z)
    else
        return distance_from_center(center_x, center_y, center_z, coordinate_x, coordinate_y, coordinate_z)
    end
end


"""
    Change_State_of_Sphere_Coordinates()

Assigns a sphere state (5) and a unique identifier to all coordinates stored in `Possible_Sphere_Coordinates_Set`.
This finalizes the creation of a sphere after its coordinates are determined.
"""
function Change_State_of_Sphere_Coordinates()
    println("Total number of coordinates assigned to sphere: ", length(Possible_Sphere_Coordinates_Set))

    Sphere_Unique_Number = Randomly_Choose_Unique_Number_Sphere()
    println("Sphere assigned unique identifier: $Sphere_Unique_Number")

    for coordinate in keys(Possible_Sphere_Coordinates_Set)
        Assigns_State_Monomer_Sphere(coordinate, Sphere_Unique_Number)
    end
end


"""
    Randomly_Choose_Unique_Number_Sphere() -> Int

Selects and removes a unique identifier from the available sphere number pool.

# Returns
- Unique integer ID for the newly created sphere.
"""
function Randomly_Choose_Unique_Number_Sphere()
    if isempty(Sphere_Unique_Numbers)
        error("No more unique sphere identifiers available!")
    end

    idx = rand(1:length(Sphere_Unique_Numbers))
    unique_code = Sphere_Unique_Numbers[idx]
    deleteat!(Sphere_Unique_Numbers, idx)
    
    return unique_code
end

"""
    Assigns_State_Monomer_Sphere(Location::Tuple, Sphere_Unique_Number::Int)

Assigns the Sphere state (5) and a unique identifier to a specified lattice coordinate.

# Arguments
- `Location`: Tuple of (x, y, z) coordinate.
- `Sphere_Unique_Number`: Unique identifier assigned to this sphere.
"""
function Assigns_State_Monomer_Sphere(Location, Sphere_Unique_Number)
    global Locations_and_States_Dict
    Locations_and_States_Dict[Location] = (SphereState_Value, Sphere_Unique_Number)
end


"""
    Empty_Possible_Sphere_Coordinates()

Clears the temporary dictionary that stores coordinates belonging to a sphere under construction.
"""
function Empty_Possible_Sphere_Coordinates()
    empty!(Possible_Sphere_Coordinates_Set)
end


"""
    Push_Coordinate_Sphere(X::Float64, Y::Float64, Z::Float64)

Adds a coordinate to the `Possible_Sphere_Coordinates_Set` for temporary sphere assembly.

# Arguments
- `X`, `Y`, `Z`: Coordinate components to store.
"""
function Push_Coordinate_Sphere(X, Y, Z)
    Possible_Sphere_Coordinates_Set[(X, Y, Z)] = nothing
end

"""
    determine_growth_direction()

Analyzes and prints the primary growth axis (X, Y, or Z) of the most recently generated sphere,
based on the spread of coordinates in `Possible_Sphere_Coordinates_Set`.
"""
function determine_growth_direction()
    x_values = [coord[1] for coord in keys(Possible_Sphere_Coordinates_Set)]
    y_values = [coord[2] for coord in keys(Possible_Sphere_Coordinates_Set)]
    z_values = [coord[3] for coord in keys(Possible_Sphere_Coordinates_Set)]

    x_range = maximum(x_values) - minimum(x_values)
    y_range = maximum(y_values) - minimum(y_values)
    z_range = maximum(z_values) - minimum(z_values)

    if x_range > y_range && x_range > z_range
        println("Sphere growth is primarily along the X-axis.")
    elseif y_range > x_range && y_range > z_range
        println("Sphere growth is primarily along the Y-axis.")
    elseif z_range > x_range && z_range > y_range
        println("Sphere growth is primarily along the Z-axis.")
    else
        println("Sphere growth is approximately isotropic across axes.")
    end

    println("X range: $x_range, Y range: $y_range, Z range: $z_range")
end


"""
    Calculate_Target_Number_of_Spheres() -> Int

Estimates how many spheres should be generated based on the desired crowding concentration and sphere volume.

# Returns
- Target number of spheres to place.
"""
function Calculate_Target_Number_of_Spheres()
    total_locations = length(keys(Locations_and_States_Dict))
    occupied_spaces = Crowder_Concentration_Spheres * total_locations
    target_occupied_spaces = round(Int, occupied_spaces)

    target_number_of_spheres = round(Int, target_occupied_spaces / Return_Sphere_Volume())

    println("Total lattice locations: $total_locations")
    println("Target occupied spaces based on desired concentration: $target_occupied_spaces")
    println("Estimated number of spheres needed: $target_number_of_spheres")

    return target_number_of_spheres
end


"""
    Count_Number_Coordinates_Spheres() -> Int

Counts the number of lattice coordinates currently assigned to sphere crowders.

# Returns
- Integer number of coordinates assigned to spheres.
"""
function Count_Number_Coordinates_Spheres()
    return count(state == SphereState_Value for (state, _) in values(Locations_and_States_Dict))
end

"""
    Return_Sphere_Volume() -> Int

Returns the number of FCC lattice coordinates that make up a spherical crowder based on the specified obstacle radius.

# Returns
- Integer representing the volume (in number of lattice sites) occupied by a single spherical crowder.
"""

function Return_Sphere_Volume()
    if Obstacle_Radius == 0
        return 1
    elseif Obstacle_Radius == 1
        return 20
    elseif Obstacle_Radius == 2
        return 130
    elseif Obstacle_Radius == 3
        return 400
    elseif Obstacle_Radius == 4
        return 920
    elseif Obstacle_Radius == 5
        return 1850
    elseif Obstacle_Radius == 6
        return 3300
    end
end

#-------------------------------------------------

Generate_Coordinates(Lattice_Size)
Number_of_Coordinates_Sphere = Count_Number_Coordinates_Spheres()
println("The number of coordinates that are spheres: $Number_of_Coordinates_Sphere")
