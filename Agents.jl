### Agents.jl — FAIR-Enhanced Version ###

"""
Module: Agents

Simulates the initial state of protein aggregation by initializing a 3D lattice with monomers and optionally crowder spheres.
Adheres to FAIR principles via structured metadata, global state transparency, and traceable logic.

* Findable: Clearly defined structures with traceable global dictionaries.
* Accessible: Self-contained data generation, extensible to export formats.
* Interoperable: Uses standard Julia types (Dict, Tuple, Int, Float64).
* Reusable: Functions are modular, each serving a single purpose.

Author: [Your Lab/Project Name]
Creation Date: [Date]
Last Modified: [Date]
Dependencies:
  - Random, Plots, CSV, DataFrames
Outputs:
  - Global dictionaries with spatial monomer/sphere data
License: http://www.apache.org/licenses/LICENSE-2.0
Date: $(Dates.format(now(), "yyyy-mm-dd"))
"""

using Random
using Plots
using CSV
using DataFrames

# State values: Enumeration of possible states for lattice locations.
const NativeState_Value = 1  # Native monomer state
const AmyloidProne_Value = 2 # Amyloid-prone monomer state
const OligomerState_Value = 3 # Oligomer state
const FibrilState_Value = 4  # Fibril state
const SphereState_Value = 5  # Sphere (crowder) state

# Data structures
# Locations_and_States_Dict: Dictionary mapping 3D coordinates (Tuple{Float64, Float64, Float64}) to a Tuple containing the state (Int) and a unique identifier (Int).
# Possible_Sphere_Coordinates_Set: Set to store coordinates of potential sphere (crowder) locations during generation.
# Initial_Locations_and_States_Dict: Dictionary to store initial locations and states of monomers (Native and Amyloid-prone).
# used_centers: Set to keep track of sphere center coordinates to avoid overlap.
# Sphere_Unique_Numbers: Array of unique identifiers for sphere crowders.
# Monomer_Unique_Numbers: Array of unique identifiers for monomers.
Locations_and_States_Dict = Dict{Tuple{Float64, Float64, Float64}, Tuple{Int, Int}}()
global Possible_Sphere_Coordinates_Set = Dict{Tuple{Float64, Float64, Float64}, Nothing}()
global Initial_Locations_and_States_Dict = Dict{Tuple{Float64, Float64, Float64}, Tuple{Int, Int}}()
global used_centers = Set{Tuple{Float64, Float64, Float64}}()
Sphere_Unique_Numbers = collect(500000:2000000)
Monomer_Unique_Numbers = collect(200001:499999)

# Simulation parameters
const Lattice_Size = 30           # Size of the cubic lattice (Lattice_Size x Lattice_Size x Lattice_Size)
const Max_NumberMonomers_Native = 500   # Maximum number of native monomers
const Max_NumberMonomers_Amyloid = 500  # Maximum number of amyloid-prone monomers
const Obstacle_Radius = 1         # Radius of spherical crowders (if Obstacle is true)
const Crowder_Concentration_Spheres = 0.1 # Concentration of spherical crowders
const Obstacle = false            # Boolean to enable/disable spherical crowders
const Sphere_Volume = 1           # Volume of a single sphere (in lattice units)

"""
Generate_Coordinates(Lattice_Size)

Generates the 3D lattice coordinates and initializes the Locations_and_States_Dict.
If `Obstacle` is true, it also generates spherical crowders. Finally, it randomly
assigns locations to native and amyloid-prone monomers and copies the initial
locations to Initial_Locations_and_States_Dict.

# Arguments
- `Lattice_Size::Int`: The size of the cubic lattice.

# Global variables modified
- `Locations_and_States_Dict`
- `Possible_Sphere_Coordinates_Set` (conditionally)
- `Initial_Locations_and_States_Dict`

# Calls
- `Add_Position`
- `Differentiate_Sphere_Crowder_Radius` (conditionally)
- `Randomly_Assigns_Location_Monomers_Native`
- `Randomly_Assigns_Location_Monomers_Amyloid`
- `Copy_Original_Location`
"""
function Generate_Coordinates(Lattice_Size)
    global Locations_and_States_Dict

    for X in 0:(Lattice_Size - 1)
        for Y in 0:(Lattice_Size - 1)
            for Z in 0:(Lattice_Size - 1)
                State = 0           # 0 indicates an empty location
                Unique_Number = 0   # 0 indicates no monomer/crowder

                # Add corner and face positions for each unit cell
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

Copies the locations and states of native and amyloid-prone monomers from
Locations_and_States_Dict to Initial_Locations_and_States_Dict. This is used
to store the initial conditions of the simulation.

# Global variables modified
- `Initial_Locations_and_States_Dict`

# Assumptions
- `Locations_and_States_Dict` is populated with monomer locations and states.
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

Adds corner and face positions of a unit cell to the Locations_and_States_Dict.

# Arguments
- `X::Int`: X-coordinate of the unit cell.
- `Y::Int`: Y-coordinate of the unit cell.
- `Z::Int`: Z-coordinate of the unit cell.
- `State::Int`: State value to assign to the position (0 for empty).
- `Unique_Number::Int`: Unique identifier to assign to the position (0 for empty).

# Global variables modified
- `Locations_and_States_Dict`
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
First_Corner_Position(X, Y, Z)

Returns the coordinates of the first corner of a unit cell.

# Arguments
- `X::Int`: X-coordinate of the unit cell.
- `Y::Int`: Y-coordinate of the unit cell.
- `Z::Int`: Z-coordinate of the unit cell.

# Returns
- `Tuple{Int, Int, Int}`: The (X, Y, Z) coordinates.
"""
function First_Corner_Position(X, Y, Z)
    return X, Y, Z
end

# (Similar documentation for other coordinate functions: Second_Corner_Position, Third_Corner_Position, etc.)

"""
Randomly_Assigns_Location_Monomers_Native()

Randomly assigns locations to native monomers within the lattice.

# Global variables modified
- `Locations_and_States_Dict`

# Calls
- `Assigns_State_Monomer_Native`
- `Randomly_Choosing_Unique_Number_Monomer`

# Assumptions
- `Locations_and_States_Dict` is initialized.
- `Max_NumberMonomers_Native` defines the maximum number of native monomers.
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

# (Similar documentation for Randomly_Assigns_Location_Monomers_Amyloid, Assigns_State_Monomer_Native, etc.)

"""
Randomly_Choosing_Unique_Number_Monomer()

Randomly selects a unique number from the available list for a monomer.

# Global variables modified
- `Monomer_Unique_Numbers`

# Returns
- `Int`: A unique identifier for the monomer.

# Throws
- `ErrorException`: If there are no more unique numbers available.
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

# (Similar documentation for Differentiate_Sphere_Crowder_Radius, Generate_Spherical_Crowders_Radius_1, etc.)

"""
Assigns_State_Monomer_Native(Location)

Assigns the NativeState_Value and a unique number to a specified location.

# Arguments
- `Location::Tuple{Float64, Float64, Float64}`: The (X, Y, Z) coordinates of the location.

# Global variables modified
- `Locations_and_States_Dict`

# Calls
- `Randomly_Choosing_Unique_Number_Monomer`
"""
function Assigns_State_Monomer_Native(Location)
    global Locations_and_States_Dict
    Unique_Number = Randomly_Choosing_Unique_Number_Monomer()

    # Assign state as Native and keep the unique number
    Locations_and_States_Dict[Location] = (NativeState_Value, Unique_Number)
end

"""
Assigns_State_Monomer_Amyloid(Location)

Assigns the AmyloidProne_Value and a unique number to a specified location.

# Arguments
- `Location::Tuple{Float64, Float64, Float64}`: The (X, Y, Z) coordinates of the location.

# Global variables modified
- `Locations_and_States_Dict`

# Calls
- `Randomly_Choosing_Unique_Number_Monomer`
"""
function Assigns_State_Monomer_Amyloid(Location)
    global Locations_and_States_Dict
    Unique_Number = Randomly_Choosing_Unique_Number_Monomer()
    # Assign state as Amyloid-prone and keep the unique number
    Locations_and_States_Dict[Location] = (AmyloidProne_Value, Unique_Number)
end

"""
Randomly_Choosing_Unique_Number_Monomer()

Randomly selects a unique number from the available list for a monomer.

# Global variables modified
- `Monomer_Unique_Numbers`

# Returns
- `Int`: A unique identifier for the monomer.

# Throws
- `ErrorException`: If there are no more unique numbers available.
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

"""
Differentiate_Sphere_Crowder_Radius()

Determines whether to generate spherical crowders with radius 1 or a different radius.

# Global variables read
- `Obstacle_Radius`

# Calls
- `Generate_Spherical_Crowders_Radius_1` (if Obstacle_Radius == 1)
- `Generate_Spherical_Crowders` (if Obstacle_Radius != 1)
"""
function Differentiate_Sphere_Crowder_Radius()
    if Obstacle_Radius == 1
        Generate_Spherical_Crowders_Radius_1()
    else
        Generate_Spherical_Crowders()
    end
end

"""
Generate_Spherical_Crowders_Radius_1()

Generates spherical crowders with a radius of 1.

# Global variables modified
- `Locations_and_States_Dict`

# Calls
- `Calculate_Target_Number_of_Spheres`
- `Making_Spheres_Radius_1`

# Prints
- Messages indicating sphere creation success/failure.
"""
function Generate_Spherical_Crowders_Radius_1()
    println("We are in the function Generate_Spherical_Crowders_Radius_1")

    # Determine the target number of spheres to create
    target_spheres = Calculate_Target_Number_of_Spheres()
    generated_spheres = 0

    # Generate spheres until we reach the target count
    while generated_spheres < target_spheres
        # Attempt to generate a new sphere
        success = Making_Spheres_Radius_1()

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
Making_Spheres_Radius_1()

Attempts to create a single spherical crowder with a radius of 1.

# Global variables modified
- `Locations_and_States_Dict`

# Calls
- `Randomly_Choose_Unique_Number_Sphere`
- `Assigns_State_Monomer_Sphere`

# Returns
- `Bool`: True if a sphere was created successfully, false otherwise.
"""
function Making_Spheres_Radius_1()
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

Generates spherical crowders with a radius other than 1.

# Global variables modified
- `Locations_and_States_Dict`

# Calls
- `Calculate_Target_Number_of_Spheres`
- `Calling_Sphere_Coordinate_Functions`

# Prints
- Messages indicating sphere creation success/failure.
"""
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
Randomly_Decide_Point()

Randomly selects a coordinate from Locations_and_States_Dict with a state of 0
and that hasn't been used as a sphere center before.

# Global variables read
- `Locations_and_States_Dict`
- `used_centers`

# Global variables modified
- `used_centers`

# Returns
- `Tuple{Float64, Float64, Float64}`: The (centerX, centerY, centerZ) coordinates of the selected point.
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
Calling_Sphere_Coordinate_Functions()

Attempts to create a sphere by calling coordinate calculation functions in each direction.

# Global variables modified
- `Possible_Sphere_Coordinates_Set`
- `Locations_and_States_Dict`

# Calls
- `Randomly_Decide_Point`
- `calculate_sphere_coordinates_left_center_x`
- `calculate_sphere_coordinates_right_center_x`
- `calculate_sphere_coordinates_forward_center_y`
- `calculate_sphere_coordinates_backward_center_y`
- `calculate_sphere_coordinates_upward_center_z`
- `calculate_sphere_coordinates_downward_center_z`
- `determine_growth_direction`
- `Change_State_of_Sphere_Coordinates`
- `Empty_Possible_Sphere_Coordinates`

# Returns
- `Bool`: True if a sphere was created successfully, false otherwise.
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
calculate_sphere_coordinates_left_center_x(center_x, center_y, center_z)

Calculates and stores coordinates to the left of the sphere's center along the x-axis.

# Global variables modified
- `Possible_Sphere_Coordinates_Set`

# Global variables read
- `Locations_and_States_Dict`
- `Obstacle_Radius`

# Calls
- `Push_Coordinate_Sphere`
- `X_Coordinate_Left`
- `Empty_Possible_Sphere_Coordinates`

# Returns
- `Bool`: True if the coordinates are successfully calculated, false if an invalid state is encountered.
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
calculate_sphere_coordinates_right_center_x(center_x, center_y, center_z)

Calculates and stores coordinates to the right of the sphere's center along the x-axis.

# Global variables modified
- `Possible_Sphere_Coordinates_Set`

# Global variables read
- `Locations_and_States_Dict`
- `Obstacle_Radius`

# Calls
- `Push_Coordinate_Sphere`
- `X_Coordinate_Right`
- `Empty_Possible_Sphere_Coordinates`

# Returns
- `Bool`: True if the coordinates are successfully calculated, false if an invalid state is encountered.
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
calculate_sphere_coordinates_forward_center_y(center_x, center_y, center_z)

Calculates and stores coordinates forward of the sphere's center along the y-axis.

# Global variables modified
- `Possible_Sphere_Coordinates_Set`

# Global variables read
- `Locations_and_States_Dict`
- `Obstacle_Radius`

# Calls
- `Push_Coordinate_Sphere`
- `Y_Coordinate_Forward`
- `Empty_Possible_Sphere_Coordinates`

# Returns
- `Bool`: True if the coordinates are successfully calculated, false if an invalid state is encountered.
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
calculate_sphere_coordinates_backward_center_y(center_x, center_y, center_z)

Calculates and stores coordinates backward of the sphere's center along the y-axis.

# Global variables modified
- `Possible_Sphere_Coordinates_Set`

# Global variables read
- `Locations_and_States_Dict`
- `Obstacle_Radius`

# Calls
- `Push_Coordinate_Sphere`
- `Y_Coordinate_Backward`
- `Empty_Possible_Sphere_Coordinates`

# Returns
- `Bool`: True if the coordinates are successfully calculated, false if an invalid state is encountered.
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
calculate_sphere_coordinates_upward_center_z(center_x, center_y, center_z)

Calculates and stores coordinates upward of the sphere's center along the z-axis.

# Global variables modified
- `Possible_Sphere_Coordinates_Set`

# Global variables read
- `Locations_and_States_Dict`
- `Obstacle_Radius`

# Calls
- `Push_Coordinate_Sphere`
- `Z_Coordiante_Upward`
- `Empty_Possible_Sphere_Coordinates`

# Returns
- `Bool`: True if the coordinates are successfully calculated, false if an invalid state is encountered.
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
calculate_sphere_coordinates_downward_center_z(center_x, center_y, center_z)

Calculates and stores coordinates downward of the sphere's center along the z-axis.

# Global variables modified
- `Possible_Sphere_Coordinates_Set`

# Global variables read
- `Locations_and_States_Dict`
- `Obstacle_Radius`

# Calls
- `Push_Coordinate_Sphere`
- `Z_Coordinate_Downward`
- `Empty_Possible_Sphere_Coordinates`

# Returns
- `Bool`: True if the coordinates are successfully calculated, false if an invalid state is encountered.
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

# (Documentation for filter_coordinates functions remains - they are relatively simple)

"""
distance_from_center(center_x, center_y, center_z, coordinate_x, coordinate_y, coordinate_z)

Calculates the Euclidean distance between two points in 3D space.

# Arguments
- `center_x::Float64`: X-coordinate of the center point.
- `center_y::Float64`: Y-coordinate of the center point.
- `center_z::Float64`: Z-coordinate of the center point.
- `coordinate_x::Float64`: X-coordinate of the other point.
- `coordinate_y::Float64`: Y-coordinate of the other point.
- `coordinate_z::Float64`: Z-coordinate of the other point.

# Returns
- `Float64`: The Euclidean distance between the two points.
"""
function distance_from_center(center_x, center_y, center_z, coordinate_x, coordinate_y, coordinate_z)
    distance = sqrt((coordinate_x - center_x)^2 + (coordinate_y - center_y)^2 + (coordinate_z - center_z)^2)
    return distance
end

"""
distance_from_center_X_Coordinate_Exception(center_x, center_y, center_z, distance_x, coordinate_y, coordinate_z)

Calculates the Euclidean distance with a given x-distance, and y and z coordinates.

# Arguments
- `center_x::Float64`: X-coordinate of the center point.
- `center_y::Float64`: Y-coordinate of the center point.
- `center_z::Float64`: Z-coordinate of the center point.
- `distance_x::Float64`: Pre-calculated distance along the x-axis.
- `coordinate_y::Float64`: Y-coordinate of the other point.
- `coordinate_z::Float64`: Z-coordinate of the other point.

# Returns
- `Float64`: The Euclidean distance.
"""
function distance_from_center_X_Coordinate_Exception(center_x, center_y, center_z, distance_x, coordinate_y, coordinate_z)
    distance = sqrt((distance_x)^2 + (coordinate_y - center_y)^2 + (coordinate_z - center_z)^2)
    return distance
end

# (Similar documentation for distance_from_center_Y_Coordinate_Exception, distance_from_center_Z_Coordinate_Exception)

"""
X_Coordinate_Right(center_x, center_y, center_z, coordinate_x, coordinate_y, coordinate_z)

Calculates the x-coordinate for the right direction, considering periodic boundary conditions.

# Arguments
- `center_x::Float64`: X-coordinate of the center point.
- `center_y::Float64`: Y-coordinate of the center point.
- `center_z::Float64`: Z-coordinate of the center point.
- `coordinate_x::Float64`: X-coordinate of the other point.
- `coordinate_y::Float64`: Y-coordinate of the other point.
- `coordinate_z::Float64`: Z-coordinate of the other point.

# Calls
- `distance_from_center`
- `distance_from_center_X_Coordinate_Exception`

# Returns
- `Float64`: The calculated x-coordinate.
"""
function X_Coordinate_Right(center_x, center_y, center_z, coordinate_x, coordinate_y, coordinate_z)
    if coordinate_x < center_x #MAKE X COORDINATE OUT OF BOUNDS
        distance = (Lattice_Size - center_x) + coordinate_x + 1
        return distance_from_center_X_Coordinate_Exception(center_x, center_y, center_z, distance, coordinate_y, coordinate_z)
    else
        return distance_from_center(center_x, center_y, center_z, coordinate_x, coordinate_y, coordinate_z)
    end
end

# (Similar documentation for X_Coordinate_Left, Y_Coordinate_Forward, Y_Coordinate_Backward, Z_Coordiante_Upward, Z_Coordinate_Downward)

"""
Change_State_of_Sphere_Coordinates()

Changes the state of all coordinates in Possible_Sphere_Coordinates_Set to SphereState_Value.

# Global variables modified
- `Locations_and_States_Dict`
- `Sphere_Unique_Numbers`

# Global variables read
- `Possible_Sphere_Coordinates_Set`

# Calls
- `Randomly_Choose_Unique_Number_Sphere`
- `Assigns_State_Monomer_Sphere`
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
Randomly_Choose_Unique_Number_Sphere()

Randomly selects a unique number from the available list for a sphere.

# Global variables modified
- `Sphere_Unique_Numbers`

# Returns
- `Int`: A unique identifier for the sphere.

# Throws
- `ErrorException`: If there are no more unique numbers available.
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
Assigns_State_Monomer_Sphere(Location, Sphere_Unique_Number)

Assigns the SphereState_Value and a unique number to a specified location.

# Arguments
- `Location::Tuple{Float64, Float64, Float64}`: The (X, Y, Z) coordinates of the location.
- `Sphere_Unique_Number::Int`: The unique identifier for the sphere.

# Global variables modified
- `Locations_and_States_Dict`
"""
function Assigns_State_Monomer_Sphere(Location, Sphere_Unique_Number)
    global Locations_and_States_Dict
    _, Unique_Number = Locations_and_States_Dict[Location]
    Locations_and_States_Dict[Location] = (SphereState_Value, Sphere_Unique_Number)
end

"""
Empty_Possible_Sphere_Coordinates()

Clears the Possible_Sphere_Coordinates_Set.

# Global variables modified
- `Possible_Sphere_Coordinates_Set`
"""
function Empty_Possible_Sphere_Coordinates()
    empty!(Possible_Sphere_Coordinates_Set)
end

"""
Push_Coordinate_Sphere(X, Y, Z)

Adds a coordinate to the Possible_Sphere_Coordinates_Set.

# Global variables modified
- `Possible_Sphere_Coordinates_Set`

# Arguments
- `X::Float64`: X-coordinate of the location.
- `Y::Float64`: Y-coordinate of the location.
- `Z::Float64`: Z-coordinate of the location.
"""
function Push_Coordinate_Sphere(X, Y, Z)
    Possible_Sphere_Coordinates_Set[(X, Y, Z)] = nothing
end

"""
determine_growth_direction()

Determines the primary growth direction of the sphere.

# Global variables read
- `Possible_Sphere_Coordinates_Set`

# Prints
- Message indicating the primary growth direction of the sphere.
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
Calculate_Target_Number_of_Spheres()

Calculates the target number of spheres to generate based on the crowder concentration.

# Global variables read
- `Locations_and_States_Dict`
- `Crowder_Concentration_Spheres`
- `Sphere_Volume`

# Returns
- `Int`: The target number of spheres.

# Prints
- Information about total lattice locations, target occupied spaces, and estimated number of spheres.
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
Count_Number_Coordinates_Spheres()

Counts the number of coordinates with the SphereState_Value.

# Global variables read
- `Locations_and_States_Dict`

# Returns
- `Int`: The number of coordinates with the SphereState_Value.
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
println("This is the number of coordinates that are spheres made: ", Count_Number_Coordinates_Spheres())
#println(Calculate_Target_Number_of_Spheres())
#print(Locations_and_States_Dict)
Total_Number_Locations = length(keys(Locations_and_States_Dict))
