using Random
using Plots
using CSV
using DataFrames

# State values
NativeState_Value = 1
AmyloidProne_Value = 2
OligomerState_Value = 3
FibrilState_Value = 4
SphereState_Value = 5

# Declare the dictionary to store locations, states, and unique numbers
Locations_and_States_Dict = Dict{Tuple{Float64, Float64, Float64}, Tuple{Int, Int}}()
global Possible_Sphere_Coordinates_Set = Dict{Tuple{Float64, Float64, Float64}, Nothing}()
global Initial_Locations_and_States_Dict = Dict{Tuple{Float64, Float64, Float64}, Tuple{Int, Int}}()
global used_centers = Set{Tuple{Float64, Float64, Float64}}()
Sphere_Unique_Numbers = collect(500000:2000000)
Monomer_Unique_Numbers = collect(200001:499999)

# Assign individual variables
Lattice_Size = 30
Max_NumberMonomers_Native = 500
Max_NumberMonomers_Amyloid = 500
Obstacle_Radius = 1
Crowder_Concentration_Spheres = .1
Obstacle = false
Sphere_Volume = 1


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

function Copy_Original_Location()
    global Initial_Locations_and_States_Dict  # Ensure it modifies the global variable

    for (location, (state, unique_number)) in Locations_and_States_Dict
        if state == NativeState_Value || state == AmyloidProne_Value
            Initial_Locations_and_States_Dict[location] = (state, unique_number)
        end
    end
end


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

function First_Corner_Position(X, Y, Z)
    return X, Y, Z
end

function Second_Corner_Position(X, Y, Z)
    return X + 1, Y, Z
end

function Third_Corner_Position(X, Y, Z)
    return X, Y + 1, Z
end

function Fourth_Corner_Position(X, Y, Z)
    return X, Y, Z + 1
end

function Fifth_Corner_Position(X, Y, Z)
    return X + 1, Y + 1, Z
end

function Sixth_Corner_Position(X, Y, Z)
    return X + 1, Y, Z + 1
end

function Seventh_Corner_Position(X, Y, Z)
    return X, Y + 1, Z + 1
end

function Eighth_Corner_Position(X, Y, Z)
    return X + 1, Y + 1, Z + 1
end

function First_Face_Position(X, Y, Z) #Face-centered along each axis
    return X + 0.5, Y + .5, Z
end

function Second_Face_Position(X, Y, Z) #Face-centered along each axis
    return X, Y + 0.5, Z + .5
end
 
function Third_Face_Position(X, Y, Z) #Face-centered along each axis
    return X + .5, Y, Z + 0.5
end

function Fourth_Face_Position(X, Y, Z) #Opposite Faces
    return X + .5, Y + 0.5, Z + 1
end

function Fifth_Face_Position(X, Y, Z)  #Opposite Faces
    return X + 0.5, Y + 1, Z + .5
end

function Sixth_Face_Position(X, Y, Z)  #Opposite Faces
    return X + 1, Y + .5, Z + 0.5
end

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

# Functions to assign states to specific coordinates
function Assigns_State_Monomer_Native(Location)
    global Locations_and_States_Dict
    Unique_Number = Randomly_Choosing_Unique_Number_Monomer()
    
    # Assign state as Native and keep the unique number
    Locations_and_States_Dict[Location] = (NativeState_Value, Unique_Number)
end

function Assigns_State_Monomer_Amyloid(Location)
    global Locations_and_States_Dict
    Unique_Number = Randomly_Choosing_Unique_Number_Monomer()
    # Assign state as Native and keep the unique number
    Locations_and_States_Dict[Location] = (AmyloidProne_Value, Unique_Number)
end

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
function Differentiate_Sphere_Crowder_Radius()
    if Obstacle_Radius == 1 
        Generate_Spherical_Crowders_Radius_1()
    else
        Generate_Spherical_Crowders()
    end

end

function Generate_Spherical_Crowders_Radius_1()
    println("We are in the funciton Generate_Spherical_Crowders_Radius_1")

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

function Making_Spheres_Radius_1()
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


# LEFT OF CENTER_X
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

# RIGHT OF CENTER_X
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

# FORWARD OF CENTER_Y
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

# BACKWARD OF CENTER_Y
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

# UPWARD OF CENTER_Z
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

# DOWNWARD OF CENTER_Z
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




#0. Filter coordinats to the left of center_x
function filter_coordinates_left_of_center(center_x)
    # Return only coordinates with x-value less than center_x
    return filter(coord -> coord[1] <= center_x, keys(Locations_and_States_Dict))
end

# 1. Filter coordinates to the right of center_x
function filter_coordinates_right_of_center(center_x)
    return filter(coord -> coord[1] > center_x, keys(Locations_and_States_Dict))
end

# 2. Filter coordinates forward of center_y
function filter_coordinates_forward_of_center(center_y)
    return filter(coord -> coord[2] > center_y, keys(Locations_and_States_Dict))
end

# 3. Filter coordinates backward of center_y
function filter_coordinates_backward_of_center(center_y)
    return filter(coord -> coord[2] <= center_y, keys(Locations_and_States_Dict))
end

# 4. Filter coordinates upward of center_z
function filter_coordinates_upward_of_center(center_z)
    return filter(coord -> coord[3] > center_z, keys(Locations_and_States_Dict))
end

# 5. Filter coordinates downward of center_z
function filter_coordinates_downward_of_center(center_z)
    return filter(coord -> coord[3] <= center_z, keys(Locations_and_States_Dict))
end

function distance_from_center(center_x, center_y, center_z, coordinate_x, coordinate_y, coordinate_z)
    distance = sqrt((coordinate_x - center_x)^2 + (coordinate_y - center_y)^2 + (coordinate_z - center_z)^2)
    return distance
end

function distance_from_center_X_Coordinate_Exception(center_x, center_y, center_z, distance_x, coordinate_y, coordinate_z)
    distance = sqrt((distance_x)^2 + (coordinate_y - center_y)^2 + (coordinate_z - center_z)^2)
    return distance
end

function distance_from_center_Y_Coordinate_Exception(center_x, center_y, center_z, coordinate_x, distance_y, coordinate_z)
    distance = sqrt((coordinate_x - center_x)^2 + (distance_y)^2 + (coordinate_z - center_z)^2)
    return distance
end

function distance_from_center_Z_Coordinate_Exception(center_x, center_y, center_z, coordinate_x, coordinate_y, distance_z)
    distance = sqrt((coordinate_x - center_x)^2 + (coordinate_y - center_y)^2 + (distance_z)^2)
    return distance
end

function X_Coordinate_Right(center_x, center_y, center_z, coordinate_x, coordinate_y, coordinate_z)
    if coordinate_x < center_x #MAKE X COORDINATE OUT OF BOUNDS
        distance = (Lattice_Size - center_x) + coordinate_x + 1
        distance_from_center_X_Coordinate_Exception(center_x, center_y, center_z, distance, coordinate_y, coordinate_z)
    else
        distance_from_center(center_x, center_y, center_z, coordinate_x, coordinate_y, coordinate_z)
    end
end

function X_Coordinate_Left(center_x, center_y, center_z, coordinate_x, coordinate_y, coordinate_z)
    if coordinate_x > center_x #MAKE X COORDINATE OUT OF BOUNDS
        distance = center_x + 1 + (Lattice_Size - coordinate_x)
        #println("For coordinate_x: $coordinate_x the new distance is: $distance")
        distance_from_center_X_Coordinate_Exception(center_x, center_y, center_z, distance, coordinate_y, coordinate_z)
    else
        distance_from_center(center_x, center_y, center_z, coordinate_x, coordinate_y, coordinate_z)
    end 
end

function Y_Coordinate_Forward(center_x, center_y, center_z, coordinate_x, coordinate_y, coordinate_z)
    if coordinate_y < center_y #MAKE Y COORDINATE OUT OF BOUNDS
        distance = (Lattice_Size - center_y) + coordinate_y + 1
        distance_from_center_Y_Coordinate_Exception(center_x, center_y, center_z, coordinate_x, distance, coordinate_z)
    else
        distance_from_center(center_x, center_y, center_z, coordinate_x, coordinate_y, coordinate_z)
    end

end

function Y_Coordinate_Backward(center_x, center_y, center_z, coordinate_x, coordinate_y, coordinate_z)
    if coordinate_y > center_y #MAKE Y COORDINATE OUT OF BOUNDS
        distance = center_y + 1 + (Lattice_Size - coordinate_y)
        distance_from_center_Y_Coordinate_Exception(center_x, center_y, center_z, coordinate_x, distance, coordinate_z)
    else
        distance_from_center(center_x, center_y, center_z, coordinate_x, coordinate_y, coordinate_z)
    end
    
end

function Z_Coordiante_Upward(center_x, center_y, center_z, coordinate_x, coordinate_y, coordinate_z)
    if coordinate_z < center_z #MAKE Z COORDINATE OUT OF BOUNDS
        distance = (Lattice_Size - center_z) + coordinate_z + 1
        distance_from_center_Z_Coordinate_Exception(center_x, center_y, center_z, coordinate_x, coordinate_y, distance)
    else
        distance_from_center(center_x, center_y, center_z, coordinate_x, coordinate_y, coordinate_z)
    end
end

function Z_Coordinate_Downward(center_x, center_y, center_z, coordinate_x, coordinate_y, coordinate_z)
    if coordinate_z > center_z #MAKE Z COORDINATE OUT OF BOUNDS
        distance =  center_z + 1 + (Lattice_Size - coordinate_z)
        distance_from_center_Z_Coordinate_Exception(center_x, center_y, center_z, coordinate_x, coordinate_y, distance)
    else
        distance_from_center(center_x, center_y, center_z, coordinate_x, coordinate_y, coordinate_z)
    end

end



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

function Assigns_State_Monomer_Sphere(Location, Sphere_Unique_Number)
    global Locations_and_States_Dict
    _, Unique_Number = Locations_and_States_Dict[Location]
    Locations_and_States_Dict[Location] = (SphereState_Value, Sphere_Unique_Number)
end

function Empty_Possible_Sphere_Coordinates()
    empty!(Possible_Sphere_Coordinates_Set)  
end

function Push_Coordinate_Sphere(X, Y, Z)
    Possible_Sphere_Coordinates_Set[(X, Y, Z)] = nothing
end


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
println("This is the number of coordinates that are spheres made: ",Count_Number_Coordinates_Spheres())
#println(Calculate_Target_Number_of_Spheres())
#print(Locations_and_States_Dict)
Total_Number_Locations = length(keys(Locations_and_States_Dict))