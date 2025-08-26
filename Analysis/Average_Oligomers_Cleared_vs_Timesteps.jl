using CSV
using Statistics
using Plots
using DataFrames
using Dates


"""
Module: Oligomer_Clearance_Mean_vs_Timesteps.jl

Computes the **mean number of oligomers cleared per timestep** across multiple simulation
runs and generates a timestamped plot for downstream analysis.

This script reads the aggregated comparison table
`Compare_Simulations/Appending_Oligomers_Clearance_Count.csv`, computes the row-wise mean
(across simulations) for each timestep, and plots **Oligomers Cleared vs Timesteps**.

Implements FAIR principles:

- **Findable**: Inputs/outputs use clear, consistent filenames under `Compare_Simulations/`.
- **Accessible**: Data is handled as CSV and DataFrames; plots saved as PNG.
- **Interoperable**: Uses common Julia packages and file formats.
- **Reusable**: Small functions for import, processing, and plotting.

Key Features:
- Skips missing values when computing per-timestep means.
- Supports a user-defined `Total_Timesteps`.
- Saves a timestamped figure to the `Directory` folder.

Authors: Santiago Schnell; Conner Sandefur; Isabella Gimon  
Dependencies:
- CSV
- DataFrames
- Dates
- Plots
- Statistics

License: http://www.apache.org/licenses/LICENSE-2.0
"""


Mean_Values = Float64[]
Timestep_Values = Int[]


"""
    Import_Data()

Loads the aggregated oligomer-clearance comparison CSV, computes per-timestep means,
and generates/saves the summary plot.

# Behavior
- Constructs path: `Directory * "/Appending_Oligomers_Clearance_Count.csv"`.
- Reads the CSV into a `DataFrame`.
- Builds the `Timestep_Values` vector (`1:Total_Timesteps`).
- Extracts and averages oligomer-clearance values across simulations per timestep.
- Calls `Mean_vs_Timesteps_Graph` to save and display the figure.

# Calls
- `Timesteps`
- `Extract_Data`
- `Mean_vs_Timesteps_Graph`

# Global variables read
- `Directory`
- `Mean_Values`
- `Timestep_Values`
- `Total_Timesteps`
"""

function Import_Data(Directory, Total_Timesteps)
    # Set the working directory to the CSV file
    Working_Directory = Directory * "/Compare_Simulations/Appending_Oligomers_Clearance_Count.csv"
    
    # Read the CSV file into a DataFrame
    data = CSV.read(Working_Directory, DataFrame)
    
    # Process timesteps and extract data
    Timesteps(Total_Timesteps)
    Extract_Data(data)
    
    println("This is Mean_Values: ", Mean_Values)
    println("This is Mean_Values size: ", size(Mean_Values))
    
    # Generate the plot and save it
    Mean_vs_Timesteps_Graph(Directory)
end

"""
    Extract_Data(data::DataFrame)

Computes the mean oligomer-clearance value for each timestep (row-wise),
excluding the first column (`Timesteps`) and skipping missing values.

# Arguments
- `data`: Aggregated comparison table, where column 1 is `Timesteps` and
  columns 2..end are per-simulation oligomer-clearance counts.

# Behavior
- Iterates over rows with `nrow(data)`.
- For each row, selects columns `2:end` and forwards them to `Mean_Per_Timestep`.

# Calls
- `Mean_Per_Timestep`

# Global variables written
- `Mean_Values` (appends a mean per row)
"""

function Extract_Data(data::DataFrame)
    for Row in 1:nrow(data)  # Change here to use nrow(data)
        Row_Data = data[Row, 2:end]
        Mean_Per_Timestep(Row_Data)
    end
end

"""
    Mean_Per_Timestep(Row_Data)

Computes the mean of a single timestep across simulations and appends it to `Mean_Values`.

# Arguments
- `Row_Data`: A row slice (columns 2..end) containing oligomer-clearance values
  for one timestep across simulations.

# Behavior
- Uses `mean(skipmissing(Row_Data))` to ignore missing entries.
- Appends the result to `Mean_Values`.

# Global variables written
- `Mean_Values`
"""

function Mean_Per_Timestep(Row_Data)
    Row_Mean = mean(skipmissing(Row_Data))  # Skip missing values
    push!(Mean_Values, Row_Mean)
end

"""
    Timesteps()

Populates `Timestep_Values` with integers from `1` to `Total_Timesteps`.

# Behavior
- Pushes sequential integers into `Timestep_Values`.

# Global variables read
- `Total_Timesteps`

# Global variables written
- `Timestep_Values`
"""

function Timesteps(Total_Timesteps)
    for Number in 0:Total_Timesteps
        push!(Timestep_Values, Number)
    end
end

"""
    Mean_vs_Timesteps_Graph()

Creates and saves a line plot of **Oligomers Cleared vs Timesteps** using
`Timestep_Values` (x) and `Mean_Values` (y). Saves the figure with a timestamped
filename in `Directory`.

# Behavior
- Labels axes and sets a descriptive title.
- Saves a PNG as `Directory * "/Oligomers_Cleared_vs_Timesteps_<timestamp>.png"`.
- Displays the plot in the current session.

# Global variables read
- `Directory`
- `Timestep_Values`
- `Mean_Values`
"""

function Mean_vs_Timesteps_Graph(Directory)
    # Generate the plot
    Plot = plot(Timestep_Values, Mean_Values, xlabel="Timesteps", ylabel="Total number of oligomers cleared",
    title="Oligomers Cleared vs Timesteps", linewidth=2)

    # Save the plot with a timestamped filename
    timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
    file_name = Directory * "/Compare_Simulations/Oligomers_Cleared_vs_Timesteps_$timestamp.png"
    savefig(Plot, file_name)
    println("Graph saved as: $file_name")
    
    # Display the plot
    display(Plot)
end

