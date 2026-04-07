using CSV
using DataFrames
using Statistics
using Plots
using Measures

# Resolve data directory relative to this script's location.
const SCRIPT_DIR = @__DIR__
const DATA_DIR   = joinpath(SCRIPT_DIR, "..", "data", "Figure_4")
const OUT_DIR    = SCRIPT_DIR

# Load the two CSV files
df_no_clearance = CSV.read(
    joinpath(DATA_DIR, "No_Oligomer_Clearance",
             "Appending_Aggregate_Count_No_Oligomer_Clearance.csv"),
    DataFrame)
df_clearance = CSV.read(
    joinpath(DATA_DIR, "Oligomer_Clearance",
             "Appending_Aggregate_Count_Oligomer_Clearance.csv"),
    DataFrame)

# Keep only timesteps up to 300
df_no_clearance = filter(row -> row.Timesteps <= 300, df_no_clearance)
df_clearance    = filter(row -> row.Timesteps <= 300, df_clearance)

# Row-wise mean across all simulation columns
function row_averages(df::DataFrame)
    sim_matrix = Matrix(df[:, 2:end])
    return vec(mean(sim_matrix, dims=2))
end

timesteps        = df_no_clearance.Timesteps
avg_no_clearance = row_averages(df_no_clearance)
avg_clearance    = row_averages(df_clearance)

plot(
    timesteps, avg_no_clearance,
    label      = "No Oligomer Clearance",
    xlabel     = "Timestep",
    ylabel     = "Average Aggregate Count (n=5)",
    lw         = 1.8,
    xlims      = (0, 300),
    legend     = :bottomright,
    guidefont  = font(13),
    tickfont   = font(10),
    legendfont = font(9),
    size       = (900, 600),
    dpi        = 300,
    left_margin = 2mm
)
plot!(timesteps, avg_clearance,
      label = "Oligomer Clearance (p=0.05)",
      lw    = 1.8)

savefig(joinpath(OUT_DIR, "Figure4.png"))
println("Saved Figure4.png to ", OUT_DIR)
