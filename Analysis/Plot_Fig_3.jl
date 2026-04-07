using CSV
using DataFrames
using Statistics
using Plots
using Measures

# Resolve data directory relative to this script's location.
const SCRIPT_DIR = @__DIR__
const DATA_DIR   = joinpath(SCRIPT_DIR, "..", "data", "Figure_3")
const OUT_DIR    = SCRIPT_DIR

# Load the four CSV files
df_aggprone  = CSV.read(joinpath(DATA_DIR, "Appending_AggregateProne_Count.csv"),  DataFrame)
df_native    = CSV.read(joinpath(DATA_DIR, "Appending_Native_Count.csv"),           DataFrame)
df_oligomer  = CSV.read(joinpath(DATA_DIR, "Appending_Oligomer_Count.csv"),         DataFrame)
df_aggregate = CSV.read(joinpath(DATA_DIR, "Appending_Aggregate_Count.csv"),        DataFrame)

# Keep only timesteps up to 300
df_aggprone  = filter(row -> row.Timesteps <= 300, df_aggprone)
df_native    = filter(row -> row.Timesteps <= 300, df_native)
df_oligomer  = filter(row -> row.Timesteps <= 300, df_oligomer)
df_aggregate = filter(row -> row.Timesteps <= 300, df_aggregate)

# Row-wise mean across all simulation columns
function row_averages(df::DataFrame)
    sim_matrix = Matrix(df[:, 2:end])
    return vec(mean(sim_matrix, dims=2))
end

timesteps     = df_aggprone.Timesteps
avg_aggprone  = row_averages(df_aggprone)
avg_native    = row_averages(df_native)
avg_oligomer  = row_averages(df_oligomer)
avg_aggregate = row_averages(df_aggregate)

plot(
    timesteps, avg_native,
    label      = "Native",
    xlabel     = "Timestep",
    ylabel     = "Average count (n=5)",
    lw         = 1.8,
    xlims      = (0, 300),
    legend     = :topright,
    guidefont  = font(13),
    tickfont   = font(10),
    legendfont = font(10),
    size       = (900, 600),
    dpi        = 300,
    left_margin = 5mm
)
plot!(timesteps, avg_aggprone,  label = "AggregateProne", lw = 1.8)
plot!(timesteps, avg_oligomer,  label = "Oligomer",       lw = 1.8)
plot!(timesteps, avg_aggregate, label = "Aggregate",      lw = 1.8)

savefig(joinpath(OUT_DIR, "Figure3.png"))
println("Saved Figure3.png to ", OUT_DIR)
