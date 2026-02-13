using Test
using FAIR_Implementation_ABM

@testset "FAIR_Implementation_ABM basic checks" begin
    # Timestamp should be filesystem-safe (notably: no ':' for Windows)
    @test !occursin(":", safe_timestamp())

    # Directory creation should work in a temporary location and write provenance files
    mktempdir() do tmp
        run_dir = Make_Directory(output_root=tmp, run_id="TEST_RUN")
        @test isdir(run_dir)
        @test isfile(joinpath(run_dir, "Input_Parameters_used.csv"))
        @test isfile(joinpath(run_dir, "Simulation_Information.csv"))
    end
end
