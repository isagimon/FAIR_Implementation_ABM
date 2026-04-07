using Test
using FAIR_Implementation_ABM

@testset "FAIR_Implementation_ABM QA checks" begin
    # Timestamp should be filesystem-safe (notably: no ':' for Windows)
    @test !occursin(":", safe_timestamp())

    @testset "Parameter validation" begin
        params = Dict{String,Any}(
            "Lattice_Size" => 3,
            "MAX_NumberMovements" => 10,
            "Max_NumberMonomers_Native" => 15,
            "Max_NumberMonomers_AggregateProne" => 15,
            "Spheres?" => false,
            "Obstacle_Radius" => 0,
            "Crowder_Concentration_Spheres" => 0.0,
            "Native_to_AggregateProne" => 0.2,
            "AggregateProne_to_Native" => 0.2,
            "Oligomer_Formation" => 0.05,
            "Oligomer_Dissociation_rate" => 0.005,
            "Fibril_Formation" => 0.1,
            "Fibril_Growth" => 0.9,
            "Probability_of_Oligomer_Removal" => 0.0
        )

        @test FAIR_Implementation_ABM.validate_parameters!(copy(params)) isa Dict

        bad = copy(params)
        delete!(bad, "Lattice_Size")
        @test_throws ErrorException FAIR_Implementation_ABM.validate_parameters!(bad)

        bad2 = copy(params)
        bad2["Native_to_AggregateProne"] = 2.0
        @test_throws ErrorException FAIR_Implementation_ABM.validate_parameters!(bad2)

        bad3 = copy(params)
        bad3["Obstacle_Radius"] = 7
        @test_throws ErrorException FAIR_Implementation_ABM.validate_parameters!(bad3)
    end

    @testset "FAIR_Implementation_ABM.Return_Sphere_Volume guard" begin
        old = FAIR_Implementation_ABM.Obstacle_Radius
        try
            FAIR_Implementation_ABM.Obstacle_Radius = 7
            @test_throws ErrorException FAIR_Implementation_ABM.Return_Sphere_Volume()
        finally
            FAIR_Implementation_ABM.Obstacle_Radius = old
        end
    end

    @testset "Output directory precedence + provenance copy" begin
        mktempdir() do tmp
            # Minimal parameter CSV used only for provenance copying
            paramfile = joinpath(tmp, "params.csv")
            open(paramfile, "w") do io
                write(io, "Input_Parameters,Values,Instructions\n")
                write(io, "Directory,Data_Collection,\n")
            end

            run_dir = Make_Directory(output_root=tmp, run_id="TEST_RUN", parameter_file=paramfile)
            @test isdir(run_dir)
            @test isfile(joinpath(run_dir, "Simulation_Information.csv"))
            @test isfile(joinpath(run_dir, "Input_Parameters_used.csv"))
        end
    end
end
