using Test
using FAIR_Implementation_ABM

const FABM = FAIR_Implementation_ABM

@testset "FAIR_Implementation_ABM" begin

    @testset "Basic utilities" begin
        @test !occursin(":", FABM.safe_timestamp())
    end

    @testset "Movement dispatch consistency" begin
        for move in FABM.Possible_Movement_Options
            if move != "None"
                @test haskey(FABM.MovementFunctions, move)
            end
        end
    end

    @testset "Output directory precedence" begin
        repo_root = abspath(joinpath(@__DIR__, ".."))
        data_root = abspath(joinpath(repo_root, "Data_Collection"))

        old_env = get(ENV, "FAIR_ABM_OUTPUT_DIR", nothing)
        old_parameters = copy(FABM.Parameters)

        try
            ENV["FAIR_ABM_OUTPUT_DIR"] = "env_out"
            @test FABM.Resolve_Output_Directory(nothing, repo_root, data_root) == abspath("env_out")

            pop!(ENV, "FAIR_ABM_OUTPUT_DIR", nothing)
            @test FABM.Resolve_Output_Directory("arg_out", repo_root, data_root) == abspath("arg_out")

            empty!(FABM.Parameters)
            FABM.Parameters["Directory"] = "csv_out"
            @test FABM.Resolve_Output_Directory(nothing, repo_root, data_root) == abspath("csv_out")

            empty!(FABM.Parameters)
            @test FABM.Resolve_Output_Directory(nothing, repo_root, data_root) == abspath(data_root)
        finally
            empty!(FABM.Parameters)
            merge!(FABM.Parameters, old_parameters)

            if old_env === nothing
                pop!(ENV, "FAIR_ABM_OUTPUT_DIR", nothing)
            else
                ENV["FAIR_ABM_OUTPUT_DIR"] = old_env
            end
        end
    end

    @testset "Per-run provenance files" begin
        mktempdir() do tmp
            run_dir = FABM.Make_Directory(
                output_root=tmp,
                run_id="TEST_RUN",
                parameter_file=FABM.PARAMETER_FILE
            )

            @test isdir(run_dir)
            @test isfile(joinpath(run_dir, "Simulation_Information.csv"))
            @test isfile(joinpath(run_dir, "Input_Parameters_used.csv"))
        end
    end

    @testset "Return_Sphere_Volume invalid radius" begin
        old_radius = FABM.Obstacle_Radius
        try
            FABM.Obstacle_Radius = 7
            @test_throws Exception FABM.Return_Sphere_Volume()
        finally
            FABM.Obstacle_Radius = old_radius
        end
    end

end