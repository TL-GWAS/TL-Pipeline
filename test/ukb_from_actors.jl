module TestUKBFromActors

using Test
using JLD2
using Arrow
using DataFrames
using TMLE
using CSV
using TargetedEstimation
using Serialization

# "local" profile assumes singularity is installed
args = length(ARGS) > 0 ? ARGS : ["-profile", "local", "-resume"] 

include("utils.jl")

@testset "Test ukb_from_actors.config" begin
    cmd = `nextflow run main.nf -c test/configs/ukb_from_actors.config $args`
    @info string("The following command will be run:\n", cmd)

    r = run(cmd)
    @test r.exitcode == 0

    ## Checking main output
    # Results
    hdf5_results_file = jldopen(joinpath("results", "results.hdf5"))
    results_from_hdf5 = vcat((hdf5_results_file[key] for key in keys(hdf5_results_file))...)
    results_from_json = TMLE.read_json(joinpath("results", "results.json"))
    @test length(results_from_json) == length(results_from_hdf5) > 300

    failed_results = retrieve_failed_results(results_from_hdf5; expected_keys=(:TMLE, :OSE))
    # All fails are due to fluctuation failure due non positive definite matrix
    # This does not affect the OSE
    @test isempty(failed_results.OSE)
    # Less than 1/3 affected: this is still quite significant
    @test length(failed_results.TMLE) / length(results_from_hdf5) < 1/3
    @test all(startswith(x.msg, "Could not fluctuate") for x ∈ failed_results.TMLE)

    dataset = Arrow.Table(joinpath("results", "dataset.arrow")) |> DataFrame

    check_fails_are_extremely_rare_traits(failed_results.TMLE, dataset; ncases=3)

    ## Checking parameter files correspond to either bQTL only or bQTL/eQTL
    bQTLs = Symbol.(CSV.read(joinpath("test", "assets", "actors", "bqtls.csv"), DataFrame).ID)

    config_1 = deserialize(joinpath("results", "estimands", "final.estimands_1.jls"))
    @test length(config_1.estimands) == 400
    config_2 = deserialize(joinpath("results", "estimands", "final.estimands_2.jls"))
    @test 0 < length(config_2.estimands) <= 400
    for Ψ in vcat(config_1.estimands, config_2.estimands)
        @test length(intersect(keys(Ψ.treatment_values), bQTLs)) == 1
    end
end

@testset "Test PERMUTATION_TEST" begin
    cmd = `nextflow run main.nf -entry PERMUTATION_TEST -c test/configs/negcontrol.config $args`
    @info string("The following command will be run:\n", cmd)

    r = run(cmd)
    @test r.exitcode == 0

    data = Arrow.Table(joinpath("results", "permutation_tests", "permutation_dataset.arrow")) |> DataFrame
    
    n_permuted_cols = 0
    for colname in names(data)
        if endswith(colname, "permuted")
            n_permuted_cols +=1
        end
    end
    @test n_permuted_cols > 20

    permutation_config = deserialize(joinpath("results", "permutation_tests", "estimands", "permutation_estimands_1.jls"))
    @test length(permutation_config.estimands) == 100
    @test Set(typeof(Ψ) for Ψ in permutation_config.estimands) == Set([TMLE.StatisticalATE, TMLE.StatisticalIATE])
    
    results = jldopen(joinpath("results", "permutation_results.hdf5"))["Batch_1"]
    @test length(results) == 100
    failed_results = retrieve_failed_results(results; expected_keys=(:TMLE, :OSE))
    @test failed_results == (TMLE=[], OSE=[])
end

@testset "Test RANDOMIZATION_TEST" begin
    cmd = `nextflow run main.nf -entry RANDOMIZATION_TEST -c test/configs/negcontrol.config $args`
    @info string("The following command will be run:\n", cmd)

    r = run(cmd)
    @test r.exitcode == 0

    random_config = deserialize(joinpath("results", "random_variants_estimands.jls"))
    @test length(random_config.estimands) > 200
    @test Set(typeof(Ψ) for Ψ in random_config.estimands) == Set([TMLE.StatisticalATE, TMLE.StatisticalIATE])
end

end
true