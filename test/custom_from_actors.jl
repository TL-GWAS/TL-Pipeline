module TestCustomDatasetFromActors

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

@testset "Test custom_from_actors.config" begin
    cmd = `nextflow run main.nf -c test/configs/custom_from_actors.config $args`
    @info string("The following command will be run:\n", cmd)

    r = run(cmd)
    @test r.exitcode == 0

    # Sieve Variance Plateau File
    svp = jldopen(joinpath("results", "svp.hdf5"))
    @test length(svp["taus"]) == 10
    @test haskey(svp, "variances")
    @test length(svp["results"]) > 200
    @test all(x.TMLE isa TMLE.TMLEstimate for x in svp["results"])

    # Results
    results_file = jldopen(joinpath("results", "results.hdf5"))
    results = vcat((results_file[key] for key in keys(results_file))...)
    @test length(results) > 300
    failed_results = retrieve_failed_results(results; expected_keys=(:TMLE, :OSE, :SAMPLE_IDS))
    # All fails are due to fluctuation failure due non positive definite matrix
    # This does not affect the OSE
    @test isempty(failed_results.OSE)
    # Less than 1/3 affected: this is still quite significant
    @test length(failed_results.TMLE) / length(results) < 1/3
    @test all(startswith(x.msg, "Could not fluctuate") for x ∈ failed_results.TMLE)

    dataset = Arrow.Table(joinpath("results", "dataset.arrow")) |> DataFrame

    check_fails_are_extremely_rare_traits(failed_results.TMLE, dataset; ncases=7)

    # Here we test that the process IIDGenotypes has been run once for each chromosome
    # There should be 4 files for each chromosome
    n_chr_files = filter(x -> startswith(x, "LDpruned.filtered.ukb_chr"), readdir(joinpath("results", "ld_pruned_chromosomes")))
    @test length(n_chr_files) == 4*3

    ## Checking parameter files correspond to either bQTL only or bQTL/eQTL
    bQTLs = Symbol.(CSV.read(joinpath("test", "assets", "actors", "bqtls.csv"), DataFrame).ID)

    estimands_tf1 = deserialize(joinpath("results", "estimands", "final.TF1.estimands_1.jls"))
    estimands_tf2 = deserialize(joinpath("results", "estimands", "final.TF2.estimands_1.jls"))
    @test size(estimands_tf1.estimands, 1) == 296
    @test size(estimands_tf2.estimands, 1) == 148
    for Ψ in vcat(estimands_tf1.estimands, estimands_tf2.estimands)
        @test length(intersect(keys(Ψ.treatment_values), bQTLs)) == 1
    end
end

end

true