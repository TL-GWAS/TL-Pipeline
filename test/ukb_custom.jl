module TestFromParamFile

using Test
using DataFrames
using TMLE
using Arrow
using JLD2
using TargetedEstimation

# "local" profile assumes singularity is installed
args = length(ARGS) > 0 ? ARGS : ["-profile", "local", "-resume"]

include("utils.jl")

@testset "Test ukb_from_param_files.config" begin
    cmd = `nextflow run main.nf -c test/configs/ukb_custom.config $args`
    @info string("The following command will be run:\n", cmd)

    r = run(cmd)
    @test r.exitcode == 0
    
    result_file = jldopen(joinpath("results", "results.hdf5"))
    results = vcat(result_file["Batch_1"], result_file["Batch_2"])
    dataset = DataFrame(Arrow.Table(joinpath("results", "dataset.arrow")))

    failed_results = retrieve_failed_results(results)

    # All fails are due to fluctuation failure due non positive definite matrix
    # This does not affect the OSE
    @test isempty(failed_results.OSE)
    # Less than 1/3 affected: this is still quite significant
    @test length(failed_results.TMLE) / length(results) < 1/3
    @test all(startswith(x.msg, "Could not fluctuate") for x ∈ failed_results.TMLE)

    check_fails_are_extremely_rare_traits(failed_results.TMLE, dataset; ncases=3)
end

end

true