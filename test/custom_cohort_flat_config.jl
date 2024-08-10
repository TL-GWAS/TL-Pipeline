module TestCustomCohortFlatConfig

using Test
using Arrow
using DataFrames
using JLD2
using TargetedEstimation
using CSV

args = length(ARGS) > 0 ? ARGS : ["-profile", "local", "-resume"] 

@testset "Test custom_cohort_flat.config" begin
    cmd = `nextflow run main.nf -c test/configs/custom_cohort_flat.config $args`
    @info string("The following command will be run:\n", cmd)

    r = run(cmd)
    @test r.exitcode == 0

    dataset_origin = CSV.read(joinpath("test", "assets", "traits.csv"), DataFrame)
    dataset = DataFrame(Arrow.Table(joinpath("results", "dataset.arrow")))
    @test issubset(names(dataset_origin), names(dataset))
    
    results = jldopen(io -> io["results"], "results/results.hdf5")
    @test length(results) > 100

    svp_results = jldopen(io -> io["results"], "results/svp.hdf5")
    @test length(svp_results) > 100

    # Check properly resumed
    resume_time = @elapsed run(cmd)
    @test resume_time < 1000
end

end

true