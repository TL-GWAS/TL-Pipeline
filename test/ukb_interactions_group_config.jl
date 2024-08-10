module TestUKBAlleleIndependent

using Test
using JLD2
using TMLE
using TargetedEstimation
using Serialization

# "local" profile assumes singularity is installed
args = length(ARGS) > 0 ? ARGS : ["-profile", "local", "-resume"] 

include("utils.jl")

@testset "Test ukb_interactions_group.config" begin
    cmd = `nextflow run main.nf -c test/configs/ukb_interactions_group.config $args`
    @info string("The following command will be run:\n", cmd)

    r = run(cmd)
    @test r.exitcode == 0

    ## Checking main output
    # Results
    results_from_hdf5 = jldopen(io -> io["results"], joinpath("results", "results.hdf5"))
    nresults = length(results_from_hdf5)
    @test nresults > 20
    nfails = 0
    treatment_combinations = Set([])
    for result in results_from_hdf5
        Ψ̂ = result.TMLE_GLM_GLM
        if Ψ̂ isa TargetedEstimation.FailedEstimate
            nfails += 1
        else
            Ψc = first(Ψ̂.estimand.args)
            push!(treatment_combinations, keys(Ψc.treatment_values))
        end
    end
    @test nfails / nresults < 1/4
    @test treatment_combinations ==  Set([ 
        (Symbol("2:14983:G:A"), Symbol("3:3502414:T:C")),
        (Symbol("1:238411180:T:C"), Symbol("2:14983:G:A"))
    ])

    # Dataset
    dataset = TargetedEstimation.instantiate_dataset("results/dataset.arrow")
    @test Set(names(dataset)) == Set(vcat("SAMPLE_ID", TRAITS, PCS, ["2:14983:G:A", "3:3502414:T:C", "1:238411180:T:C"]))

    # QQ plot
    @test isfile("results/QQ.png")

     # Check properly resumed
    resume_time = @elapsed run(cmd)
    @test resume_time < 1000
end

end