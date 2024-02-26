module TestUKBAlleleIndependent

using Test
using JLD2
using TMLE
using TargetedEstimation
using Serialization

# "local" profile assumes singularity is installed
args = length(ARGS) > 0 ? ARGS : ["-profile", "local", "-resume"] 

@testset "Test ukb_allele_independent.config" begin
    cmd = `nextflow run main.nf -c test/configs/ukb_allele_independent.config $args`
    @info string("The following command will be run:\n", cmd)

    r = run(cmd)
    @test r.exitcode == 0

    ## Checking main output
    # Results
    results_from_hdf5 = jldopen(joinpath("results", "results.hdf5"))["Batch_1"]
    nresults = length(results_from_hdf5)
    @test nresults > 50
    nfails = 0
    treatment_combinations = Set([])
    for result in results_from_hdf5
        Ψ̂ = result.TMLE
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
    # QQ plot
    @test isfile("results/QQ.png")
end

end