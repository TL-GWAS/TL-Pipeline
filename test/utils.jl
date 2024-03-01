
"""
    check_fails_are_extremely_rare_traits(output, dataset)

The failures are due to traits with only very few cases (<=3) which is due to the fact that 
the dataset is quite degenerate.
"""
function check_fails_are_extremely_rare_traits(results, dataset; ncases=3)
    for result ∈ results
        Ψ = result.estimand
        @test eltype(dataset[!, Ψ.outcome]) == Bool
        nomissing = dropmissing(dataset, Symbol.([Ψ.outcome, keys(Ψ.treatment_values)...]))
        @test sum(nomissing[!, Ψ.outcome]) <= ncases
    end
end

function retrieve_failed_results(results; expected_keys=(:TMLE, :OSE))
    failed_results = (TMLE = [], OSE = [])
    for result ∈ results
        @test keys(result) == expected_keys
        @test result.TMLE isa Union{TMLE.TMLEstimate, TargetedEstimation.FailedEstimate}
        @test result.OSE isa Union{TMLE.OSEstimate, TargetedEstimation.FailedEstimate}
        if result.TMLE isa TargetedEstimation.FailedEstimate
            push!(failed_results.TMLE, result.TMLE)
        end
        if result.OSE isa TargetedEstimation.FailedEstimate
            push!(failed_results.OSE, result.OSE)
        end
    end
    return failed_results
end

function write_custom_configuration()
    config = Configuration(estimands=[
        ATE(
            outcome = "ALL",
            treatment_values = NamedTuple{(Symbol("1:238411180:T:C"), Symbol("3:3502414:T:C"))}([(control = "TT", case = "TC"), (control = "CT", case = "TT")]),
            treatment_confounders = []
        ),
      IATE(
        outcome = "ALL",
        treatment_values = NamedTuple{(Symbol("1:238411180:T:C"), Symbol("3:3502414:T:C"))}([(control = "TT", case = "TC"), (control = "CT", case = "TT")]),
        treatment_confounders = []
      ),
      CM(
        outcome = "ALL",
        treatment_values = NamedTuple{(Symbol("1:238411180:T:C"), Symbol("3:3502414:T:C"))}(["TC", "CT"]),
        treatment_confounders = []
      ),
        ATE(
        outcome = "ALL",
        treatment_values = NamedTuple{(Symbol("2:14983:G:A"),)}([(control = "AG", case = "GG")]),
        treatment_confounders = []
        ),
      CM(
        outcome = "ALL",
        treatment_values = NamedTuple{(Symbol("2:14983:G:A"),)}(["AG"]),
        treatment_confounders = []
      )
    ])
    TMLE.write_yaml(joinpath("test", "assets", "parameters", "parameters.yaml"), config)
end