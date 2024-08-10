
const TRAITS = [
  "J40-J47 Chronic lower respiratory diseases",
  "O20-O29 Other maternal disorders predominantly related to pregnancy",
  "bacterial infection",
  "other fractures",
  "D41 Neoplasm of uncertain or unknown behaviour of urinary organs",
  "C34 Malignant neoplasm of bronchus and lung",
  "C50-C50 Malignant neoplasm of breast",
  "Cheese intake",
  "Number of vehicles in household",
  "Skin colour",
  "Pork intake"
]

const PCS = [
  "PC1",
  "PC2",
  "PC3",
  "PC4",
  "PC5",
  "PC6"
]

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

function retrieve_failed_results(results; expected_keys=(:TMLE_GLM_GLM, :OSE_GLM_GLM))
    failed_results = (TMLE_GLM_GLM = [], OSE_GLM_GLM = [])
    for result ∈ results
        @test keys(result) == expected_keys
        @test result.TMLE_GLM_GLM isa Union{TMLE.TMLEstimate, TargetedEstimation.FailedEstimate}
        @test result.OSE_GLM_GLM isa Union{TMLE.OSEstimate, TargetedEstimation.FailedEstimate}
        if result.TMLE_GLM_GLM isa TargetedEstimation.FailedEstimate
            push!(failed_results.TMLE_GLM_GLM, result.TMLE_GLM_GLM)
        end
        if result.OSE_GLM_GLM isa TargetedEstimation.FailedEstimate
            push!(failed_results.OSE_GLM_GLM, result.OSE_GLM_GLM)
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
    TMLE.write_yaml(joinpath("test", "assets", "estimands.yaml"), config)
end

function save_custom_configuration()
  configuration = TMLE.Configuration(estimands = [
    factorialEstimand(
        ATE,
        (;(Symbol("3:3502414:T:C")=>["TT", "TC", "CC"],)...),
        "ALL",
        confounders = ["Number of vehicles in household"],
        outcome_extra_covariates = ["Skin colour"]
        ),
    IATE(
        outcome="other fractures",
        treatment_values=(;(Symbol("3:3502414:T:C")=>(case="TT", control="CT"), Symbol("1:238411180:T:C")=>(case="TC", control="TT"))...),
        treatment_confounders = (;(Symbol("3:3502414:T:C") => ["Number of vehicles in household"], Symbol("1:238411180:T:C")=>[])...)
    )
  ]
  )
  TMLE.write_yaml("test/assets/simulation_estimands.yaml", configuration)
end