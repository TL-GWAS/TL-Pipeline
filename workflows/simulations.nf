include { EXTRACT_TRAITS } from './traits.nf'
include { PCA } from './pca.nf'
include { EstimationInputs } from '../modules/estimation_inputs.nf'
include { RealisticSimulationInputs; NullSimulationEstimation; RealisticSimulationEstimation; AggregateRealisticSimulationResults; AggregateNullSimulationResults; DensityEstimation } from '../modules/simulations.nf'

workflow NULL_SIMULATION {
    // Workflow specific channels
    estimators = Channel.fromPath(params.ESTIMATORS)
    estimands_files = Channel.value(file("$params.ESTIMANDS_FILE"))
    sample_sizes = Channel.fromList(params.SAMPLE_SIZES)
    rngs = Channel.fromList(params.RNGS)
    bgen_files = Channel.fromPath("$params.BGEN_FILES", checkIfExists: true).collect()

    // PCA
    PCA()

    // generate main dataset and estimand configuration files
    EstimationInputs(
        bgen_files,
        PCA.out.traits,
        PCA.out.pcs,
        estimands_files
    )

    // Null Simulation Estimation
    validated_estimands = EstimationInputs.out.estimands.flatten()
    validated_dataset = EstimationInputs.out.dataset
    bootstrap_grid = estimators.combine(validated_estimands).combine(sample_sizes).combine(rngs)
    simulation_results = NullSimulationEstimation(validated_dataset, bootstrap_grid)

    // Aggregation of Estimation Results
    AggregateNullSimulationResults(simulation_results.collect())
}

workflow REALISTIC_SIMULATION {
    // Workflow specific channels
    ga_trait_table = Channel.value(file(params.GA_TRAIT_TABLE, checkIfExists: true))
    estimators = Channel.fromPath(params.ESTIMATORS)
    estimands_files = Channel.value(file("$params.ESTIMANDS_FILE"))
    sample_sizes = Channel.fromList(params.SAMPLE_SIZES)
    rngs = Channel.fromList(params.RNGS)
    bgen_files = Channel.fromPath("$params.BGEN_FILES", checkIfExists: true).collect()
    
    // PCA
    PCA()

    // Realistic Simulation Inputs
    simulation_inputs = RealisticSimulationInputs(
        estimands_files.collect(),
        bgen_files,
        PCA.out.traits,
        PCA.out.pcs,
        ga_trait_table
    )
    dataset = simulation_inputs.dataset

    // Density Estimation
    density_estimates = DensityEstimation(
        dataset,
        simulation_inputs.conditional_densities.flatten(),
    ).collect()

    // Estimation
    validated_estimands = simulation_inputs.estimands.flatten()
    bootstrap_grid = estimators.combine(validated_estimands).combine(sample_sizes).combine(rngs)
    simulation_results = RealisticSimulationEstimation(
        dataset,
        density_estimates,
        bootstrap_grid
    )

    // Aggregation of Estimation Results
    AggregateRealisticSimulationResults(simulation_results.collect(), dataset, density_estimates)
}