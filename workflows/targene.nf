include { PCA } from './pca.nf'
include { EstimationInputs } from '../modules/estimation_inputs.nf'
include { EstimationWorkflow } from '../subworkflows/estimation.nf'
include { SVPWorkflow } from '../subworkflows/svp.nf'

workflow TARGENE {
    // Define Parameters
    bgen_files = Channel.fromPath("$params.BGEN_FILES", checkIfExists: true).collect()
    estimands_file = Channel.value(file("$params.ESTIMANDS_FILE"))
    estimator_config = Channel.value(file("${params.ESTIMATOR_FILE}"))

    // Extract Traits
    PCA()

    // generate main dataset and estimand configuration files
    EstimationInputs(
        bgen_files,
        PCA.out.traits,
        PCA.out.pcs,
        estimands_file
    )

    // generate estimates
    EstimationWorkflow(
        EstimationInputs.out.dataset,
        EstimationInputs.out.estimands.flatten(),
        estimator_config,
    )

    // Generate sieve variance plateau estimates
    if (params.SVP == true){
        sieve_results = SVPWorkflow(
            EstimationWorkflow.out.hdf5_result, 
            PCA.out.iid_genotypes,
        )
    }
}