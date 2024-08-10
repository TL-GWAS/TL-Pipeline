#!/usr/bin/env nextflow
nextflow.enable.dsl = 2
include { PCA } from './pca.nf'
include { GWASEstimationInputs } from '../modules/estimation_inputs.nf'
include { EstimationWorkflow } from '../subworkflows/estimation.nf'
include { SVPWorkflow } from '../subworkflows/svp.nf'

workflow LOCOGWAS {
    // Define Parameters
    bed_files = Channel.fromPath("$params.BED_FILES", checkIfExists: true).collect()
    estimands_file = Channel.value(file("$params.ESTIMANDS_FILE"))
    estimator_config = Channel.value(file("$params.ESTIMATOR_FILE"))

    PCA()

    GWASEstimationInputs(
        PCA.out.traits,
        PCA.out.confounders,
        estimands_file
    )

    GWASEstimationInputs.out.transpose().set{inputs}

    // generate estimates
    EstimationWorkflow(
        inputs.map{it[0]}, // dataset
        inputs.map{it[1]},  // estimands
        estimator_config,
    )

    // Generate sieve estimates
    if (params.SVP == true){
        sieve_results = SVPWorkflow(
            EstimationWorkflow.out.hdf5_result, 
            PCA.out.iid_genotypes,
        )
    }
}