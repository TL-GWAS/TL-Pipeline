#!/usr/bin/env nextflow
nextflow.enable.dsl = 2
include { ExtractTraits } from "../subworkflows/extract_traits.nf"
include { LOCOGenotypes; LOCOConfounders } from "../subworkflows/confounders.nf"
include { EstimationInputs } from "../subworkflows/estimation_inputs.nf"
include { EstimationWorkflow } from "../subworkflows/estimation.nf"

workflow LOCOGWAS{
    // Define Parameters
    bed_files = Channel.fromPath("$params.BED_FILES", checkIfExists: true).collect()
    bgen_files = Channel.fromPath("$params.BGEN_FILES", checkIfExists: false).collect()
    estimands_file = Channel.value(file("$params.ESTIMANDS_FILE"))
    bqtls_file = Channel.value(file("$params.BQTLS"))
    transactors_files = Channel.fromPath("$params.TRANS_ACTORS").collect()
    extra_confounders = Channel.value(file("$params.EXTRA_CONFOUNDERS"))
    extra_treatments = Channel.value(file("$params.ENVIRONMENTALS"))
    extra_covariates = Channel.value(file("$params.EXTRA_COVARIATES"))

    traits_dataset = Channel.value(file("$params.TRAITS_DATASET", checkIfExists: true))
    traits_config = Channel.value(file("$params.UKB_CONFIG", checkIfExists: true))
    ukb_withdrawal_list = Channel.value(file("$params.UKB_WITHDRAWAL_LIST", checkIfExists:false))
    ukb_encoding_file = params.UKB_ENCODING_FILE

    loco_bed_files = Channel.fromFilePairs("$params.BED_FILES", size: 3, checkIfExists: true)
    qc_file = Channel.value(file("$params.QC_FILE", checkIfExists: true))
    ld_blocks = Channel.value(file("$params.LD_BLOCKS", checkIfExists: false))

    estimator_config = Channel.value(file("$params.ESTIMATOR_FILE"))

    // Specify root directory where bed files are stored
    root = Channel.fromPath("$params.ROOT")

    if (params.COHORT == "UKBB") {
        ExtractTraits(
            traits_dataset, 
            traits_config, 
            ukb_withdrawal_list, 
            ukb_encoding_file)

        LOCOGenotypes(
            loco_bed_files, 
            ExtractTraits.out,
            qc_file,
            ld_blocks)

        LOCOConfounders(LOCOGenotypes.out)
        
        // Merge bed files with confounder output on chr_id
        loco_files = LOCOConfounders.out.cross(loco_bed_files).map{[it[0][0], it[0][1]]}
    
        EstimationInputs(
            bgen_files,
            ExtractTraits.out,
            LOCOConfounders.out,
            estimands_file,
            bqtls_file,
            transactors_files,
            extra_confounders,
            extra_treatments,
            extra_covariates,
            loco_files,
            root
        )
    
        EstimationInputs.out.transpose().set{inputs}

        //generate estimates
        EstimationWorkflow(
            inputs.map{it[0]}, // dataset
            inputs.map{it[1]},  // estimands
            estimator_config,
        )
    
        // Generate sieve estimates
        if (params.SVP == true){
            sieve_results = SVPWorkflow(
                EstimationWorkflow.out.hdf5_result, 
                IIDGenotypes.out,
            )
        }
    }

    
}