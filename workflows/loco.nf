#!/usr/bin/env nextflow
nextflow.enable.dsl = 2
include { ExtractTraits } from "../subworkflows/extract_traits.nf"
include {LOCOGenotypes; LOCOConfounders } from "../subworkflows/confounders.nf"

workflow LOCOGWAS{
    traits_dataset = Channel.value(file("$params.TRAITS_DATASET", checkIfExists: true))
    traits_config = Channel.value(file("$params.TRAITS_CONFIG", checkIfExists: true))
    ukb_withdrawal_list = Channel.value(file("$params.UKB_WITHDRAWAL_LIST", checkIfExists:false))
    ukb_encoding_file = params.UKB_ENCODING_FILE

    loco_bed_files = Channel.fromFilePairs("$params.BED_FILES", size: 3, checkIfExists: true)
    qc_file = Channel.value(file("$params.QC_FILE", checkIfExists: true))
    ld_blocks = Channel.value(file("$params.LD_BLOCKS", checkIfExists: false))
    
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
        )

        // generate estimates
        EstimationWorkflow(
            EstimationInputs.out.aggregated_dataset,
            EstimationInputs.out.estimands.flatten(),
            estimator_config,
        )

        // Generate sieve estimates
        if (params.SVP == true){
            sieve_results = SVPWorkflow(
                EstimationWorkflow.out.hdf5_result, 
                IIDGenotypes.out,
            )
        }

    } else {
        LOCOGenotypes(
            loco_bed_files, 
            traits_dataset,
            qc_file,
            ld_blocks)
    }
    
    

    
}