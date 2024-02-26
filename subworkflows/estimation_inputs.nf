include { TMLEInputsFromActors; TMLEInputsFromParamFile } from '../modules/estimation_inputs.nf'

workflow EstimationInputs {
    take:
        bgen_files
        traits
        genetic_confounders
        estimands_file
        bqtls_file
        transactors_files
        extra_confounders
        extra_treatments
        extra_covariates

    main:
        if (params.STUDY_DESIGN == "FROM_ACTORS") {
            tmle_inputs = TMLEInputsFromActors(
                bgen_files,
                traits,
                genetic_confounders,
                extra_confounders,
                extra_treatments,
                extra_covariates,
                bqtls_file,
                transactors_files,
                )
        }
        else if (params.STUDY_DESIGN == "CUSTOM"){
            tmle_inputs = TMLEInputsFromParamFile(
                bgen_files,
                traits,
                genetic_confounders,
                estimands_file,
                "from-param-file"
                )
        }
        else if (params.STUDY_DESIGN == "ALLELE_INDEPENDENT"){
            tmle_inputs = TMLEInputsFromParamFile(
                bgen_files,
                traits,
                genetic_confounders,
                estimands_file,
                "allele-independent"
                )
        }
        else { 
            throw new Exception("This STUDY_DESIGN is not available.")
        }
    
    emit:
        aggregated_dataset = tmle_inputs.dataset
        estimands = tmle_inputs.estimands
}