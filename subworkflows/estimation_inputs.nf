include { TMLEInputsFromActors; TMLEInputsFromParamFile; TMLEInputsFromLOCOGWAS } from '../modules/estimation_inputs.nf'

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
        loco_files


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
        else if (params.STUDY_DESIGN == "LOCO_GWAS"){
            tmle_inputs = TMLEInputsFromLOCOGWAS(
                traits,
                loco_files,
                estimands_file,
                "loco-gwas",
            )
        }
        else { 
            throw new Exception("This STUDY_DESIGN is not available.")
        }
    
    emit:
        tmle_inputs
        // aggregated_dataset = tmle_inputs.dataset
        // estimands = tmle_inputs.estimands
}