def longest_prefix(files){
    // Only one file, strangely it is not passed as a list
    if (files instanceof Collection == false) {
        return files.getName()
    }
    // More than one file
    index = 0
    while(true){
        current_prefix = files[0].getName()[0..index]
        for (file in files){
            if(file.getName()[0..index] != current_prefix){
                return current_prefix[0..-2]
            }
        }
        index++
    }
}

process TMLE {
    container "olivierlabayle/targeted-estimation:cv_tmle"
    publishDir "$params.OUTDIR/csvs",  mode: 'symlink', pattern: "*.csv"
    publishDir "$params.OUTDIR/hdf5files/inf_curves",  mode: 'symlink', pattern: "*.hdf5"
    label "bigmem"
    label "multithreaded"

    input:
        path data
        path estimandsfile
        path estimatorfile
    
    output:
        path "${hdf5out}"
    
    script:
        basename = "result." + estimandsfile.getName().take(estimandsfile.getName().lastIndexOf('.'))
        hdf5out = basename + ".hdf5"
        pval_threshold = KEEP_IC == true ? "--outputs.hdf5.pval_threshold=${params.PVAL_THRESHOLD}" : ""
        sample_ids = SVP == true ? "--outputs.hdf5.sample_ids=true" : ""
        """
        TEMPD=\$(mktemp -d)
        JULIA_DEPOT_PATH=\$TEMPD:/opt julia --project=/TargetedEstimation.jl --threads=${task.cpus} --startup-file=no tmle tmle \
        $data \
        --estimands=$estimandsfile \
        --estimators=$estimatorfile \
        --outputs.hdf5.filename=$hdf5out \
        $pval_threshold \
        $sample_ids \
        --chunksize=$params.TMLE_SAVE_EVERY \
        """
}

process TMLEInputsFromParamFile {
    container "olivierlabayle/tl-core:cvtmle"
    publishDir "$params.OUTDIR/estimands", mode: 'symlink', pattern: "*.yaml"
    publishDir "$params.OUTDIR", mode: 'symlink', pattern: "*.arrow"
    label "bigmem"

    input:
        path bgenfiles
        path traits
        path genetic_confounders
        path parameter

    output:
        path "${params.ARROW_OUTPUT}", emit: traits
        path "final.*.jls", emit: estimands

    script:
        bgen_prefix = longest_prefix(bgenfiles)
        batch_size = params.BATCH_SIZE == 0 ? "" : "--batch-size ${params.BATCH_SIZE}"
        """
        TEMPD=\$(mktemp -d)
        JULIA_DEPOT_PATH=\$TEMPD:/opt julia --project=/TargeneCore.jl --startup-file=no /TargeneCore.jl/bin/tmle_inputs.jl \
        --traits $traits \
        --bgen-prefix $bgen_prefix \
        --call-threshold ${params.CALL_THRESHOLD} \
        --pcs $genetic_confounders \
        $batch_size \
        --positivity-constraint ${params.POSITIVITY_CONSTRAINT} \
        from-param-file $parameter
        """
}

process TMLEInputsFromActors {
    container "olivierlabayle/tl-core:cvtmle"
    publishDir "$params.OUTDIR/estimands", mode: 'symlink', pattern: "*.yaml"
    publishDir "$params.OUTDIR", mode: 'symlink', pattern: "*.arrow"
    label "bigmem"

    input:
        path bgenfiles
        path traits
        path genetic_confounders
        path extra_confounders
        path extra_treatments
        path extra_covariates
        path bqtls
        path trans_actors

    output:
        path "${params.ARROW_OUTPUT}", emit: traits
        path "final.*.jls", emit: estimands

    script:
        bgen_prefix = longest_prefix(bgenfiles)
        trans_actors_prefix = longest_prefix(trans_actors)
        batch_size = params.BATCH_SIZE == 0 ? "" :  "--batch-size ${params.BATCH_SIZE}"
        extra_confounders = extra_confounders.name != 'NO_EXTRA_CONFOUNDER' ? "--extra-confounders $extra_confounders" : ''
        extra_treatments = extra_treatments.name != 'NO_EXTRA_TREATMENT' ? "--extra-treatments $extra_treatments" : ''
        extra_covariates = extra_covariates.name != 'NO_EXTRA_COVARIATE' ? "--extra-covariates $extra_covariates" : ''
        """
        TEMPD=\$(mktemp -d)
        JULIA_DEPOT_PATH=\$TEMPD:/opt julia --project=/TargeneCore.jl --startup-file=no /TargeneCore.jl/bin/tmle_inputs.jl \
        --traits $traits \
        --bgen-prefix $bgen_prefix \
        --call-threshold ${params.CALL_THRESHOLD} \
        --pcs $genetic_confounders \
        --positivity-constraint ${params.POSITIVITY_CONSTRAINT} \
        $batch_size \
        from-actors $bqtls $trans_actors_prefix $extra_confounders $extra_treatments $extra_covariates --orders ${params.ORDERS}
        """
}
