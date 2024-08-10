include { longest_prefix } from './utils.nf'

process GWASEstimationInputs {
    publishDir "$params.OUTDIR/estimands", mode: 'symlink', pattern: "*.jls"
    publishDir "$params.OUTDIR", mode: 'symlink', pattern: "*.arrow", saveAs: { filename -> "${params.ARROW_OUTPUT}" }
    label "largemem"
    label 'targenecore_image'

    input:
        path traits
        tuple val(chr), path(genetic_confounders), path(bed), path(bim), path(fam)
        path config_file

    output:
        tuple(path("${chr}_final.data.arrow"), path("${chr}_final.*.jls"))

    script:
        batch_size = params.BATCH_SIZE == 0 ? "" :  "--batchsize ${params.BATCH_SIZE}"
        call_threshold = params.CALL_THRESHOLD == null ? "" : "--call-threshold ${params.CALL_THRESHOLD}"
        """ 
        TEMPD=\$(mktemp -d)
        JULIA_DEPOT_PATH=\$TEMPD:/opt julia --threads=${task.cpus} --project=/TargeneCore.jl --startup-file=no --sysimage=/TargeneCore.jl/TargeneCoreSysimage.so /TargeneCore.jl/targenecore.jl \
        estimation-inputs ${config_file} \
         --genotypes-prefix=${chr} \
        --traits-file=${traits} \
        --pcs-file=${genetic_confounders} \
        --outprefix=${chr}_final \
        ${batch_size} \
        ${call_threshold} \
        --positivity-constraint=${params.POSITIVITY_CONSTRAINT} \
        --verbosity=${params.VERBOSITY}
        """
}

process EstimationInputs {
    publishDir "$params.OUTDIR/estimands", mode: 'symlink', pattern: "*.jls"
    publishDir "$params.OUTDIR", mode: 'symlink', pattern: "*.arrow", saveAs: { filename -> "${params.ARROW_OUTPUT}" }
    label "bigmem"
    label 'targenecore_image'

    input:
        path genotypes_prefix
        path traits
        path genetic_confounders
        path config_file

    output:
        path "final.data.arrow", emit: dataset
        path "final.*.jls", emit: estimands

    script:
        genotypes_prefix = longest_prefix(genotypes_prefix)
        batch_size = params.BATCH_SIZE == 0 ? "" :  "--batchsize ${params.BATCH_SIZE}"
        call_threshold = params.CALL_THRESHOLD == null ? "" : "--call-threshold ${params.CALL_THRESHOLD}"
        """
        TEMPD=\$(mktemp -d)
        JULIA_DEPOT_PATH=\$TEMPD:/opt julia --project=/TargeneCore.jl --startup-file=no --sysimage=/TargeneCore.jl/TargeneCoreSysimage.so /TargeneCore.jl/targenecore.jl \
        estimation-inputs ${config_file} \
        --genotypes-prefix=${genotypes_prefix} \
        --traits-file=${traits} \
        --pcs-file=${genetic_confounders} \
        --outprefix=final \
        ${batch_size} \
        ${call_threshold} \
        --positivity-constraint=${params.POSITIVITY_CONSTRAINT} \
        --verbosity=${params.VERBOSITY}
        """
}