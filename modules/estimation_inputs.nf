include { longest_prefix } from './utils.nf'

def removeLastPartOfPath(String filePath) {
    // Create a File object from the file path
    def file = new File(filePath)
    
    // Get the parent directory of the file
    def parentDirectory = file.parent
    
    // Return the parent directory
    return parentDirectory
}

process GWASEstimationInputs {
    label "bigmem" // this will need 10 cores so add new label
    label 'targenecore_image'
    publishDir "$params.OUTDIR/estimands", mode: 'symlink', pattern: "*.jls"
    publishDir "$params.OUTDIR", mode: 'symlink', pattern: "*.arrow", saveAs: { filename -> "${params.ARROW_OUTPUT}" }

    input:
        path traits
        tuple val(chr), path(genetic_confounders)
        path parameter
        val command

    output:
        tuple(path("${chr}_final.data.arrow"), path("${chr}_final.*.jls"))

    script:
        batch_size = params.BATCH_SIZE == 0 ? "" :  "--batch-size ${params.BATCH_SIZE}"

        """ 
        TEMPD=\$(mktemp -d)
        JULIA_DEPOT_PATH=\$TEMPD:/opt julia --threads=${task.cpus} --project=/TargeneCore.jl --startup-file=no --sysimage=/TargeneCore.jl/TargeneCoreSysimage.so /TargeneCore.jl/targenecore.jl \
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