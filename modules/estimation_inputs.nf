include { longest_prefix } from './utils.nf'

def removeLastPartOfPath(String filePath) {
    // Create a File object from the file path
    def file = new File(filePath)
    
    // Get the parent directory of the file
    def parentDirectory = file.parent
    
    // Return the parent directory
    return parentDirectory
}

process TMLEInputsGWAS {
    container "olivierlabayle/tl-core:loco-gwas"
    publishDir "$params.OUTDIR/estimands", mode: 'symlink', pattern: "*.jls"
    publishDir "$params.OUTDIR", mode: 'symlink', pattern: "*.arrow", saveAs: { filename -> "${params.ARROW_OUTPUT}" }
    label "bigmem"

    input:
        path traits
        tuple val(chr), path(genetic_confounders)
        path parameter
        val command
        path root
    output:
        tuple(path("${chr}_final.data.arrow"), path("${chr}_final.*.jls"))

    script:
        batch_size = params.BATCH_SIZE == 0 ? "" :  "--batch-size ${params.BATCH_SIZE}"

        """ 
        dir=\$( echo "$params.BED_FILES" | rev | cut -f2 -d'/' - | rev)
        bed_prefix=\$( echo "\$dir/$chr")
        TEMPD=\$(mktemp -d)
        JULIA_DEPOT_PATH=\$TEMPD:/opt julia --project=/TargeneCore.jl --startup-file=no --sysimage=/TargeneCore.jl/TargeneCoreSysimage.so /TargeneCore.jl/bin/generate_tl_inputs.jl \
        --positivity-constraint ${params.POSITIVITY_CONSTRAINT} \
        $batch_size \
        --out-prefix=${chr}_final \
        --verbosity=${params.VERBOSITY} \
        $command $parameter\
        --traits $traits \
        --genotype-prefix \$bed_prefix \
        --pcs $genetic_confounders \
        """
}

process TMLEInputsFromParamFile {
    container "olivierlabayle/tl-core:loco-gwas"
    publishDir "$params.OUTDIR/estimands", mode: 'symlink', pattern: "*.jls"
    publishDir "$params.OUTDIR", mode: 'symlink', pattern: "*.arrow", saveAs: { filename -> "${params.ARROW_OUTPUT}" }
    label "bigmem"

    input:
        path bgenfiles
        path traits
        path genetic_confounders
        path parameter
        val command

    output:
        path "final.data.arrow", emit: dataset
        path "final.*.jls", emit: estimands

    script:
        bgen_prefix = longest_prefix(bgenfiles)
        batch_size = params.BATCH_SIZE == 0 ? "" :  "--batch-size ${params.BATCH_SIZE}"
        """
        TEMPD=\$(mktemp -d)
        JULIA_DEPOT_PATH=\$TEMPD:/opt julia --project=/TargeneCore.jl --startup-file=no --sysimage=/TargeneCore.jl/TargeneCoreSysimage.so /TargeneCore.jl/bin/generate_tl_inputs.jl \
        --positivity-constraint ${params.POSITIVITY_CONSTRAINT} \
        $batch_size \
        --out-prefix=final \
        --verbosity=${params.VERBOSITY} \
        $command $parameter \
        --traits $traits \
        --bgen-prefix $bgen_prefix \
        --call-threshold ${params.CALL_THRESHOLD} \
        --pcs $genetic_confounders \
        """
}

process TMLEInputsFromActors {
    container "olivierlabayle/tl-core:loco-gwas"
    publishDir "$params.OUTDIR/estimands", mode: 'symlink', pattern: "*.jls"
    publishDir "$params.OUTDIR", mode: 'symlink', pattern: "*.arrow", saveAs: { filename -> "${params.ARROW_OUTPUT}" }
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
        path "final.data.arrow", emit: dataset
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
        JULIA_DEPOT_PATH=\$TEMPD:/opt julia --project=/TargeneCore.jl --startup-file=no --sysimage=/TargeneCore.jl/TargeneCoreSysimage.so /TargeneCore.jl/bin/generate_tl_inputs.jl \
        --positivity-constraint ${params.POSITIVITY_CONSTRAINT} \
        $batch_size \
        --out-prefix=final \
        --verbosity=${params.VERBOSITY} \
        from-actors $bqtls $trans_actors_prefix \
        $extra_confounders \
        $extra_treatments \
        $extra_covariates \
        --orders ${params.ORDERS} \
        --traits $traits \
        --bgen-prefix $bgen_prefix \
        --call-threshold ${params.CALL_THRESHOLD} \
        --pcs $genetic_confounders \
        """
}