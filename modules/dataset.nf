include { longest_prefix } from './utils.nf'

process MakeDataset {
    label 'targenecore_image'
    publishDir "${params.OUTDIR}", mode: 'symlink'
    label "bigmem"

    input:
        path bgenfiles
        path traits
        path pcs_file
        path variants_file
    
    output:
        path "dataset.arrow"
    
    script:
        bgenprefix = longest_prefix(bgenfiles)
        call_threshold = params.CALL_THRESHOLD == null ? "" : "--call-threshold ${params.CALL_THRESHOLD}"
        """
        TEMPD=\$(mktemp -d)
        JULIA_DEPOT_PATH=\$TEMPD:/opt julia --project=/TargeneCore.jl --startup-file=no --sysimage=/TargeneCore.jl/TargeneCoreSysimage.so /TargeneCore.jl/targenecore.jl \
            make-dataset \
            --genotypes-prefix=${bgenprefix} \
            --traits-file=${traits} \
            --pcs-file=${pcs_file} \
            --variants-file=${variants_file} \
            --out=dataset.arrow \
            ${call_threshold} \
            --verbosity=${params.VERBOSITY}
        """
}