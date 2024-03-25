include { longest_prefix } from './utils.nf'

process MakeDataset {
    container "olivierlabayle/tl-core:0.8"
    publishDir "${params.OUTDIR}", mode: 'symlink'
    label "bigmem"

    input:
        path bgenfiles
        path traits
        path confounders
        path variants_list
    
    output:
        path "dataset.arrow"
    
    script:
        bgenprefix = longest_prefix(bgenfiles)
        """
        TEMPD=\$(mktemp -d)
        JULIA_DEPOT_PATH=\$TEMPD:/opt julia --project=/TargeneCore.jl --startup-file=no --sysimage=/TargeneCore.jl/TargeneCoreSysimage.so /TargeneCore.jl/bin/generate_dataset.jl \
        ${bgenprefix} ${traits} ${confounders} ${variants_list} \
        --call-threshold ${params.CALL_THRESHOLD} \
        --out dataset.arrow
        """
}