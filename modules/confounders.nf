process FlashPCA {
    label "multithreaded"
    container "ktetleycampbell/flashpca:1.0"

    input:
        tuple path(bedfiles), val(chr)    
    output:
        tuple path("pcs.txt"), val(chr)
    
    script:
        prefix = bedfiles[0].toString().minus('.bed')
        "/home/flashpca-user/flashpca/flashpca --bfile $prefix --ndim $params.NB_PCS --numthreads $task.cpus"
}

process AdaptFlashPCA {
    container "olivierlabayle/tl-core:0.6"
    publishDir "$params.OUTDIR/covariates/exc_${chr}", mode: 'symlink'
    label 'bigmem'
    
    input:
        tuple path(flashpca_out), val(chr)
    
    output:
        path "pcs.csv"
    
    script:
        """
        TEMPD=\$(mktemp -d)
        JULIA_DEPOT_PATH=\$TEMPD:/opt julia --project=/TargeneCore.jl --startup-file=no /TargeneCore.jl/bin/prepare_confounders.jl --input $flashpca_out --output pcs.csv adapt
        """
}
