include { longest_prefix } from './utils.nf'

process filterBED{
    label 'bigmem'
    container "olivierlabayle/tl-core:0.7"
    publishDir "$params.OUTDIR/qc_filtered_chromosomes", mode: 'symlink'

    input:
        path bedfiles
        path qcfile
        path ld_blocks
        path traits

    output:
        path "filtered.*", emit: filtered_bedfiles

    script:
        prefix = bedfiles[0].toString().minus('.bed')
        qc_file = qcfile.getName() != 'NO_QC_FILE' ? "--qcfile $qcfile" : '' 
        ld_blocks = ld_blocks.getName() != 'NO_LD_BLOCKS' ? "--ld-blocks $ld_blocks" : ''
        """
        TEMPD=\$(mktemp -d)
        JULIA_DEPOT_PATH=\$TEMPD:/opt julia --project=/TargeneCore.jl --startup-file=no /TargeneCore.jl/bin/prepare_confounders.jl \
                    --input $prefix --output filtered.$prefix $qc_file --maf-threshold ${params.MAF_THRESHOLD} \
                    $ld_blocks --traits $traits filter
        """
}


process thinByLD{
    label 'bigmem'
    container "olivierlabayle/plink2:0.1.0"
    publishDir "$params.OUTDIR/ld_pruned_chromosomes", mode: 'symlink'

    input:
        path flashpca_excl_reg
        path bedfiles

    output:
        path "LDpruned.*"

    script:
        prefix = bedfiles[0].toString().minus('.bed')
        """
        plink2 --memory ${task.memory.toMega()} --bfile $prefix --indep-pairwise 1000 50 0.05 --exclude range $flashpca_excl_reg
        plink2 --memory ${task.memory.toMega()} --bfile $prefix --extract plink2.prune.in --make-bed --out LDpruned.$prefix
        """
}


process mergeBEDS{
    label 'bigmem'
    container "olivierlabayle/tl-core:0.7"
    publishDir "$params.OUTDIR/merged_genotypes", mode: 'symlink'
    
    input:
        path files
    
    output:
        path "ukbb_merged*"

    script:
        """
        TEMPD=\$(mktemp -d)
        JULIA_DEPOT_PATH=\$TEMPD:/opt julia --project=/TargeneCore.jl --startup-file=no /TargeneCore.jl/bin/prepare_confounders.jl \
        --input LDpruned. \
        --output ukbb_merged merge
        """

}

process LocoMergeBEDS {
    label 'bigmem'
    container "olivierlabayle/tl-core:0.6"

    input:
        tuple val(chr), path(files)
    output:
        tuple path("${chr}_excluded*"), val(chr)

    script:
        prefix = longest_prefix(files)
        //len = get_length(files)

        """
        #!/bin/bash
    
        TEMPD=\$(mktemp -d)
        JULIA_DEPOT_PATH=\$TEMPD:/opt julia --project=/TargeneCore.jl --startup-file=no /TargeneCore.jl/bin/prepare_confounders.jl --input ${prefix} --output ${chr}_excluded merge
        
        """
}

process SampleQCFilter {
    label 'bigmem'
    container "olivierlabayle/plink2:0.1.0"
    publishDir "$params.OUTDIR/iid_genotypes", mode: 'symlink'

    input:
        tuple path(merged_bed_files), val(chr)
    
    output:
        tuple path("qc*"), val(chr)

    script:
        "plink2 --bfile ${chr}_excluded --make-bed --hwe 1e-10 --geno --mind --out qc_exc_${chr}"
}

process FlashPCA {
    label "multithreaded"
    container "ktetleycampbell/flashpca:1.0"

    input:
        tuple path(bedfiles), val(chr)    
    output:
        tuple path("pcs.txt"), val(chr)
    
    script:
        prefix = bedfiles[0].toString().minus('.bed')
        "/home/flashpca-user/flashpca/flashpca --bfile $prefix --ndim ${params.NB_PCS} --numthreads $task.cpus"
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