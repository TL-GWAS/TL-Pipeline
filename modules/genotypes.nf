process filterBED{
    label 'bigmem'
    container "olivierlabayle/tl-core:0.6"
    publishDir "$params.OUTDIR/qc_filtered_chromosomes", mode: 'symlink'

    input:
        tuple val(chr_id), file(bedfiles)
        path qcfile
        path ld_blocks
        path traits

    output:
        path "filtered.*", emit: filtered_bedfiles

    script:
        prefix = bedfiles[0].toString().minus('.bed')
        qc_file = params.QC_FILE !== 'NO_QC_FILE' ? "--qcfile $qcfile" : '' 
        ld_blocks = params.LD_BLOCKS !== 'NO_LD_BLOCKS' ? "--ld-blocks $ld_blocks" : ''

        """
        TEMPD=\$(mktemp -d)
        JULIA_DEPOT_PATH=\$TEMPD:/opt julia --project=/TargeneCore.jl --startup-file=no /TargeneCore.jl/bin/prepare_confounders.jl \
                    --input $prefix --output filtered.$prefix $qc_file --maf-threshold $params.MAF_THRESHOLD \
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
    container "olivierlabayle/tl-core:0.6"
    publishDir "$params.OUTDIR/merged_genotypes", mode: 'symlink'
    
    input:
        path files
    
    output:
        path "ukbb_merged*"

    script:
        """
        TEMPD=\$(mktemp -d)
        JULIA_DEPOT_PATH=\$TEMPD:/opt julia --project=/TargeneCore.jl --startup-file=no /TargeneCore.jl/bin/prepare_confounders.jl --input LDpruned. --output ukbb_merged merge
        """

}

// Changed for LOCO
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


workflow IIDGenotypes{
    take:
        flashpca_excl_reg
        ld_blocks
        bed_files
        qc_file
        sample_ids

    main:
        filtered_bedfiles = filterBED(bed_files, qc_file, ld_blocks, sample_ids)
        ld_pruned = thinByLD(flashpca_excl_reg, filtered_bedfiles)
        mergeBEDS(ld_pruned.collect())
        SampleQCFilter(mergeBEDS.out.collect())

    emit:
        SampleQCFilter.out
}
