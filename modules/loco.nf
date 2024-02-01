#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

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

// NF parameter for PCA step
params.NB_PCS = 6
params.OUTDIR = "$launchDir/test/output"

// Access the genotype process from TarGene
include { mergeBEDS } from "$launchDir/modules/genotypes.nf"
//  PCA processes from TarGene
include { FlashPCA; AdaptFlashPCA } from "$launchDir/modules/confounders.nf"
// Create channel that stores genotype info
input_channel_pairs = Channel.fromFilePairs("$launchDir/test/data/phased_bed/ukb_chr*.{bed,bim,fam}", size: 3)

process MergeBEDS {
    label 'bigmem'
    container "olivierlabayle/tl-core:0.6"
    publishDir "$params.OUTDIR/genotypes", mode: 'symlink'

    input:
        path files
    output:
        path "ukbb_merged*"
    script:
        prefix = longest_prefix(files)
        """
        TEMPD=\$(mktemp -d)
        JULIA_DEPOT_PATH=\$TEMPD:/opt julia --project=/TargeneCore.jl --startup-file=no /TargeneCore.jl/bin/prepare_confounders.jl --input ${prefix} --output ukbb_merged merge
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
        """
        #!/bin/bash
    
        TEMPD=\$(mktemp -d)
        JULIA_DEPOT_PATH=\$TEMPD:/opt julia --project=/TargeneCore.jl --startup-file=no /TargeneCore.jl/bin/prepare_confounders.jl --input ${prefix} --output ${chr}_excluded merge
        
        """
}

workflow RunPCALoco {
    take:
        bed_files
    main:
        merged_bed_ch = LocoMergeBEDS(bed_files)
        pcs_txt = FlashPCA(merged_bed_ch)
        pcs_csv = AdaptFlashPCA(pcs_txt)

    emit:
        pcs_csv
        
}

workflow {
    // Test to extract merged .bed
    genotype_files = Channel.fromPath("$launchDir/test/data/phased_bed/ukb_chr*.{bed,bim,fam}").collect()
    genotype_files.view()
    MergeBEDS(genotype_files)
    
    // DO NOT DELETE
    all_files = input_channel_pairs.map { it[1] }.collect().toList()
    
    input_channel_pairs
        .combine(all_files)
        .map{ [it[0], it[2].findAll{f->!f.normalize().toString().contains(it[0])}] }
        .set{exclusion_set}
    
    RunPCALoco(exclusion_set)

}