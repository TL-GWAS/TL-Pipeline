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

// parameter for root directory (where all genotype data is stored)
params.ROOT = "$launchDir/test/data/phased_bed"

// olvier merge
include { mergeBEDS } from "$launchDir/modules/genotypes.nf"
//  PCA processes
include { FlashPCA; AdaptFlashPCA } from "$launchDir/modules/confounders.nf"

// Create channel that stores genotype info
input_channel_pairs = Channel.fromFilePairs("$launchDir/test/data/phased_bed/ukb_chr*.{bed,bim,fam}", size: 3)
// Creats channel that collects all files in root directory and stores them in a list
root_channel = Channel.fromPath("$params.ROOT")

process LocoMergeBEDS {
    label 'bigmem'
    container "olivierlabayle/tl-core:0.6"
    input:
        path files
        val chr
    output:
        path "${chr}_excluded*"

    script:
        
        prefix = longest_prefix(files)
        """
        #!/bin/bash
    
        TEMPD=\$(mktemp -d)
        JULIA_DEPOT_PATH=\$TEMPD:/opt julia --project=/TargeneCore.jl --startup-file=no /TargeneCore.jl/bin/prepare_confounders.jl --input ${prefix} --output ${chr}_excluded merge
        
        """
}

// Write new subworkflow that will take that chr channel and will run the this workflow once per
// chromosome
// workflow Loco {
//     take:
//         chr_channel
//         exc_chr
//         files

//     chr_channel
//         .filter{it != exc_chr}
//         .join(input_channel_pairs)
//         .map{prefix, files -> files}
//         .collect()
//         .set{merged_ch}
//     exc_chr_ch = Channel.of(exc_chr) // channel to manipulate to add multiple emissions
//     LocoMergeBEDS(merged_ch, exc_chr_ch)
// }

workflow {
    
    // obtain the chromosome numbers and associated files
    chr_channel = input_channel_pairs.map { it[0] }
    // chr_channel.view()

    exc_chr = 'ukb_chr1'
    chr_channel
        .filter{it != exc_chr}
        .join(input_channel_pairs)
        .map{prefix, files -> files}
        .collect()
        .set{merged_ch}
    exc_chr_ch = Channel.of(exc_chr) // channel to manipulate to add multiple emissions
    merged_bed_ch = LocoMergeBEDS(merged_ch, exc_chr_ch)
    pcs_txt = FlashPCA(merged_bed_ch)
    pcs_csv = AdaptFlashPCA(pcs_txt)
    pcs_csv.view()


}

// clean up unused code and commit to branch
// this version works for one chr