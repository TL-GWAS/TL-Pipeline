#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

def longest_prefix(files){
    // Only one file, strangely it is not passed as a list
    if (files instanceof Collection == false) {
        return files.getName()
    }
    // More than one file
    def index = 0
    while(true){
        def current_prefix = files[0].getName()[0..index]
        for (file in files){
            if(file.getName()[0..index] != current_prefix){
                return current_prefix[0..-2]
            }
        }
        index++
    }
}

def get_length(files) {
    def len = 0
    for (file in files){
        len++
    }
    return len
}

// NF parameter for PCA step
params.NB_PCS = 20
params.OUTDIR = "$launchDir/test/output"
params.BED_FILES = "/exports/igmm/eddie/UK-BioBank-53116/genotypes/ukb_*.{bed,bim,fam}"
params.QC_FILE = "/exports/igmm/eddie/UK-BioBank-53116/imputed/ukb_snp_qc.txt"
params.LD_BLOCKS = "NO_LD_BLOCKS"
params.TRAITS_CONFIG = "/exports/igmm/eddie/khamseh-lab/jslaughter/targene_dev/targene-pipeline/loco_gwas/test/data/ukbconfig_small.yaml"
params.WITHDRAWAL_LIST = 'NO_WITHDRAWAL_LIST'
params.COHORT = "UKBB"

// imports from main.nf

//  PCA processes from TarGene
include { FlashPCA; AdaptFlashPCA } from "$launchDir/modules/confounders.nf"
//QC
include { SampleQCFilter; filterBED} from "$launchDir/modules/genotypes.nf"

// Create channel that stores genotype info
input_channel_pairs = Channel.fromFilePairs(params.BED_FILES, size: 3, checkIfExists: true)
// "/exports/igmm/eddie/UK-BioBank-53116/genotypes/ukb_*.{bed,bim,fam}"
// $launchDir/test/data/phased_bed/*.{bed,bim,fam}
// process MergeBEDS {
//     label 'bigmem'
//     container "olivierlabayle/tl-core:0.6"
//     publishDir "$params.OUTDIR/genotypes", mode: 'symlink'

//     input:
//         path files
//     output:
//         path "ukbb_merged*"
//     script:
//         prefix = longest_prefix(files)
//         """
//         TEMPD=\$(mktemp -d)
//         JULIA_DEPOT_PATH=\$TEMPD:/opt julia --project=/TargeneCore.jl --startup-file=no /TargeneCore.jl/bin/prepare_confounders.jl --input ${prefix} --output ukbb_merged merge
//         """

// }

process LocoMergeBEDS {
    label 'bigmem'
    container "olivierlabayle/tl-core:0.6"

    input:
        tuple val(chr), path(files)
    output:
        tuple path("${chr}_excluded*"), val(chr)

    script:
        prefix = longest_prefix(files)
        len = get_length(files)
        // println('Chr: ' + chr + ', Prefix: ' + prefix + ", Length: " + len)

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
        // Maybe have a filter bed step here before the merge???
        merged_bed_ch = LocoMergeBEDS(bed_files)
        qc_filtered = SampleQCFilter(merged_bed_ch)
        pcs_txt = FlashPCA(qc_filtered)
        pcs_csv = AdaptFlashPCA(pcs_txt)

    emit:
        pcs_csv
        
}

workflow extractTraits {
    traits_config = Channel.value(file("$params.TRAITS_CONFIG"))
    withdrawal_list = Channel.value(file("$params.WITHDRAWAL_LIST"))
    if (params.DECRYPTED_DATASET == "NO_FILE") {
        encrypted_dataset = Channel.value(file("$params.ENCRYPTED_DATASET"))
        encoding_file = Channel.value(file("$params.ENCODING_FILE"))
        UKBFieldsList(traits_config)
        decrypted_dataset = UKBConv(UKBFieldsList.out, encrypted_dataset, encoding_file)
    }
    else {
        decrypted_dataset = Channel.value(file("$params.DECRYPTED_DATASET"))
    }

    if (params.COHORT == "UKBB") {
        extracted_traits = TraitsFromUKB(decrypted_dataset, traits_config, withdrawal_list)
    } 
    else {
        extracted_traits = decrypted_dataset 
    }

    emit:
        extracted_traits
}

// Place in subworkflow 
// workflow LOCO_Config {
    
// }
workflow {
    // Test to extract merged .bed
    // genotype_files = Channel.fromPath("/exports/igmm/eddie/UK-BioBank-53116/genotypes/ukb_*").collect()
    // genotype_files.view()
    // MergeBEDS(genotype_files)
    
    // DO NOT DELETE
    all_files = input_channel_pairs.map { it[1] }.collect().toList()

    // input_channel_pairs
    //     .combine(all_files)
    //     .map{ [it[0], it[2].findAll{f->!f.normalize().toString().contains(it[0])}] }
    //     .set{exclusion_set}

    input_channel_pairs
        .combine(all_files)
        .map{ [it[0], it[2].findAll{f-> ! it[1].collect{x -> x.normalize().toString()}.contains(f.normalize().toString())}] }
        .set{exclusion_set}
    RunPCALoco(exclusion_set)
    extractTraits()
}