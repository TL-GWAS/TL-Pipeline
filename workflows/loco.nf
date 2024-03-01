#!/usr/bin/env nextflow
nextflow.enable.dsl = 2
include { ExtractTraits } from "../subworkflows/extract_traits.nf"
include { FlashPCA; AdaptFlashPCA; SampleQCFilter; filterBED; LocoMergeBEDS} from "../modules/confounders.nf"

// Create channel that stores genotype info
input_channel_pairs = Channel.fromFilePairs(params.BED_FILES, size: 3, checkIfExists: true)

workflow LOCO_Input_Channels {
    take:
        input_channel_pairs
        traits
    main:
        filter_prep = input_channel_pairs.map {it[1]}
        filtered = filterBED(filter_prep, params.QC_FILE, params.LD_BLOCKS, traits)
        filtered = filtered.collect().toList()

        input_channel_pairs
            .combine(filtered)
            .map {[it[0], it[2].findAll { filePath -> 
                def filePathString = filePath.normalize().toString() 
                ! filePathString.contains(it[0].normalize().toString()+".") }] }
            .set {qc_exclusion_set}
    
    emit:
        qc_exclusion_set
}

workflow RunPCALoco {
    take:
        bed_files
    main:
        merged_bed_ch = LocoMergeBEDS(bed_files)
        qc_filtered = SampleQCFilter(merged_bed_ch)
        pcs_txt = FlashPCA(qc_filtered)
        pcs_csv = AdaptFlashPCA(pcs_txt)

    emit:
        pcs_csv
        
}

workflow {
    ExtractTraits(params.TRAITS_DATASET, params.TRAITS_CONFIG, params.UKB_WITHDRAWAL_LIST, params.UKB_ENCODING_FILE)
    LOCO_Input_Channels(input_channel_pairs, ExtractTraits.out)
    RunPCALoco(LOCO_Input_Channels.out)
    
}