include { FilterBED; ThinByLD; MergeBEDS; MergeBEDSLOCO; SampleQCFilter; SampleQCFilterLOCO; FlashPCA; FlashPCALOCO; AdaptFlashPCA; AdaptFlashPCALOCO } from '../modules/confounders.nf'

workflow LOCOGenotypes {
    take:
        loco_bed_files
        traits
        qc_file
        ld_blocks
    main:
        filter_prep = loco_bed_files.map {it[1]}
        filtered = FilterBED(filter_prep, qc_file, ld_blocks, traits)
        filtered = filtered.collect().toList()

        loco_bed_files
            .combine(filtered)
            .map {[it[0], it[2].findAll { filePath -> 
                def filePathString = filePath.normalize().toString() 
                ! filePathString.contains(it[0].normalize().toString()+".") }] }
            .set {qc_exclusion_set}
    
    emit:
        qc_exclusion_set
}

workflow LOCOConfounders{
    take:
        bed_files
    main:
        merged_bed_ch = MergeBEDSLOCO(bed_files)
        qc_filtered = SampleQCFilterLOCO(merged_bed_ch)
        pcs_txt = FlashPCALOCO(qc_filtered)
        pcs_csv = AdaptFlashPCALOCO(pcs_txt)

    emit:
        pcs_csv // this is a tuple [0] = chr, [1] = pcs
        
}

workflow IIDGenotypes{
    take:
        flashpca_excl_reg
        ld_blocks
        bed_files
        qc_file
        traits

    main:
        filtered_bedfiles = FilterBED(bed_files, qc_file, ld_blocks, traits)
        ld_pruned = ThinByLD(flashpca_excl_reg, filtered_bedfiles)
        MergeBEDS(ld_pruned.collect())
        SampleQCFilter(MergeBEDS.out.collect())

    emit:
        SampleQCFilter.out
}

workflow GeneticConfounders {
    take:
        iid_genotypes

    main:
        FlashPCA(iid_genotypes)
        AdaptFlashPCA(FlashPCA.out)

    emit:
        AdaptFlashPCA.out
}