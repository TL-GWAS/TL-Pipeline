include { IIDGenotypes; GeneticConfounders; LOCOGenotypes; LOCOConfounders } from '../subworkflows/confounders.nf'
include { ExtractTraits } from '../subworkflows/extract_traits.nf'

workflow PCA {
    // Define Parameters
    ukb_encoding_file = params.UKB_ENCODING_FILE
    ukb_config = Channel.value(file("$params.UKB_CONFIG", checkIfExists: true))
    ukb_withdrawal_list = Channel.value(file("$params.UKB_WITHDRAWAL_LIST", checkIfExists: true))
    traits_dataset = Channel.value(file("$params.TRAITS_DATASET", checkIfExists: true))

    qc_file = Channel.value(file("$params.QC_FILE", checkIfExists: true))
    flashpca_excl_reg = Channel.value(file("$params.FLASHPCA_EXCLUSION_REGIONS", checkIfExists: true))
    ld_blocks = Channel.value(file("$params.LD_BLOCKS", checkIfExists: true))
    bed_files = Channel.fromFilePairs("$params.BED_FILES", size: 3, checkIfExists: true){ file -> file.baseName }
    loco_bed_files = Channel.fromFilePairs("$params.BED_FILES", size: 3, checkIfExists: true)

    // Extract Traits
    ExtractTraits(
        traits_dataset,
        ukb_config,
        ukb_withdrawal_list,
        ukb_encoding_file,
    )
    
    if (params.STUDY_DESIGN == "GWAS") {
        LOCOGenotypes(
            loco_bed_files, 
            ExtractTraits.out,
            qc_file,
            ld_blocks
        )
        genotypes = LOCOGenotypes.out

        // Genetic confounders
        LOCOConfounders(LOCOGenotypes.out)

        // Merge bed files with confounder output on chr_id
        confounders = LOCOConfounders.out.cross(loco_bed_files).map{[it[0][0], it[0][1], it[1][1][0], it[1][1][1], it[1][1][2]]}
    }
    else{
        // IID Genotypes
        IIDGenotypes(
            flashpca_excl_reg,
            ld_blocks,
            bed_files,
            qc_file,
            ExtractTraits.out,
        )
        genotypes = IIDGenotypes.out

        // Genetic confounders
        confounders = GeneticConfounders(IIDGenotypes.out)
    }  
    emit:
        traits = ExtractTraits.out
        iid_genotypes = genotypes
        confounders = confounders
}