include { ExtractTraits } from '../subworkflows/extract_traits.nf'

workflow EXTRACT_TRAITS {
    ukb_encoding_file = params.UKB_ENCODING_FILE
    ukb_config = Channel.value(file("$params.UKB_CONFIG", checkIfExists: true))
    ukb_withdrawal_list = Channel.value(file("$params.UKB_WITHDRAWAL_LIST", checkIfExists: true))
    traits_dataset = Channel.value(file("$params.TRAITS_DATASET", checkIfExists: true))

    ExtractTraits(
        traits_dataset,
        ukb_config,
        ukb_withdrawal_list,
        ukb_encoding_file,
    )

    emit:
        ExtractTraits.out
}