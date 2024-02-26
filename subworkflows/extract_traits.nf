include { UKBConv; TraitsFromUKB; UKBFieldsList } from '../modules/extract_traits.nf'

workflow ExtractTraits {
    take:
        traits_dataset
        ukb_config
        ukb_withdrawal_list
        ukb_encoding_file
        
    main:
        if (params.COHORT == "UKBB") {
            if (ukb_encoding_file != "NO_UKB_ENCODING_FILE") {
                UKBFieldsList(ukb_config)
                decrypted_dataset = UKBConv(
                    UKBFieldsList.out, 
                    traits_dataset, 
                    Channel.value(file("${ukb_encoding_file}", checkIfExists: true))
                )
            }
            else {
                decrypted_dataset = traits_dataset
            }
            extracted_traits = TraitsFromUKB(decrypted_dataset, ukb_config, ukb_withdrawal_list)
        } 
        else {
            extracted_traits = traits_dataset 
        }

    emit:
        extracted_traits
}