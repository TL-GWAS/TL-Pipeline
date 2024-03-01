process UKBFieldsList {
    container "olivierlabayle/ukbmain:0.4"

    input:
        path ukb_config

    output:
        path "fields_list.txt"

    script:
        ukb_config = ukb_config.getName() != 'NO_UKB_TRAIT_CONFIG' ? "--conf $ukb_config" : ''
        """
        TEMPD=\$(mktemp -d)
        JULIA_DEPOT_PATH=\$TEMPD:/opt julia --project=/UKBMain.jl --startup-file=no /UKBMain.jl/scripts/build_fields_list.jl $ukb_config
        """
}

process UKBConv {
    container "olivierlabayle/ukbmain:0.4"

    input:
        path fields_list
        path encrypted_dataset
        path encoding_file

    output:
        path "ukb_dataset.csv"

    script:
        "ukbconv $encrypted_dataset csv -i$fields_list -e$encoding_file -oukb_dataset"
}

process TraitsFromUKB {
    publishDir "$params.OUTDIR/traits", mode: 'symlink'
    container "olivierlabayle/ukbmain:0.4"
    label "bigmem"
    label "multithreaded"

    input:
        path dataset
        path ukbconfig
        path withdrawal_list
    
    output:
        path 'traits.csv'

    script:
        ukbconfig = ukbconfig.getName() != 'NO_UKB_TRAIT_CONFIG' ? "--conf $ukbconfig" : ''
        withdrawal_list = withdrawal_list.getName() != 'NO_WITHDRAWAL_LIST' ? "--withdrawal-list $withdrawal_list" : ''
        """
        TEMPD=\$(mktemp -d)
        JULIA_DEPOT_PATH=\$TEMPD:/opt julia --project=/UKBMain.jl --startup-file=no --threads $task.cpus /UKBMain.jl/scripts/process_main_dataset.jl \
        $dataset $ukbconfig $withdrawal_list
        """
}