
include { longest_prefix } from '../modules/utils.nf'

process GeneratePermutationTestsData {
    container "olivierlabayle/tl-core:0.8"
    publishDir "${params.OUTDIR}/permutation_tests", mode: 'symlink', pattern: '*.arrow'
    publishDir "${params.OUTDIR}/permutation_tests/estimands", mode: 'symlink', pattern: '*.jls'
    label "bigmem"
    
    input:
        path dataset
        path results

    output:
        path "permutation_dataset.arrow", emit: dataset
        path "*.jls", emit: estimands

    script:
        limit = params.MAX_PERMUTATION_TESTS == "" ? "" : "--limit=${params.MAX_PERMUTATION_TESTS}"
        """
        TEMPD=\$(mktemp -d)
        JULIA_DEPOT_PATH=\$TEMPD:/opt julia --project=/TargeneCore.jl --startup-file=no --sysimage=/TargeneCore.jl/TargeneCoreSysimage.so /TargeneCore.jl/bin/generate_tl_inputs.jl \
        --positivity-constraint=${params.POSITIVITY_CONSTRAINT} \
        --verbosity=${params.VERBOSITY} \
        --batch-size=${params.BATCH_SIZE} \
        --out-prefix="." \
        permutation-tests ${dataset} ${results} \
        ${limit} \
        --pval-threshold=${params.PVAL_THRESHOLD} \
        --estimator-key=${params.ESTIMATOR_KEY} \
        --orders=${params.PERMUTATION_ORDERS} \
        --rng=${params.RNG} \
        --max-attempts=${params.MAX_PERMUTATION_ATTEMPTS} \
        """
}

process GenerateRandomVariantsTestsData {
    container "olivierlabayle/tl-core:0.8"
    publishDir "${params.OUTDIR}", mode: 'symlink'
    label "bigmem"
    
    input:
        path variants_to_randomize
        path bgenfiles
        path results

    output:
        path "random_variants_estimands.jls"

    script:
        bgen_prefix = longest_prefix(bgenfiles)
        """
        TEMPD=\$(mktemp -d)
        JULIA_DEPOT_PATH=\$TEMPD:/opt julia --project=/TargeneCore.jl --startup-file=no --sysimage=/TargeneCore.jl/TargeneCoreSysimage.so /TargeneCore.jl/bin/generate_random_variant_estimands.jl \
        ${variants_to_randomize} ${results} ${bgen_prefix} \
        --out=random_variants_estimands.jls \
        --p=${params.N_RANDOM_VARIANTS} \
        --reltol=${params.MAF_MATCHING_RELTOL} \
        --estimator-key=${params.ESTIMATOR_KEY} \
        --pval-threshold=${params.PVAL_THRESHOLD} \
        --rng=${params.RNG} \
        --verbosity=${params.VERBOSITY}
        """
}