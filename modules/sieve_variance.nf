process SieveVarianceEstimation {
    container "olivierlabayle/tl-core:sample_filtering"
    publishDir "$params.OUTDIR/hdf5files/sieve_variances", mode: 'symlink'

    input:
        tuple val(rsprefix), file(estimate_file)
        path GRM_ids
        path GRM_matrix

    output:
        tuple val(rsprefix), file("${rsprefix}_sieve_variance.hdf5")
    
    script:
        "julia --project=/TargeneCore.jl --startup-file=no /TargeneCore.jl/bin/sieve_variance.jl $rsprefix GRM ${rsprefix}_sieve_variance.hdf5 --nb-estimators=$params.NB_VAR_ESTIMATORS --max-tau=$params.MAX_TAU --pval=$params.PVAL_SIEVE"
}

process NoSieveEstimation {
    input:
        tuple val(rsprefix), file(estimate_file)

    output:
        tuple val(rsprefix), val("NOFILE")
        
    script:
        "exit 0"
}