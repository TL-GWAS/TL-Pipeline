process Summary {
    container "olivierlabayle/tl-core:sample_filtering"
    publishDir "$params.OUTDIR/summaries", mode: 'symlink'
    label "bigmem"

    input:
        tuple val(rsids), file(tmle), file(sieve)
    
    output:
        path "${rsids}_summary.csv"
    
    script:
        sieve_option = sieve.getName() != 'NO_FILE' ? "--sieve" : ''
        """
        julia --project=/TMLEEpistasis.jl --startup-file=no \
        /TMLEEpistasis.jl/bin/summarize.jl $rsids ${rsids}_summary.csv $sieve_option
        """
}