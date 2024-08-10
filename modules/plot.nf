process GenerateSummaryPlots {
    publishDir "${params.OUTDIR}", mode: 'symlink'
    label "bigmem"
    label 'targenecore_image'

    input:
        path results_file

    output:
        path "*.png"

    script:
        """
        TEMPD=\$(mktemp -d)
        JULIA_DEPOT_PATH=\$TEMPD:/opt julia --project=/TargeneCore.jl --startup-file=no --sysimage=/TargeneCore.jl/TargeneCoreSysimage.so /TargeneCore.jl/targenecore.jl \
        summary-plots ${results_file} \
        --outprefix="." \
        --verbosity=${params.VERBOSITY}
        """
}