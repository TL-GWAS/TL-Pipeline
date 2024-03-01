module TestUKBAlleleIndependent

using Test
using Arrow
using DataFrames

args = length(ARGS) > 0 ? ARGS : ["-profile", "local", "-resume"] 

@testset "Test MAKE_DATASET" begin
    cmd = `nextflow run main.nf -entry MAKE_DATASET -c test/configs/dataset.config $args`
    @info string("The following command will be run:\n", cmd)

    r = run(cmd)
    @test r.exitcode == 0

    dataset = Arrow.Table("results/dataset.arrow") |> DataFrame
    @test size(dataset, 1) == 500
    # PCs, SAMPLE_ID and the variant from the list
    @test all(n ∈ names(dataset) for n in ("SAMPLE_ID", "2:14983:G:A", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6",))
end

end

true