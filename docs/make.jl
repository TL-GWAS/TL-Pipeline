using Documenter

# DocMeta.setdocmeta!(TarGene, :DocTestSetup, :(using TarGene); recursive=true)

makedocs(;
    # modules=[TarGene],
    authors="Olivier Labayle",
    repo="https://github.com/TARGENE/targene-pipeline/blob/{commit}{path}#{line}",
    sitename="TarGene",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://TARGENE.github.io/targene-pipeline",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Workflows" => [
            "overview.md",
            "The TarGene Workflow " => [
                joinpath("targene", "overview.md"),
                joinpath("targene", "data_sources.md"),
                joinpath("targene", "confounding_adjustment.md"),
                joinpath("targene", "study_designs.md"),
                joinpath("targene", "tmle.md"),
                joinpath("targene", "sieve_variance.md"),
                joinpath("targene", "miscellaneous.md"),
                joinpath("targene", "configuration.md"),
                joinpath("targene", "runtime_considerations.md")
            ],
            "The Simulation Workflows" => [
                joinpath("simulations", "overview.md"),
                joinpath("simulations", "null_simulation.md"),
            ],
            "Negative Control" => [
                joinpath("negative_control", "overview.md"),
                joinpath("negative_control", "permutation_tests.md"),
                joinpath("negative_control", "randomized_tests.md")
            ],
            joinpath("secondary_workflows", "pca.md"),
            joinpath("secondary_workflows", "make_dataset.md")
            ],
        "Developper Guide" => [
            joinpath("developer_guide", "project_organization.md"),
            joinpath("developer_guide", "contribution_guide.md")
            ],
        "Associated Softwares" => "associated_softwares.md",
        "Related Publications" => "publications.md"
    ],
)

deploydocs(;
    repo="github.com/TARGENE/targene-pipeline",
    devbranch="main",
    push_preview=true
)
