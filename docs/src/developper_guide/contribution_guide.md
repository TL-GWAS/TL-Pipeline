# Contributing

Contributions, whether bug fixes, new features or documentation improvements are very welcome!

## Raise an issue

In order to discuss and track the evolution of the project, please first raise an issue on the [targene-pipeline](https://github.com/TARGENE/targene-pipeline/issues) repository. If a change is agreed upon, the discussion should identify the relevant repositories that are concerned by the change and open an issue on each of the repository. For instance, if one wishes to improve the extraction of traits from the UK-Biobank, the [UKBMain.jl](https://github.com/TARGENE/UKBMain.jl) would surely be impacted and a new release for that package necessary.

## Workflow

Following our previous UKBMain.jl example, there are two repositories that need to be updated, the current workflow is as follows:

1. Develop
    - UKBMain.jl
        - Create a new git branch for your change
        - Develop and test
        - Release an image for your branch by selecting it after clicking the [Run workflow button](https://github.com/TARGENE/UKBMain.jl/actions/workflows/Release.yml). If the tests pass, a new docker image will be generated and hosted on Docker hub with your branch's name (see [Note on Docker images](@ref)).
    - targene-pipeline
        - Create a new git branch for your change
        - For each Nextflow process using the UKBMain.jl's docker image, update to the branch's image name.
        - Develop further required changes and run/add the tests (see [Note on the pipeline's tests](@ref)).
2. Review: When everything is working, ask for a review
3. Release UKBMain.jl:
    - Merge your branch into main
    - Create a new Github release following semantic versioning, this will create a new docker image with your release name.
4. Release TarGene
    - For each Nextflow process using the UKBMain.jl's docker image, update to the released image name (as before).
    - Create a new Github release following semantic versioning

## Note on Docker images

Currently, all TarGene building blocks (executables) are provided as docker images. The following table provides a map linking each TarGene repository to the associated Docker image tags.

| Repository | Docker tag |
| --- | --- |
| [TargeneCore.jl](https://github.com/TARGENE/TargeneCore.jl) | [tl-core](https://hub.docker.com/r/olivierlabayle/tl-core/tags) |
| [UKBMain.jl](https://github.com/TARGENE/UKBMain.jl) | [ukbmain](https://hub.docker.com/r/olivierlabayle/ukbmain/tags) |
| [TargetedEstimation.jl](https://github.com/TARGENE/TargetedEstimation.jl) | [targeted-estimation](https://hub.docker.com/r/olivierlabayle/targeted-estimation/tags) |

## Note on the pipeline's tests

The pipeline is automatically tested for every push/pull-request made to the github repository. To run the tests locally, you will need Julia and the container engine of your choice to be installed (see below). Each test corresponds to a pipeline run and a file in the `test` directory is associated with it. Each test run can be launched locally as follows:

```bash
julia --project=test --startup-file=no test/TESTFILE -profile PROFILE -resume
```

where `TESTFILE` is the corresponding file and `PROFILE` is one of: 'local' (local with singularity engine), 'ci' (local with singularity engine), 'eddie' (SGE with singularity engine.).

The tests can also be run interactively via the Julia REPL.
