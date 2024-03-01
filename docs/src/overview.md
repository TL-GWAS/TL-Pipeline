# Overview

## Running the Workflows

Since TarGene uses [Nextflow](https://www.nextflow.io/), all workflows can be run in the same way from the command line:

```bash
nextflow run https://github.com/TARGENE/targene-pipeline/ -r TAG -entry WORKFLOW_NAME -profile P -resume
```

where:

- `TAG` is the latest TarGene version, e.g. `v0.9.0`
- `WORKFLOW_NAME` is any of the [TarGene workflows](@ref "Project Configuration")
- `P` is a [Nextflow profile](https://www.nextflow.io/docs/latest/config.html) describing your run environment (see [Environment Configuration](@ref)).

Additional Nextflow command line arguments can be found in their documentation, for example:

- `-resume`: Tells Nextflow to try to resume the pipeline if an error occurred during the execution (if you forgot to specify a parameter for instance)
- `-with-trace` and `-with-report` will generate additional report files.

## Workflows Configurations

There are mainly two parts to configuring a workflow run. The first part describes the execution and environment and the second part describes your actual project. When running the `nextflow run` command, Nextflow will look for a `nextflow.config` configuration file in your current directory. If you are new to Nextflow, we suggest you use this file to setup both the environment and project configurations. As your project grows you may want to benefit from splitting this file in more [modular components](https://www.nextflow.io/docs/latest/config.html#configuration).

### Environment Configuration

It is likely that you will run TarGene on a HPC platform, in particular the [Executors](https://www.nextflow.io/docs/latest/executor.html) and [Singularity](https://www.nextflow.io/docs/latest/container.html#singularity) configurations are required. Since Nextflow is so widespread, it is probable that such a configuration file is already available from your HPC administrators. Since this configuration only describes de computing platform and not your project, it is often described as a [Profile](https://www.nextflow.io/docs/latest/config.html#config-profiles). If your HPC uses the SGE executor, the `-profile eddie` may work with no, or minor adjustment (it can also serve as a template for other executors [see file](https://github.com/TARGENE/targene-pipeline/blob/main/conf/eddie.config)).

### Project Configuration

These are the configuration details associated with your project, this is usually done in a `nextflow.config` file living at the root of your project's directory. The configuration parameters are specific to each workflow and described in the following sections. There are currently two main workflows and two secondary workflows within TarGene.

#### Main Workflows

- [The TarGene Workflow](@ref) (`WORKFLOW_NAME: TARGENE`): It is the main workflow for the targeted estimation of genetic effects.
- [The Negative Control Workflows](@ref "Negative Control Overview"): These workflows enable the control of the false discovery rate by using the results obtained from a previous TarGene discovery run. There are currently two of them:
  - [The Permutation Test Workflow](@ref) (`WORKFLOW_NAME: PERMUTATION_TEST`): Performs permutation tests by independently shuffling the individuals in the columns of an aggregated dataset.
  - [The Randomized Variants Workflow](@ref) (`WORKFLOW_NAME: RANDOMIZATION_TEST`): When there are multiple genetic variants of interest, e.g. in an interaction study, one can replace one of the variant at random by another variant and the effect is expected to be 0 in average.

#### Secondary Workflows

- [The PCA Workflow](@ref) (`WORKFLOW_NAME: PCA`): This workflow computes principal components.
- [The Make Dataset Workflow](@ref) (`WORKFLOW_NAME: MAKE_DATASET`): This workflow generates an aggregated dataset from traits and genetic data.
