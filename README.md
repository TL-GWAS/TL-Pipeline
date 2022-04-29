# TL-Pipeline

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://olivierlabayle.github.io/UKBBEpistasisPipeline.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://olivierlabayle.github.io/UKBBEpistasisPipeline.jl/dev)
[![Build Status](https://github.com/olivierlabayle/UKBBEpistasisPipeline.jl/workflows/CI/badge.svg)](https://github.com/olivierlabayle/UKBBEpistasisPipeline.jl/actions)
[![Coverage](https://codecov.io/gh/olivierlabayle/UKBBEpistasisPipeline.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/olivierlabayle/UKBBEpistasisPipeline.jl)

Here we provide the main workflow for estimating genetic variants interactions responsible of traits in the UK Biobank. For that purpose, we rely on [Singularity](https://sylabs.io/guides/3.0/user-guide/quick_start.html) and [NextFlow](https://www.nextflow.io/).

## Description of the Workflow

The workflow is divided into 4 steps, where the first 3 steps generate the inputs for the final TMLE process: 
1. The generation of the queries files. Each query file specifies a set of potentially interacting variants together with the chromosome they are located in and for each variant which allele is the `control` and which allele is the `treatment` value.
2. The generation of the phenotypes file. From the UK-Biobank main dataset, generates a set of traits of interest for the study.
3. The generation of the confounders file. To adjust for confounding effects in the final estimation step.
4. The TMLE estimation step, requiring all 3 previous steps.

We now describe the pipeline arguments for each sub-workflow:

### Queries Generation

There are currently two possible procedures which are specified by the `QUERIES_MODE` argument.
- `QUERIES_MODE` = "ASBxTransActors". It is assumed that an initial set of variants has been pre-identified from a previous allele-specific binding (ASB) study as output by the [ball-nf](https://git.ecdf.ed.ac.uk/oalmelid/baal-nf) pipeline. It is also assumed that another set of potential trans-actors if given in a csv file. The cross product of those two sets of variants is taken to generate the queries. For more information on the procedure and expected format of those files see: `julia bin/generate_queries.jl --help`. Two arguments thus have to be passed:
    - `ASB_FILES`: Path to all files output by the ball-nf pipeline.
    - `TRANS_ACTORS_FILE`: Path to the trans-actors files.
    - Additionaly, the `THRESHOLD` argument can be overided to specify the hard-calling threshold to convert probabilities to genotypes values.

- `QUERIES_MODE` = "given". This is useful if only a few query files have to be generated and can be written "by hand". Then the following argument has to be provided:
    - `QUERY_FILES`: Path to the query files.

### Covariates Generation

To account for potential confounding effect due to population stratification, we extract principal components from the genetic data using [flashpca](https://github.com/gabraham/flashpca). We follow the recommended procedure for this tool which implies some preprocessing and filtering. The following arguments are compulsory:
- `QC_FILE`: A path to the UK-biobank SNP quaility control `ukb_snp_qc.txt` file.
- `LD_BLOCKS`: A path to pre-identified linkage desequlibrium blocks around the variants that will be queried for causal effect estimation. Those will be removed from the data.
- `FLASHPCA_EXCLUSION_REGIONS`: A path to the flashpca special exclusion regions which is provided in their repository.
- `NB_PCS`: The number of PCA components to extract.
- `UKBB_BED_FILES`: PCA components are built from the UK-biobankk bed files which must be provided.

### Phenotypes Generation

For now, we rely on the phenotypes that have been considered and generated by the [GeneAtlas](http://geneatlas.roslin.ed.ac.uk/) study. Because each study is independent, the sample ids in each study are unique. This means a bridge file has to be provided to link the data from the GeneAtlas and your study. The arguments are:
- `BINARY_PHENOTYPES`: A path to the binary phenotypes file from the GeneAtlas.
- `CONTINUOUS_PHENOTYPES`: A path to the continuous phenotypes file from the GeneAtlas.
- `GENEATLAS_BRIDGE`: A path to the bridge file linking your study to the GeneAtlas study.
- `WITHDRAWAL_LIST`: A path to the withdrawal sample list to exclude removed participants from the study.
- `PHENOTYPES_LIST`: A file containing a list of the phenotypes from the GeneAtlas to use, one line for each phenotype

### Genetic Relationship Matrix

This is used solely by the Sieve Variance correction step and uses the GCTA software under the hood.

- `GRM_NSPLITS`: The number of sub grm parts to be computed since the full GRM requires more than 1TB of memory.
### TMLE

Almost the last step of the pipeline: targeted estimation. This step uses "SuperLearning" (Stacking) which means a configuration for this learning step has to be provided:

- `ESTIMATORFILE`: This configuration will be used to build super-learners estimators
- `CROSSVAL`: A boolean to indicate weither an additional cross validation procedure is run to evaluate the models in each super learning step. Unfortunately this is currently run aside of TMLE and thus quite expensive. Only use it if you really need those figures.
- `PHENOTYPES_BATCH_SIZE`: Estimation is parallelized over queries. If the number of queries is low it can be advantageous to parallelize over phenotypes too.


### Sieve Variance Estimation

Arguments for the sieve variance correction step:

- `NB_VAR_ESTIMATORS`: Number of estimators to compute, the interval [0, MAX_TAU], will be split.
- `MAX_TAU`: Maximum distance up to which individuals will be included in the estimation. Previous analysis showed that above 0.8, some side effects seem to occur. The theoretical maximum is 2.
- `PVAL_SIEVE`: Only traits for which the IID pvalue is lower than this threshold will be considered 
for sieve variance correction. This is because in theory the Sieve Variance curve is supposed to be monotically increasing.

## Running the Workflow

If you are part of the University of Edinburgh, you can use the Eddie cluster for which the `eddie` profile in this repo can be used directly. For more detailed information, refer to [Eddie specific Nextflow documentation](https://www.wiki.ed.ac.uk/display/ResearchServices/Bioinformatics). Again, if you have access rights, you can see two example worflows at:
-   [The VDR project](https://git.ecdf.ed.ac.uk/tfomics/vdr)
-   [The FTO project](https://git.ecdf.ed.ac.uk/tfomics/uk-biobank/tmle-rs1421085)

The idea is that you can simply run:

```bash
nextflow run olivierlabayle/UKBBEpistasisPipeline.jl -profile eddie
```

If you are not part of the University and cannot access Eddie, you will have to adapt the configuration file to your specific platflorm, fill an issue to get in contact.
