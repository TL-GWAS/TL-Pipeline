# Configuration

## Main Outputs

The workflow will produce the following main outputs in the output directory (`OUTDIR`, defaults to `results`):

- An HDF5 file containing estimation results (`HDF5_OUTPUT`, default: results.hdf5)
- An optional JSON file containing estimation results (`JSON_OUTPUT`)
- An optional `svp.hdf5` containing Sieve Variance Plateau corrected variance estimates.
- A Quantile-Quantile summary plot: QQ.png

## Arguments

- **`ESTIMATOR_FILE`**: Estimator name or Julia file containing the description of the Targeted Estimators to use. To be consistent it should match the argument provided to the previous TarGene run.
- **`BGEN_FILES`**: Path to imputed BGEN files from which the variants in `VARIANTS_LIST` will be extracted.
- **`BED_FILES`**: Path expression to PLINK BED files.
- **`TRAITS_DATASET`**: Path to a traits dataset. If you are running this for a non-UKBB cohort, your sample IDs must be specified in the first column of this CSV file, with the column name `SAMPLE_ID`.

## Main Options

- **`STUDY_DESIGN` (default: CUSTOM)**: How genetic variants and effect sizes are specified. See [the study design section](@ref "Study Designs").
- **`COHORT` (default: "UKBB")**: Current default for this is UKBB. If set to a value other than UKBB, this will not run UKBB-specific trait extraction.
- **`POSITIVITY_CONSTRAINT` (default: 0.01)**: When the list of estimands is generated or validated. Treatment variables' rarest configuration should have at least that frequency. For example if the treatment variables are two variants with minor allele A and T respectively. The rarest configuration will be (AA, TT) and should have a frequency of at least `POSITIVITY_CONSTRAINT`.
- **`NB_PCS` (default: 6)**: The number of PCA components to extract.
- **`SVP` (default: false)**: Whether Sieve Variance Plateau correction should be performed.
- **`FLASHPCA_EXCLUSION_REGIONS` (default: assets/exclusion_regions_hg19.txt)**: A path to the flashpca special exclusion regions.
- **`MAF_THRESHOLD` (default: 0.01)**: Only variants with that minor allele frequency are considered
- **`LD_BLOCKS`**: A path to pre-identified linkage disequilibrium blocks to be removed from the BED files. It is good practice to specify `LD_BLOCKS`, as it will remove SNPs correlated with your variants-of-interest before running PCA.

### If `COHORT=UKBB` (default)

- **`UKB_CONFIG` (default: assets/ukbconfig.yaml)**: YAML configuration file describing which traits should be extracted and how the population should be subsetted.
- **`UKB_ENCODING_FILE`**: If the `TRAITS_DATASET` is encrypted, an encoding file must be provided.
- **`UKB_WITHDRAWAL_LIST`**: List of participants withdrawn from the study.
- **`QC_FILE`**: Genotyping quality control file from the UK-Biobank study.

### If `STUDY_DESIGN=CUSTOM` (default)

- **`ESTIMANDS_FILE`**: YAML configuration file describing the effect sizes of interest, see the [custom study design section](@ref "`CUSTOM`")

### If `STUDY_DESIGN=ALLELE_INDEPENDENT`

- **`ESTIMANDS_FILE`**: YAML configuration file describing the effect sizes of interest, see the [allele independent section](@ref "`ALLELE_INDEPENDENT`").

### If `STUDY_DESIGN=FROM_ACTORS`

- **`BQTLS`**: A CSV file containing binding quantitative trait loci (bQTLs). If multiple transcription factors (TFs) are included in a single run, you must include a column called `TF`, which specifies the TF associated with each bQTL.
- **`TRANS_ACTORS`**: A prefix to CSV files containing quantitative trait loci potentially interacting with the previous bqtls. If multiple transcription factors (TFs) are included in a single run, you must include a column called `TF`, which specifies the TF associated with each transactor.
- **`EXTRA_CONFOUNDERS`, default: nothing**: Path to additional confounders file, one per line, no header.
- **`EXTRA_COVARIATES`, default: nothing**: Path to additional covariates file, one per line, no header.
- **`ENVIRONMENTALS`, default: nothing**: Path to additional environmental treatments file, one per line, no header.
- **`ORDERS`, default: "1,2"**: Comma separated list describing the order of desired interaction estimands, 1 for the ATE (no interaction), 2 for pairwise interactions etc... e.g. "1,2"

### If `SVP=true`

- **`GRM_NSPLITS`, default: 100**: To fasten GRM computation, it is typically split in batches.
- **`NB_SVP_ESTIMATORS`, default: 100**: Number of sieve variance estimates per curve. Setting this value to 0 results in skipping sieve variance correction.
- **`MAX_SVP_THRESHOLD`, default: 0.9**: Variance estimates are computed for tau ranging from 0 to MAX_SVP_THRESHOLD
- **`PVAL_THRESHOLD`, default: 0.05**: Only results with a p-value below this threshold are considered for Sieve Plateau Variance correction.
- **`ESTIMATOR_KEY`, default: TMLE**: The p-value for `PVAL_THRESHOLD` is computed using the result from this estimator.

## Secondary Options

- **`RNG`, default: 123**: General random seed used to induce permutation.
- **`VERBOSITY`, default: 0**: Verbosity level of the the Workflow's processes.
- **`BATCH_SIZE`, default: 400**: The set of estimands to be estimated is batched and the Targeted Learning processes will run in parallel across batches.
- **`TL_SAVE_EVERY`, default: 50**: During the estimation process, results are appended to the file in chunks to free memory.
- **`KEEP_IC`, default: SVP**: To save the Influence Curves for each estimate. Depending on the size of your dataset, this can result in massive disk usage.
