# Tweaking additional behaviour

Further Nextflow parameter affecting the behaviour of the pipeline but that does not fit in any previously described category is listed here:

- `CALL_THRESHOLD` (optional, default: 0.9): For putative causal variants (listed in the parameter files described in the [Study Designs](@ref) section). If an individual's allele's probability is greater than the threshold, then it is called, otherwise it is considered missing.
- `BATCH_SIZE` (optional, default: 400): The set of estimands to be estimated is batched and the TMLE processes will run in parallel across batches on your platform.
- `SVP` (optional, default: false): To activate Sieve Variance Plateau correction. If set to true, `KEEP_IC` will automatically be set to true as well.
- `KEEP_IC` (optional, default: false): For all estimands with a p-value below `PVAL_THRESHOLD`, the influence curve is saved.
- `NB_SVP_ESTIMATORS` (optional, default: 100): The number of estimators in Sieve Variance Plateau.
- `OUTDIR` (optional, default: "results"): Output directory
