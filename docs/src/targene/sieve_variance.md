# Correcting for population relatedness

If the i.i.d. (independent and identically distributed) hypothesis is not satisfied, most of the traditional statistical inference theory falls apart. This is typically possible in population genetics where a study may contain related individuals. Here we leverage a non-parametric method called [Sieve Variance Plateau](https://biostats.bepress.com/ucbbiostat/paper322/) (SVP) estimation. The hypothesis is that the dependence between individuals is sufficiently small, so that our targeted estimator will still be asymptotically unbiased, but its variance will be under estimated.
In brief, the SVP estimator computes a variance estimate for a range of thresholds 𝜏, by considering individuals to be genetically independent if their genetic distance exceeds 𝜏. The genetic distance between a pair of individuals (𝑖, 𝑗) equals 1 − GRM𝑖,𝑗 , i.e., one minus their genetic relatedness value. As the distance threshold 𝜏 increases, fewer individuals are assumed to be genetically independent. For instance, the estimate corresponding to a distance of 𝜏 = 0 corresponds to the i.i.d. hypothesis, while a distance of 𝜏 = 1 incorporates pairs of individuals who are not genetically correlated. TarGene varies the threshold 𝜏 from 0 to 1 and fits a curve to the corresponding variance estimates. The maximum of this curve is the most conservative estimate of the variance estimator and constitutes our corrected variance estimator.

The following arguments can be changed to control the behaviour of the pipeline:

- `SVP` (default: false): Must be set to true to enable variance adjustement.
- `ESTIMATOR_KEY` (default: TMLE): The estimator from the `ESTIMATOR_FILE` that will be used for sieve variance plateau adjustment.
- `GRM_NSPLITS` (default: 100): This is a purely computational argument. The GRM is typically very large and splitting enables a good memory/parallelization tradeoff.
- `MAX_SVP_THRESHOLD` (default: 0.8): Controls the maximum genetic distance considered.
- `NB_SVP_ESTIMATORS` (default: 0): Controls the number of points in the interval [0, `MAX_SVP_THRESHOLD`]. If 0, the Sieve Variance Plateau method will not be applied.
- `PVAL_THRESHOLD` (default: 0.05): Only estimates with a p-value lower than `PVAL_THRESHOLD` will be considered for SVP correction. This is because SVP will only increase the variance of the estimator.
