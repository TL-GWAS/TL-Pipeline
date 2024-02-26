# Negative Control Overview

After running a TarGene workflow, you may want to perform some quality control checks on your positive hits. This is done via negative control. By negative control we mean inducing a transformation on the dataset that alters the relationship between the variables of interest. Typically, between a genetic variant and a trait. The expectation is that any subsequent TarGene run, on this transformed dataset, will not deviate from the null hypothesis of no association.

At the moment, there are two main strategies for negative control in TarGene:

- [The Permutation Test Workflow](@ref)
- [The Randomized Variants Workflow](@ref)
