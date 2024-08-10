# The Simulations Workflows

The **asymptotic** performance of semi-parametric estimators is theoretically optimal. However, in practice, asymptotic regimes are not necessarily achieved even for large sample sizes because some events are extremely rare. This is particularly prevalent in population genetics where some genetic variants and traits are found in less than ``1\%`` of individuals. In the absence of finite sample guarantees, simulation studies provide an effective way to validate statistical methods. An ideal simulation would both yield ground-truth values for the causal effects of interest and be representative of real data. In TarGene, we propose two types of simulations to analyse the performance of semi-parametric estimators. 

The first type of simulations, called the Null Simulation, aims at analysing the behaviour of the estimators when the Null-Hypothesis of "no effect" is true. The rationale behind this simulation is that most genetic variants are believed to have no effect, it is thus of particular importance that the type 1 error rate be controlled appropriately. 

The second type of simulations, called the Realistic Simulation, fits the generating process using data-adaptive flexible machine-learning generative models. As such, it retains important features of the true generating process, but also provides ground truth values.
