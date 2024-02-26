# Specifying a Targeted Estimator

TarGene is a flexible procedure that does not impose any constraint on the functional form of the relationship between genetic variants, environmental variables and traits. In practice, we rely on [MLJ](https://alan-turing-institute.github.io/MLJ.jl/dev/) to provide machine learning algorithms. In population genetics studies, there are two learning algorithms we need to specify:

- `E[Y|T, W, C]`: The mean outcome given the treatment, confounders and extra covariates. It is commonly denoted by `Q` in the Targeted Learning literature. In reality, we will need one specification for continuous outcomes and one specification for binary outcomes.
- `p(T|W)`: The propensity score, which enables the targeting step of the estimation procedure. It is commonly denoted by `G` in the Targeted Learning literature.

There are two main ways to define targeted estimators, from a predefined configuration or from a custom file.

## Predefined estimators

A set of predefined estimators is readily available and can be accessed by using the configuration's name as the `ESTIMATOR_FILE` parameter. For instance, we provide the following:

- G-superlearning-Q-glm
- G-superlearning-Q-glmnet
- glm-with-interactions-for-Q
- glm
- glmnet-with-interactions-for-Q
- glmnet
- superlearning-with-interactions-for-Q
- superlearning
- tuned-xgboost

all using TMLE as the meta statistical inference method. While these should cover most use cases, it may be useful to define a custom estimation strategy. This can be achieved by writing a small [Julia](https://julialang.org/) file described below.

## Custom estimators from a file

Writing a [Julia](https://julialang.org/) estimators file is the most flexible way to define an estimation strategy for your study.

This file should simply define a [NamedTuple](https://docs.julialang.org/en/v1/base/base/#Core.NamedTuple) called `ESTIMATORS` and listing estimators to be used for inference. For example, the following:

```julia
default_models = TMLE.default_models(
  # For the estimation of E[Y|W, T]: continuous outcome
  Q_continuous = LinearRegressor(),
  # For the estimation of E[Y|W, T]: binary target
  Q_binary = LogisticClassifier(lambda=0.),
  # For the estimation of p(T| W)
  G = LogisticClassifier(lambda=0.)
)

ESTIMATORS = (
  TMLE_weighted   = TMLEE(models=default_models, weighted=true),
  TMLE_unweighted = TMLEE(models=default_models, weighted=false),
  OSE             = OSE(models=default_models)
)
```

defines three estimators:

1. A weighted-fluctuation Targeted Maximum Likelihood Estimator: `TMLE`
2. An unweighted-fluctuation Targeted Maximum Likelihood Estimator: `TMLE`
3. A One-Step Estimator: `OSE`

All estimators will learn the nuisance functions `Q` and `G` with the provided `models` `NamedTuple`:

- `Q_continuous`: A MLJ model used for the estimation of `E[Y|T, W, C]` when the outcome `Y` is continuous.
- `Q_binary`: A MLJ model used for the estimation of `E[Y|T, W, C]` when the outcome `Y` is binary.
- `G`: A MLJ model used for the estimation of `p(T|W)`.

For the list of available models and resampling strategies, checkout the [TargetedEstimation documentation](https://targene.github.io/TargetedEstimation.jl/stable/models/).

For full details, on available estimators and how to specify them, visit the [TMLE.jl documentation](https://targene.github.io/TMLE.jl/stable/).
