default_models = TMLE.default_models(
  # For the estimation of E[Y|W, T]: continuous outcome
  Q_continuous = LinearRegressor(),
  # For the estimation of E[Y|W, T]: binary target
  Q_binary = LogisticClassifier(lambda=0.),
  # For the estimation of p(T| W)
  G = LogisticClassifier(lambda=0.)
)

ESTIMATORS = (
  wTMLE_GLMNET = TMLEE(models=default_models, weighted=true),
  TMLE_GLMNET = TMLEE(models=default_models, weighted=false),
  OSE_GLMNET = OSE(models=default_models)
)