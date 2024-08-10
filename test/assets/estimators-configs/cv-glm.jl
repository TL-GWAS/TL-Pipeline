resampling = JointStratifiedCV(patterns=[r"^rs[0-9]+"], resampling=StratifiedCV(nfolds=3))

default_models = TMLE.default_models(
  # For the estimation of E[Y|W, T]: continuous outcome
  Q_continuous = LinearRegressor(),
  # For the estimation of E[Y|W, T]: binary target
  Q_binary = LogisticClassifier(lambda=0.),
  # For the estimation of p(T| W)
  G = LogisticClassifier(lambda=0.)
)

ESTIMATORS = (
  CV_wTMLE_GLMNET = TMLEE(models=default_models, weighted=true, resampling=resampling),
  CV_TMLE_GLMNET = TMLEE(models=default_models, weighted=false, resampling=resampling),
  CV_OSE_GLMNET = OSE(models=default_models, resampling=resampling)
)