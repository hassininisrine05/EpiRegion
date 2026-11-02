#' Train a penalised-logistic episignature classifier
#'
#' @param X Sample x feature matrix of methylation proportions.
#' @param y Two-level factor (control, case).
#' @param alpha Elastic-net mixing (0 = ridge, 1 = lasso).
#' @param nfolds Cross-validation folds for lambda selection.
#' @return A `cv.glmnet` object with feature means stored for imputation.
#' @export
fit_episignature <- function(X, y, alpha = 0, nfolds = 5) {
  X <- as.matrix(X)
  mu <- colMeans(X, na.rm = TRUE)
  X[is.na(X)] <- mu[col(X)][is.na(X)]
  fit <- glmnet::cv.glmnet(X, factor(y), family = "binomial",
                           alpha = alpha, nfolds = nfolds)
  fit$feature_means <- mu
  fit
}

#' Predict case probability
#' @param fit Output of [fit_episignature()].
#' @param Xnew Sample x feature matrix.
#' @export
predict_episignature <- function(fit, Xnew) {
  Xnew <- as.matrix(Xnew)
  Xnew[is.na(Xnew)] <- fit$feature_means[col(Xnew)][is.na(Xnew)]
  as.numeric(stats::predict(fit, Xnew, s = "lambda.min", type = "response"))
}
