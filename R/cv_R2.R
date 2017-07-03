#' Cross-Validated R-Squared
#'
#' Returns a cross-validated estimate of R-squared for how well a fitted model
#' object is expected to explain new data.
#'
#' \code{cv_r2} estimates a \emph{predictive} R-squared (the fraction of
#' uncertainty in a new set of outcomes that is explained by a model fit to a
#' prior set of outcomes) using repeated \eqn{k}-fold cross-validation. The
#' conventional R-squared for describing models fit using ordinary least squares
#' (OLS) regression is calculated as \eqn{1 - RSS/TSS}, where \eqn{RSS} is the
#' residual sum of squares and \eqn{TSS} is the total sum of squares. The
#' predicted R-squared replaces \eqn{RSS} with \eqn{PRESS} (prediction sum of
#' squares), which, for OLS models, can be calculated formulaically using
#' leverage values but is equivalent to the sum of squared prediction errors
#' under leave-one-out cross-validation (LOOCV). However, LOOCV is known to have
#' high variance, and \eqn{k}-fold cross-validation gives more accurate
#' estimates of test error more often. Therefore, the PRESS statistic here is
#' calculated using \eqn{k}-fold cross-validation rather than LOOCV.
#'
#' \code{cv_r2} extends the predictive R-squared to generalized linear models
#' (GLMs) by using predicted deviance and null deviance in place of PRESS and
#' TSS, respectively: \eqn{1 - predicted deviance / null deviance}. Predicted
#' deviance and null deviance are equivalent to PRESS and TSS, respectively,
#' for models fit with OLS or GLMs with a Gaussian error distribution and
#' identity link function. But unlike variance-based R-squared, the deviance-
#' based R-squared generalizes to several exponential family models. (See
#' \code{\link{r2d}} for more information and references.)
#'
#' To obtain a cross-validated predictive R-squared, first fit a model as you
#' normally would using \code{\link[stats]{lm}}, \code{\link[stats]{glm}},
#' \code{\link[MASS]{glm.nb}}, or \code{\link[pscl]{zeroinfl}}, and then pass
#' the model object to \code{cv_r2}.
#'
#' @return a "\code{cv_R2}" object consisting of a list with the following
#' elements:
#' \describe{
#'   \item{cv_R2}{the mean R-squared of all folds of all repetitions of the
#'     cross-validation procedure}
#'   \item{95\% CI}{a bootstrapped 95\% confidence interval for \code{cv_R2}
#'     obtained using 1000 replicates and the adjusted bootstrap percentile
#'     method}
#'   \item{R2}{a named vector corresponding to the individual R-squared values
#'     obtained for each fold of each repetition of the cross-validation;
#'     provided in case the user wishes to manually calculate alternative
#'     summary statistics, e.g., the median R-squared or a different confidence
#'     interval}}
#'
#' @param object A model object for which a cross-validated R-squared is
#' desired.
#'
#' @param n_cores Integer value indicating the number of workers to run in
#' parallel during cross-validation and boot-strapping. By default, this will
#' be set to 2. To determine the theoretical maximum number of cores you have
#' available, see \code{\link[parallel]{detectCores}}, but note that the actual
#' number of cores available may be less. See
#' \code{\link[parallel]{parallel-package}} for more information.
#'
#' @param n_folds Integer indicating the number of folds to use for
#' cross-validation.
#'
#' @param n_repeats Integer indicating the number of times cross-validation
#' should be repeated.
#'
#' @param seed Integer used to seed the random number generator.
#'
#' @export
cv_r2 <- function(object, n_cores = 2, n_folds = 10, n_repeats = 10, seed = 42){
  y <- as_vector(object$model[1])
  if(is.null(y)) y <- as_vector(object$data[all.vars(object$terms)[1]])
  if(is.null(y)) stop(paste("Model data not found in model object.",
                            "Try refitting model with `model = TRUE` arg."))
  set.seed(seed)
  fold_ids <- caret::createMultiFolds(y, k = n_folds, times = n_repeats)
  cl <- parallel::makeCluster(n_cores)
  r2 <- parallel::parSapply(cl, fold_ids, function(i, object){
    form <- object$terms; fam <- object$family; link <- object$family$link
    if(class(object)[1] == "zeroinfl"){
      form <- object$terms$full
      link <- object$link
    }
    if(!is.null(object$model)){
      train_data <- object$model[i,]
      test_data <- object$model[-i,]
    } else {
      train_data <- object$data[i,]
      test_data <- object$data[-i,]
    }
    wts <- object$prior.weights[i]
    if(is.null(wts)){
      train_data$wts <- rep(1, nrow(train_data))
    } else {
      train_data$wts <- wts
    }
    dist <- object$dist; ctrl <- object$control; method <- object$method
    ctrst <- object$contrasts; theta <- object$theta
    fit <- switch(
      class(object)[1],
      lm = stats::lm(formula = form, data = train_data, weights = wts,
                     contrasts = ctrst),
      glm = stats::glm(formula = form, family = fam, data = train_data,
                       weights = wts, control = ctrl, method = method,
                       contrasts = ctrst),
      negbin = glm_nb(formula = form, data = train_data, weights = wts,
                      control = ctrl, method = method, contrasts = ctrst,
                      init.theta = theta, link = "log"),
      zeroinfl = pscl::zeroinfl(formula = form, data = train_data,
                                weights = wts, dist = dist, link = link,
                                control = ctrl)
    )
    stats <- try(predict_metrics(fit, test_data), silent = TRUE)
    if(class(stats) == "prediction_metrics") stats$R_squared else NA
  }, object = object)
  boot_r2 <- boot::boot(r2, function(x, i) mean(x[i], na.rm = TRUE), 1000,
                        parallel = "snow", cl = cl)
  parallel::stopCluster(cl)
  boot_CI <- boot::boot.ci(boot_r2, type = "bca")
  structure(list(cv_R2 = mean(r2, na.rm = TRUE), `95% CI` = boot_CI$bca[4:5],
            R2 = r2), class = "cv_R2")
}

