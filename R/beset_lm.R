#' Best Subset Selection for Least-Squares Regression
#'
#' \code{beset_lm} performs best subset selection using repeated
#' cross-validation to find the optimal number of linear regression parameters.
#'
#' \code{beset_lm} performs best subset selection, using
#' \code{\link[leaps]{regsubsets}} to fit a separate least-squares regression
#' for each possible combination of predictors, i.e., all models that contain
#' exactly 1 predictor, all models that contain exactly 2 predictors, and so
#' forth. For each number of predictors, \code{\link[leaps]{regsubsets}} picks
#' the best model, i.e., the one with the smallest residual sum of squares (RSS)
#' or equivalently the largest \eqn{R^2}{R-squared}. This results in a best fit
#' for every possible number of predictors. \code{beset_lm} then uses
#' \code{k}-fold cross-validation to select the "best of the best": the best
#' model with the number of predictors that minimizes the prediction RSS or
#' equivalently maximizes the prediction \eqn{R^2}{R-squared}, i.e., how well
#' the best models trained using \eqn{k - 1} folds predict the out-of-fold
#' sample. \code{beset_lm} uses \code{\link[caret]{createFolds}} to randomly
#' assign observations to \code{k} folds within subgroups based on percentiles
#' of the numeric outcome. This insures that the folds will be matched in terms
#' of the outcome's frequency distribution. \code{beset_lm} also insures the
#' reproducibility of your analysis by requiring a \code{seed} to the random
#' number generator as one of its arguments.
#'
#' @section Warnings:
#' \enumerate{
#'  \item \code{beset_lm} handles missing data by performing listwise deletion.
#'   No other options for handling missing data are provided. The user is
#'   encouraged to deal with missing values prior to running this function.
#'  \item \code{beset_lm} is intended for use with additive linear models only.
#'    While there is no prohibition against the inclusion of interaction or
#'    polynomial terms, this practice is strongly discouraged. At best, this
#'    will result in an inefficient search because
#'    \code{\link[leaps]{regsubsets}} performs an exhaustive search over all
#'    possible variable subsets, including subsets that are hierarchically
#'    incomplete, i.e. subsets that contain an interaction term but are missing
#'    one or more of the subterms that comprise it. At worst, it may return one
#'    of these hierarchically incomplete models as the best model, an
#'    undesirable result if one cares about interpretability. If one wishes the
#'    model search to include interaction and/or non-linear effects, the
#'    \href{https://cran.r-project.org/web/packages/earth/index.html}{MARS}
#'    technique is recommended instead.
#'  \item \code{beset_lm} is best suited for searching over a small number of
#'    predictors (less than 10). For a large number of predictors (more than
#'    20), \code{\link{beset_elnet}} is recommended instead.
#' }
#'
#' @seealso \code{\link[caret]{createFolds}}, \code{\link[leaps]{regsubsets}},
#' \code{\link[base]{set.seed}}, \code{\link{predict_R2}}
#'
#' @param form A model \code{\link[stats]{formula}}.
#'
#' @param train_data A \code{\link[base]{data.frame}} with the variables in
#' \code{form} and the data to be used in model training.
#'
#' @param test_data Optional \code{\link[base]{data.frame}} with the variables
#' in \code{form} and the data to be used in model testing.
#'
#' @param n_folds Integer indicating the number of folds to use for
#' cross-validation.
#'
#' @param n_repeats Integer indicating the number of times cross-validation
#' should be repeated.
#'
#' @param seed An integer used to seed the random number generator when
#' assigning observations to folds.
#'
#' @return A list with the following three components:
#' \enumerate{
#'  \item\describe{
#'    \item{best_model}{an object of class \code{\link[stats]{lm}} corresponding
#'    to the best linear model with the number of parameters with the smallest
#'    cross-validation error (largest cross-validation R-squared)}
#'    }
#'  \item\describe{
#'    \item{best_model_1SE}{an object of class \code{\link[stats]{lm}}
#'    corresponding to the best linear model with the smallest number of
#'    parameters within one standard error of the smallest cross-validation
#'    error (largest cross-validation R-squared)}
#'    }
#'  \item \code{\link[base]{data.frame}} with the following columns:
#'    \describe{
#'    \item{k_folds}{the number of k-folds used}
#'    \item{n_preds}{the number of predictors in model}
#'    \item{R2_train}{the mean in-fold (training) R-squared for each size of
#'      best model}
#'    \item{R2_train_SE}{the standard error of the mean in-fold (training)
#'      R-squared for each size of best model}
#'    \item{R2_cv}{the mean out-of-fold (test) R-squared for each size of best
#'      model}
#'    \item{R2_cv_SE}{the standard error of the mean out-of-fold (test)
#'      R-squared for each size of best model}
#'    \item{R2_test}{if \code{test_data} is provided, the mean R-squared when
#'      the model fit to each subsample of the training data is applied to the
#'      independent test set}
#'    \item{R2_test_SE}{if \code{test_data} is provided, the standard error of
#'      the mean R-squared when the model fit to each subsample of the training
#'      data is applied to the independent test set}
#'    }
#'  }
#'
#' @export

beset_lm <- function(form, train_data, test_data = NULL, n_folds = 10,
                     n_repeats = 10, seed = 42)
{
  mf <- model.frame(form, data = train_data)
  colinear_vars <- caret::findLinearCombos(mf[,2:ncol(mf)])
  if(!is.null(colinear_vars$remove)){
    stop(paste(length(colinear_vars$remove), " linear dependencies found. ",
                  "Consider removing the following predictors:\n\t",
                  paste0(names(mf)[colinear_vars$remove], collapse = "\n\t"),
               sep = ""),
            immediate. = TRUE)
    mf <- mf[, -colinear_vars$remove]
  }
  p <- ncol(mf) - 1
  if(p > 20){
    stop(paste("`beset_lm` does not allow more than 20 predictors.",
               "Use `beset_elnet` instead."))
  }
  alt_p <- nrow(mf) - 1 - nrow(mf) %/% n_folds
  if(p > alt_p){
    warning(paste("You have more predictors than your sample size will support",
                  ".\n  Setting maximum subset size to ", alt_p, ".",
                  sep = ""), immediate. = TRUE)
    p <- alt_p
  }
  y <- mf[,1]
  response <- names(mf)[1]
  set.seed(seed)
  for(n in 1:n_repeats)
  {
    folds <- caret::createFolds(y, k = n_folds, list = FALSE)
    train_R2 <- matrix(nrow = n_folds, ncol = p)
    cv_R2 <- matrix(nrow = n_folds, ncol = p)
    test_R2 <- matrix(nrow = n_folds, ncol = p)
    for (i in 1:n_folds)
    {
      if(p > 1){
        best_fit <- leaps::regsubsets(form, data = mf[folds != i,], nvmax = p)
        fit_info <- summary(best_fit)
        train_R2[i,] <- fit_info$rsq
      } else {
        best_fit <- lm(form, data = mf[folds != i,])
        fit_info <- summary(best_fit)
        train_R2[i,] <- fit_info$r.squared
      }
      for (j in 1:p)
      {
        if(p > 1){
          best_preds <- mf[folds != i, fit_info$which[j,]]
          best_lm <- lm(paste(response, ".", sep = "~"), data = best_preds)
        } else {
          best_lm <- best_fit
        }
        cv_R2[i,j] <- predict_R2(best_lm, mf[folds == i,])
        if(!is.null(test_data)) test_R2[i,j] <- predict_R2(best_lm, test_data)
      }
    }
    R2_train <- apply(train_R2, 2, mean)
    R2_train_SE <- apply(train_R2, 2, function(x) sqrt(var(x)/length(x)))
    R2_cv <- apply(cv_R2, 2, mean)
    R2_cv_SE <- apply(cv_R2, 2, function(x) sqrt(var(x)/length(x)))
    R2_test <- apply(test_R2, 2, mean)
    R2_test_SE <- apply(test_R2, 2, function(x) sqrt(var(x)/length(x)))
    temp <- data.frame(k_folds = rep(n_folds, p),
                       n_preds = 1:p,
                       R2_train = R2_train,
                       R2_train_SE = R2_train_SE,
                       R2_cv = R2_cv,
                       R2_cv_SE = R2_cv_SE,
                       R2_test = R2_test,
                       R2_test_SE = R2_test_SE)
    if(n == 1) R2 <- temp else R2 <- rbind(R2, temp)
  }
  R2 <- group_by(R2, n_preds) %>%
    summarize_each(funs(mean))
  max_R2 <- max(R2$R2_cv)
  if(max_R2 <= 0){
    best_model <- lm(paste(response, "~ 1"), data = mf)
  } else if(p == 1){
    best_model <- lm(form, data = mf)
  } else {
    n_pred <- R2$n_preds[R2$R2_cv == max_R2]
    best_fit <- leaps::regsubsets(form, data = mf, nvmax = n_pred)
    fit_info <- summary(best_fit)
    best_preds <- mf[, fit_info$which[n_pred,]]
    best_model <- lm(paste(response, "~ ."), data = best_preds)
  }
  min_R2 <- max_R2 - R2$R2_cv_SE[R2$R2_cv == max_R2]
  if(min_R2 <= 0){
    best_model_1SE <- lm(paste(response, "~ 1"), data = mf)
  } else if(p == 1){
    best_model_1SE <- lm(form, data = mf)
  } else {
    temp <- R2[R2$R2_cv > min_R2,]
    n_pred <- min(temp$n_preds)
    best_fit <- leaps::regsubsets(form, data = mf, nvmax = n_pred)
    fit_info <- summary(best_fit)
    best_preds <- mf[, fit_info$which[n_pred,]]
    best_model_1SE <- lm(paste(response, ".", sep = "~"), data = best_preds)
  }

  structure(list(R2 = R2, best_model = best_model,
                 best_model_1SE = best_model_1SE),
            class = "beset_lm")
}
