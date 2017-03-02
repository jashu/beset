#' Best Subset Selection for Generalized Linear Models
#'
#' \code{beset_glm} performs best subset selection using repeated
#' cross-validation to find the optimal number of predictors for several
#' families of generalized linear models.
#'
#' \code{beset_glm} performs best subset selection for generalized linear
#' models, fitting a separate model for each possible combination of predictors
#' (all models that contain exactly 1 predictor, all models that contain
#' exactly 2 predictors, and so forth). For each number of predictors,
#' \code{beset_glm} picks the model with the maximum likelihood and then
#' estimates how well this model predicts new data using \code{k}-fold
#' cross-validation (how well a model trained using \eqn{k - 1} folds
#' predicts the left-out fold).
#'
#' @section Cross-validation details:
#' \code{beset_glm} uses \code{\link[caret]{createMultiFolds}} to randomly
#' partition the data set into \code{n_folds} * \code{n_repeats} folds within
#' strata (factor levels for factor outcomes percentile-based groups for numeric
#' outcomes). This insures that the folds will be matched in terms of the
#' outcome's frequency distribution. \code{beset_glm} also insures the
#' reproducibility of your analysis by requiring a \code{seed} to the random
#' number generator as one of its arguments.
#'
#' @section List of available families and link functions:
#' \describe{
#'  \item{\code{"gaussian"}}{The Gaussian family accepts the links
#'    \code{"identity"} (default), \code{"log"}, and \code{"inverse"}.}
#'  \item{\code{"binomial"}}{The binomial family accepts the links
#'    \code{"logit"} (default), \code{"probit"}, \code{"cauchit"}, \code{"log"}, and
#'    \code{"cloglog"} (complementary log-log).}
#'  \item{\code{"poisson"}}{The Poisson family accepts the links \code{"log"}
#'    (default), \code{"sqrt"}, and \code{"identity"}.}
#'  \item{\code{"negbin"}}{The negative binomial family accepts the links
#'    \code{"log"} (default), \code{"sqrt"}, and \code{"identity"}.}
#'  }
#'
#' @section Warnings:
#' \enumerate{
#'  \item \code{beset_glm} handles missing data by performing listwise deletion.
#'   No other options for handling missing data are provided at this time. The
#'   user is encouraged to deal with missing values prior to running this
#'   function.
#'  \item \code{beset_glm} is intended for use with additive models only.
#'    While there is no prohibition against the inclusion of interaction or
#'    polynomial terms, this practice is strongly discouraged. At best, this
#'    will result in an inefficient search because \code{beset_glm} performs an
#'    exhaustive search over all possible variable subsets, including subsets
#'    that are hierarchically incomplete, i.e., subsets that contain an
#'    interaction term but are missing one or more of the subterms that comprise
#'    it. At worst, it may return one of these hierarchically incomplete models
#'    as the best model, an undesirable result if one cares about
#'    interpretability. If one wishes the model search to include interaction
#'    and/or non-linear effects, the
#'    \href{https://cran.r-project.org/web/packages/earth/index.html}{MARS}
#'    technique is recommended instead. \item \code{beset_glm} is best suited
#'    for searching over a small number of predictors (less than 10). For a
#'    large number of predictors (more than 20), \code{\link{beset_elnet}} is
#'    recommended instead. However, note that \code{\link{beset_elnet}} only
#'    works with a more restricted set of distributions.
#' }
#'
#' @seealso \code{\link[caret]{createFolds}}, \code{\link[stats]{glm}},
#' \code{\link[base]{set.seed}}, \code{\link[MASS]{glm.nb}}
#'
#' @param form A model \code{\link[stats]{formula}}.
#'
#' @param train_data Data frame with the variables in \code{form} and the data
#' to be used for model fitting.
#'
#' @param test_data Optional data frame with the variables in \code{form} and
#' the data to be used for model validation.
#'
#' @param family Character string naming the error distribution to be used in
#' the model. Available families are listed under 'List of available families
#' and link functions'.
#'
#' @param link Optional character string naming the link function to be used in
#' the model. Available links and their defaults differ by \code{family} and are
#' listed under 'List of available families and link functions'.
#'
#' @param ... Additional arguments to be passed to \code{\link[stats]{glm}}
#'
#' @param p_max Maximum number of predictors to attempt to fit.
#'
#' @param n_cores Integer value indicating the number of workers to run in
#' parallel during subset search and cross-validation. By default, this will
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
#' @param seed An integer used to seed the random number generator when
#' assigning observations to folds.
#'
#' @return A "beset_glm" object with the following components:
#' \enumerate{
#'  \item\describe{
#'    \item{best_AIC}{an object of class \code{\link[stats]{glm}}
#'    corresponding to the model with the lowest Akaike Information Criterion}
#'    }
#'  \item\describe{
#'    \item{stats}{a list with three data frames:
#'      \describe{
#'        \item{fit}{statistics for every possible combination of predictors:
#'          \describe{
#'            \item{n_pred}{the total number of predictors in model; note that
#'               the number of predictors for a factor variable corresponds to the
#'               number of factor levels minus 1}
#'            \item{form}{formula for model}
#'            \item{AIC}{\eqn{-2*log-likelihood + k*npar}, where \eqn{npar}
#'              represents the number of parameters in the fitted model, and
#'              \eqn{k = 2}}
#'            \item{MCE}{Mean cross entropy, estimated as
#'              \eqn{-log-likelihood/N}, where \eqn{N} is the number of
#'              observations}
#'            \item{MSE}{Mean squared error}
#'            \item{R2}{R-squared, calculated as
#'              \eqn{1 - deviance/null deviance}}
#'            }
#'          }
#'      \item{cv}{a data frame containing cross-validation statistics
#'      for the best model for each \code{n_pred} listed in \code{fit_stats}.
#'      Each metric is computed using \code{\link{predict_metrics}}, with
#'      models fit to \eqn{n-1} folds and predictions made on the left-out fold.
#'      Each metric is followed by its standard error, estimated as the standard
#'      deviation of 1000 bootstrap replicates of computing the median cross-
#'      validation statistic across all folds and repetitions. The data frame
#'      is otherwise the same as that documented for \code{fit}, except
#'      AIC is omitted.}
#'      \item{test_stats}{if \code{test_data} is provided, a data frame
#'      containing prediction metrics for the best model for each \code{n_pred}
#'      listed in \code{fit} as applied to the \code{test_data}.}
#'       }
#'     }
#'   }
#' }
#' @export
beset_glm <- function(form, train_data, test_data = NULL, p_max = 10,
                      family = "gaussian", link = NULL,  ...,
                      n_cores = 2, n_folds = 10, n_repeats = 10, seed = 42){
  #==================================================================
  # Check family argument and set up link function if specified
  #------------------------------------------------------------------
  family <- tryCatch(match.arg(family, c("binomial", "gaussian", "poisson",
                                         "negbin")),
                     error = function(c){
                       c$message <- gsub("arg", "family", c$message)
                       c$call <- NULL
                       stop(c)
                     })
  if(!is.null(link)){
    link <- tryCatch(
      if(family == "binomial"){
        match.arg(link, c("logit", "probit", "cauchit", "log", "cloglog"))
        } else if(family == "gaussian"){
          match.arg(link, c("identity", "log", "inverse"))
          } else if(family %in% c("negbin", "poisson")){
            match.arg(link, c("log", "sqrt", "identity"))
          },
      error = function(c){
        c$message <- gsub("'arg'", paste("'link' for", family, "family"),
                          c$message)
        c$call <- NULL
        stop(c)
      })
  }
  if(family == "negbin"){
    if(is.null(link)) link <- "log"
  } else if(is.null(link)){
    dist <- call(family)
  } else {
    dist <- call(family, link = link)
  }
  #==================================================================
  # Create model frame and extract response name and vector
  #------------------------------------------------------------------
  mf <- model.frame(form, data = train_data, na.action = na.omit)
  if(!is.null(test_data)){
    test_data <- model.frame(form, data = test_data, na.action = na.omit)
  }
  n_drop <- nrow(train_data) - nrow(mf)
  if(n_drop > 0)
    warning(paste("Dropping", n_drop, "rows with missing data."),
            immediate. = TRUE)
  response <- names(mf)[1]
  y <- mf[,1]
  if(grepl("binomial", family)) y <- as.factor(y)

  #==================================================================
  # Screen for linear dependencies among predictors
  #------------------------------------------------------------------
  mm <- model.matrix(form, data = train_data)
  colinear_vars <- caret::findLinearCombos(mm[, 2:ncol(mm)])
  if(!is.null(colinear_vars$remove)){
    factor_idx <- which(sapply(mf, class) == "factor")
    factor_exp <- sapply(mf[, factor_idx], function(x) length(levels(x))) - 1
    mf_to_mm <- rep(1, ncol(mf))
    mf_to_mm[factor_idx] <- factor_exp
    mf_to_mm <- cumsum(mf_to_mm) - 1
    to_remove <- names(mf)[mf_to_mm %in% colinear_vars$remove]
    stop(paste(length(to_remove), " linear dependencies found. ",
               "Consider removing the following predictors:\n\t",
               paste0(to_remove, collapse = "\n\t"),
               sep = ""))
  }
  #==================================================================
  # Check that number of predictors and cv folds is acceptable
  #------------------------------------------------------------------
  p <- min(ncol(mf) - 1, p_max)
  if(family == "binomial"){
    n <- min(sum(y == levels(y)[1]), sum(y == levels(y)[2]))
  } else {
    n <- nrow(mf)
  }
  alt_folds <- n / (n - p * 10)
  alt_p <- p
  while(!dplyr::between(alt_folds, 1, 10)){
    alt_p <- alt_p - 1
    alt_folds <- n / (n - alt_p * 10)
  }
  if(alt_p < 1){
    if(family == "binomial"){
      stop("Your sample size for the minority class is too small.")
    } else {
      stop("Your sample size is too small.")
    }
  }
  if(alt_p < p){
    p <- alt_p
    warning(paste("'p_max' argument too high for your sample size",
                  ".\n  Reducing maximum subset size to ", p, ".",
                  sep = ""), immediate. = TRUE)
  }
  if(n_folds < alt_folds){
    n_folds <- as.integer(alt_folds)
    warning(paste("'n_folds' argument too low for your sample size ",
                  "and choice of 'p_max'",
                  ".\n  Increasing number of cv folds to ", n_folds, ".",
                  sep = ""), immediate. = TRUE)
  }
  #======================================================================
  # Make list of all possible formulas with number of predictors <= p_max
  #----------------------------------------------------------------------
  pred <- lapply(1:p, function(x)
    combn(names(mf)[2:ncol(mf)], x, simplify = FALSE))
  pred <- unlist(sapply(pred, function(vars){
    sapply(vars, function(x) paste0(x, collapse = " + "))}))
  pred <- c("1", pred)
  form_list <- paste(response, "~", pred)

  #======================================================================
  # Determine number of predictors for each model in list
  #----------------------------------------------------------------------
  n_pred <- sapply(pred[-1], function(x){
    x <- unlist(strsplit(x, split = " + ", fixed = TRUE))
    sum(sapply(x, function(y){
      if(is.factor(mf[,y]))
        length(levels(mf[,y])) - 1
      else
        1
    }))
  })
  n_pred <- c(intercept = 0, n_pred)

  #======================================================================
  # Obtain fit for every model
  #----------------------------------------------------------------------
  cl <- parallel::makeCluster(n_cores)
  negative.binomial <- MASS::negative.binomial
  parallel::clusterExport(cl, c("mf", "family", "dist", "link", "test_data",
                                "glm_nb", "prediction_metrics", ...,
                                "negative.binomial"), envir=environment())
  if(family == "negbin"){
    CE <- parallel::parLapplyLB(cl, form_list, function(form){
      fit <- glm_nb(form, mf, link = link, ...)
      list(AIC = AIC(fit),
           MCE = -logLik(fit)/nrow(mf),
           MSE = mean(residuals(fit, type = "response")^2),
           R2 = 1 - fit$deviance / fit$null.deviance)
    })
  } else {
    CE <- parallel::parLapplyLB(cl, form_list, function(form){
      fit <- glm(form, eval(dist), mf, ...)
      list(AIC = AIC(fit),
           MCE = -logLik(fit)/nrow(mf),
           MSE = mean(residuals(fit, type = "response")^2),
           R2 = 1 - fit$deviance / fit$null.deviance)
    })
  }
  fit_stats <- dplyr::data_frame(
    n_pred = n_pred,
    form = form_list,
    AIC = sapply(CE, function(x) x$AIC),
    MCE = sapply(CE, function(x) x$MCE),
    MSE = sapply(CE, function(x) x$MSE),
    R2 = sapply(CE, function(x) round(x$R2,3))
    )
  #======================================================================
  # Store the fit for the model with the best AIC
  #----------------------------------------------------------------------
  fit_stats <- dplyr::arrange(fit_stats, AIC)
  best_AIC <- if(family == "negbin"){
    glm_nb(fit_stats$form[1], mf, link = link, ...)
  } else{
    glm(fit_stats$form[1], eval(dist), mf, ...)
  }
  #======================================================================
  # Obtain model with maximum likelihood for each number of parameters
  #----------------------------------------------------------------------
  xval_stats <- dplyr::group_by(fit_stats, n_pred)
  xval_stats <- dplyr::filter(xval_stats, MCE == min(MCE))
  xval_stats <- dplyr::ungroup(xval_stats)
  xval_stats <- dplyr::select(xval_stats, n_pred, form)
  xval_stats <- dplyr::arrange(xval_stats, n_pred)

  #======================================================================
  # Perform cross-validation on best models
  #----------------------------------------------------------------------
  set.seed(seed)
  fold_ids <- caret::createMultiFolds(y, k = n_folds, times = n_repeats)
  metrics <- parallel::parLapply(cl, fold_ids, function(i, form_list){
    lapply(form_list, function(form){
      fit <- if(family == "negbin"){
        glm_nb(form, mf[i,], link = link, ...)
        } else {
          glm(form, eval(dist), mf[i,], ...)
        }
      prediction_metrics(fit, test_data = mf[-i,])
    })
  }, form_list = xval_stats$form)
  parallel::stopCluster(cl)
  #======================================================================
  # Derive cross-validation statistics
  #----------------------------------------------------------------------
  MCE <- sapply(metrics, function(models)
    sapply(models, function(x) x$mean_cross_entropy))
  xval_stats$MCE <- apply(MCE, 1, median, na.rm = T)
  xval_stats$MCE_SE <- apply(MCE, 1, function(x){
    MCE_boot <- boot::boot(x, function(x, i) median(x[i], na.rm = T), 1000)
    sd(MCE_boot$t)
  })
  MSE <- sapply(metrics, function(models)
    sapply(models, function(x) x$mean_squared_error))
  xval_stats$MSE <- apply(MSE, 1, median, na.rm = T)
  xval_stats$MSE_SE <- apply(MSE, 1, function(x){
    MSE_boot <- boot::boot(x, function(x, i) median(x[i], na.rm = T), 1000)
    sd(MSE_boot$t)
  })
  R2 <- sapply(metrics, function(models)
    sapply(models, function(x) x$R_squared))
  xval_stats$R2 <- apply(R2, 1, median, na.rm = T)
  xval_stats$R2_SE <- apply(R2, 1, function(x){
    R2_boot <- boot::boot(x, function(x, i) median(x[i], na.rm = T), 1000)
    sd(R2_boot$t)
  })

  #======================================================================
  # Compute prediction statistics for independent test set
  #----------------------------------------------------------------------
  test_stats <- NULL
  if(!is.null(test_data)){
    metrics <- lapply(form_list, function(form){
      fit <- if(family == "negbin")
        glm_nb(form, mf, link = link, ...)
      else
        glm(form, eval(dist), mf, ...)
      prediction_metrics(fit, test_data = test_data)
      })
    test_stats <- dplyr::select(xval_stats, n_pred, form)
    test_stats$MCE <- sapply(metrics, function(x) x$mean_cross_entropy)
    test_stats$MSE <- sapply(metrics, function(x) x$mean_squared_error)
    test_stats$R2 <- sapply(metrics, function(x) x$R_squared)
  }
  #======================================================================
  # Construct beset_glm object
  #----------------------------------------------------------------------
  structure(list(fit_stats = fit_stats, xval_stats = xval_stats,
                 test_stats = test_stats, best_AIC = best_AIC, model_data = mf,
                 xval_params = list(n_folds = n_folds, n_repeats = n_repeats,
                                    seed = seed)),
            class = "beset_glm")
}

#' @export
beset_lm <- function(form, train_data, test_data = NULL,
                     p_max = 10, n_folds = 10, n_repeats = 10,
                     n_cores = 2, seed = 42){
  beset_glm(form, train_data, test_data = test_data, p_max = p_max,
            n_folds = n_folds, n_repeats = n_repeats, n_cores = n_cores,
            seed = seed)
}

glm_nb <- function (formula, data, weights, subset, na.action, start = NULL,
                    etastart, mustart, control = glm.control(...),
                    method = "glm.fit", model = TRUE, x = FALSE, y = TRUE,
                    contrasts = NULL, ..., init.theta, link = log){
  negative.binomial <- MASS::negative.binomial
  loglik <- function(n, th, mu, y, w){
    sum(w * (lgamma(th + y) - lgamma(th) - lgamma(y + 1) + th * log(th) + y *
               log(mu + (y == 0)) - (th + y) * log(th + mu)))
  }
  fam0 <- if (missing(init.theta))
    do.call("poisson", list(link = link))
  else do.call("negative.binomial", list(theta = init.theta, link = link))
  mf <- Call <- match.call()
  m <- match(c("formula", "data", "subset", "weights", "na.action",
               "etastart", "mustart", "offset"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval.parent(mf)
  Terms <- attr(mf, "terms")
  if (method == "model.frame")
    return(mf)
  Y <- model.response(mf, "numeric")
  X <- if (!is.empty.model(Terms))
    model.matrix(Terms, mf, contrasts)
  else matrix(nrow = NROW(Y), ncol = 0)
  w <- model.weights(mf)
  if (!length(w))
    w <- rep(1, nrow(mf))
  else if (any(w < 0))
    stop("negative weights not allowed")
  offset <- model.offset(mf)
  mustart <- model.extract(mf, "mustart")
  etastart <- model.extract(mf, "etastart")
  n <- length(Y)
  if (!missing(method)) {
    if (!exists(method, mode = "function"))
      stop(gettextf("unimplemented method: %s", sQuote(method)),
           domain = NA)
    glm.fitter <- get(method)
  }
  else {
    method <- "glm.fit"
    glm.fitter <- stats::glm.fit
  }
  if (control$trace > 1)
    message("Initial fit:")
  fit <- glm.fitter(x = X, y = Y, w = w, start = start, etastart = etastart,
                    mustart = mustart, offset = offset, family = fam0,
                    control = list(maxit = control$maxit,
                                   epsilon = control$epsilon,
                                   trace = control$trace > 1),
                    intercept = attr(Terms, "intercept") > 0)
  class(fit) <- c("glm", "lm")
  mu <- fit$fitted.values
  th <- as.vector(MASS::theta.ml(Y, mu, sum(w), w, limit = control$maxit,
                           trace = control$trace > 2))
  if (control$trace > 1)
    message(gettextf("Initial value for 'theta': %f", signif(th)),
            domain = NA)
  fam <- do.call("negative.binomial", list(theta = th, link = link))
  iter <- 0
  d1 <- sqrt(2 * max(1, fit$df.residual))
  d2 <- del <- 1
  g <- fam$linkfun
  Lm <- loglik(n, th, mu, Y, w)
  Lm0 <- Lm + 2 * d1
  while ((iter <- iter + 1) <= control$maxit &&
         (abs(Lm0 - Lm)/d1 + abs(del)/d2) > control$epsilon) {
    eta <- g(mu)
    fit <- glm.fitter(x = X, y = Y, w = w, etastart = eta, offset = offset,
                      family = fam, control = list(maxit = control$maxit,
                                                   epsilon = control$epsilon,
                                                   trace = control$trace > 1),
                      intercept = attr(Terms, "intercept") > 0)
    t0 <- th
    th <- MASS::theta.ml(Y, mu, sum(w), w, limit = control$maxit,
                   trace = control$trace > 2)
    fam <- do.call("negative.binomial", list(theta = th, link = link))
    mu <- fit$fitted.values
    del <- t0 - th
    Lm0 <- Lm
    Lm <- loglik(n, th, mu, Y, w)
    if (control$trace) {
      Ls <- loglik(n, th, Y, Y, w)
      Dev <- 2 * (Ls - Lm)
      message(sprintf("Theta(%d) = %f, 2(Ls - Lm) = %f",
                      iter, signif(th), signif(Dev)), domain = NA)
    }
  }
  if (!is.null(attr(th, "warn")))
    fit$th.warn <- attr(th, "warn")
  if (iter > control$maxit) {
    warning("alternation limit reached")
    fit$th.warn <- gettext("alternation limit reached")
  }
  if (length(offset) && attr(Terms, "intercept")) {
    null.deviance <- if (length(Terms))
      glm.fitter(X[, "(Intercept)", drop = FALSE], Y, w, offset = offset,
                 family = fam, control = list(maxit = control$maxit,
                                              epsilon = control$epsilon,
                                              trace = control$trace > 1),
                 intercept = TRUE)$deviance
    else fit$deviance
    fit$null.deviance <- null.deviance
  }
  class(fit) <- c("negbin", "glm", "lm")
  fit$terms <- Terms
  fit$formula <- as.vector(attr(Terms, "formula"))
  Call$init.theta <- signif(as.vector(th), 10)
  Call$link <- link
  fit$call <- Call
  if (model)
    fit$model <- mf
  fit$na.action <- attr(mf, "na.action")
  if (x)
    fit$x <- X
  if (!y)
    fit$y <- NULL
  fit$theta <- as.vector(th)
  fit$SE.theta <- attr(th, "SE")
  fit$twologlik <- as.vector(2 * Lm)
  fit$aic <- -fit$twologlik + 2 * fit$rank + 2
  fit$contrasts <- attr(X, "contrasts")
  fit$xlevels <- .getXlevels(Terms, mf)
  fit$method <- method
  fit$control <- control
  fit$offset <- offset
  fit
}

