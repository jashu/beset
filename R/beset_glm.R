#' Best Subset Selection for Generalized Linear Models
#'
#' \code{beset_glm} performs best subset selection using repeated
#' cross-validation to find the optimal number of predictors for several
#' families of generalized linear models.
#'
#' \code{beset_glm} performs best subset selection for generalized linear
#' models, fitting a separate model for each possible combination of predictors,
#' i.e., all models that contain exactly 1 predictor, all models that contain
#' exactly 2 predictors, and so forth. For each number of predictors,
#' \code{beset_glm} picks the model with the minimum cross-entropy (estimated
#' as the negative log of the expected probability of obtaining the observed
#' data given the model). This results in a best fit for every possible number
#' of predictors. \code{beset_glm} then uses \code{k}-fold cross-validation to
#' select the "best of the best": the best model with the number of predictors
#' that minimizes prediction cross-entropy, i.e., how well the best models
#' trained using \eqn{k - 1} folds predict the left-out fold. \code{beset_glm}
#' uses \code{\link[caret]{createFolds}} to randomly assign observations to
#' \code{k} folds within levels of the outcome when the outcome is a factor or
#' within subgroups based on percentiles when the outcome is numeric. This
#' insures that the folds will be matched in terms of the outcome's frequency
#' distribution. \code{beset_glm} also insures the reproducibility of your
#' analysis by requiring a \code{seed} to the random number generator as one of
#' its arguments.
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
#' @param train_data A \code{\link[base]{data.frame}} with the variables in
#' \code{form} and the data to be used in model training.
#'
#' @param test_data Optional \code{\link[base]{data.frame}} with the variables
#' in \code{form} and the data to be used in model testing.
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
#' @param n_folds Integer indicating the number of folds to use for
#' cross-validation.
#'
#' @param n_repeats Integer indicating the number of times cross-validation
#' should be repeated.
#'
#' @param n_cores Optional integer value indicating the number of workers to run
#' in parallel during subset search and cross-validation. By default, this will
#' be set to equal half the detectable cores on your machine. You may wish to
#' change this depending on your hardware and OS.
#' See \code{\link[parallel]{parallel-package}} for more information.
#'
#' @param seed An integer used to seed the random number generator when
#' assigning observations to folds.
#'
#' @return A "beset_glm" object with the following components:
#' \enumerate{
#'  \item\describe{
#'    \item{best_model}{an object of class \code{\link[stats]{glm}}
#'    corresponding to the best model with the number of parameters with the
#'    smallest cross-validation error}
#'    }
#'  \item\describe{
#'    \item{best_model_1SE}{an object of class \code{\link[stats]{glm}}
#'    corresponding to the best model with the smallest number of
#'    parameters within one standard error of the smallest cross-validation
#'    error}
#'    }
#'  \item\describe{
#'    \item{all_subsets}{a data frame containing fit statistics for every
#'      possible combination of predictors:
#'      \describe{
#'      \item{n_pred}{the number of predictors in model}
#'      \item{form}{formula for model}
#'      \item{train_CE}{Cross entropy between model predictions and
#'        \code{train_data}}
#'      \item{test_CE}{if \code{test_data} is provided, the cross entropy
#'        between the model fit to \code{train_data} and the \code{test_data}}
#'       }
#'    }
#'  }
#'  \item\describe{
#'    \item{best_subsets}{a data frame containing cross-validation statistics
#'      for the best model for each \code{n_pred} listed in \code{all_subsets}.
#'      In addition to the columns found in \code{all_subsets}, contains the
#'      following:
#'      \describe{
#'      \item{cv_CE}{the mean cross entropy between the predictions of models
#'        fit to \code{n-1} folds and the left-out fold}
#'      \item{cv_CE_SE}{the standard error of the mean cross entropy}
#'       }
#'    }
#'  }
#' }
#'
#' @import foreach
#' @export
beset_glm <- function(form, train_data, test_data = NULL,
                      family = "gaussian", link = NULL, ...,
                      p_max = 10, n_folds = 10, n_repeats = 10,
                      n_cores = NULL, seed = 42){
  #==================================================================
  # Check family argument and set up link function if specified
  #------------------------------------------------------------------
  family <- try(match.arg(family,
                          c("binomial", "gaussian", "poisson", "negbin")),
                silent = TRUE)
  if(class(family) == "try-error") stop("Invalid 'family' argument.")
  if(!is.null(link)){
    if(family == "binomial"){
      link <- try(match.arg(link, c("logit", "probit", "cauchit",
                                      "log", "cloglog")), silent = TRUE)
    } else if(family == "gaussian"){
      link <- try(match.arg(link, c("identity", "log", "inverse")),
                  silent = TRUE)
    } else if(family %in% c("negbin", "poisson")){
      link <- try(match.arg(link, c("log", "sqrt", "identity")),
                  silent = TRUE)
    }
    if(class(link) == "try-error")
      stop(paste("Invalid 'link' argument for", family, "family."))
  }
  if(family == "negbin"){
    if(is.null(link)) link <- "log"
  } else if(is.null(link)){
    dist <- call(family)
  } else {
    dist <- call(family, link = link)
  }

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
  # Check that number of predictors and cv folds is acceptable
  #------------------------------------------------------------------
  if(ncol(mf) > 21)
    stop("Best subsets not recommended for use with more than 20 predictors.")
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
  n_pred <- unlist(sapply(1:p, function(i) rep(i, choose(ncol(mf)-1, i))))
  n_pred <- c(0, n_pred)
  pred <- lapply(1:p, function(x)
    combn(names(mf)[2:ncol(mf)], x, simplify = FALSE))
  pred <- unlist(sapply(pred, function(vars){
    sapply(vars, function(x) paste0(x, collapse = " + "))
    }))
  pred <- c("1", pred)
  form_list <- paste(response, "~", pred)

  #======================================================================
  # Obtain cross entropy for every model
  #----------------------------------------------------------------------
  if(is.null(n_cores)) n_cores <- parallel::detectCores() %/% 2
  cl <- parallel::makeCluster(n_cores)
  parallel::clusterExport(cl, c("mf", "dist", "link", "test_data",
                                "glm_nb", "cross_entropy"),
                          envir=environment())
  if(family == "negbin"){
    catch <- parallel::clusterEvalQ(cl, library(MASS))
    CE <- parallel::parSapplyLB(cl, form_list, function(form){
      fit <- glm_nb(form, mf, link = link, ...)
      train_CE <- cross_entropy(fit, mf)
      test_CE <- NA_real_
      if(!is.null(test_data)) test_CE <- cross_entropy(fit, test_data)
      c(train_CE, test_CE)
    })
  } else {
    CE <- parallel::parSapplyLB(cl, form_list, function(form){
      fit <- glm(form, eval(dist), mf, ...)
      train_CE <- cross_entropy(fit, mf)
      test_CE <- NA_real_
      if(!is.null(test_data)) test_CE <- cross_entropy(fit, test_data)
      c(train_CE, test_CE)
    })
  }
  parallel::stopCluster(cl)
  all_subsets <- dplyr::data_frame(n_pred = n_pred,
                                   form = form_list,
                                   train_CE = CE[1,],
                                   test_CE = CE[2,])
  all_subsets <- dplyr::arrange(all_subsets, n_pred, train_CE)

  #======================================================================
  # Obtain model with best R^2 for each number of parameters
  #----------------------------------------------------------------------
  best_subsets <- dplyr::group_by(all_subsets, n_pred)
  best_subsets <- dplyr::filter(best_subsets, train_CE == min(train_CE))

  #======================================================================
  # Perform cross-validation on best models
  #----------------------------------------------------------------------
  search_grid <- expand.grid(fold = 1:n_folds, n_pred = 0:p)
  search_grid <- dplyr::left_join(search_grid, best_subsets, by = "n_pred")
  seed_seq <- seq.int(from = seed, length.out = n_repeats)
  doParallel::registerDoParallel()
  CE <- foreach(seed = seed_seq, .combine = rbind, .packages = "MASS") %dopar% {
    set.seed(seed)
    folds <- caret::createFolds(y, k = n_folds, list = FALSE)
    fits <- mapply(function(fold, form){
      if(family == "negbin"){
        fit <- glm_nb(form, mf[folds != fold,], link = link, ...)
        } else {
          fit <- glm(form, eval(dist), mf[folds != fold,], ...)
        }
      }, fold = search_grid$fold, form = search_grid$form, SIMPLIFY = FALSE)

    cv_CE <- matrix(mapply(function(fit, fold)
      cross_entropy(fit, mf[folds == fold,]),
      fit = fits, fold = search_grid$fold), nrow = n_folds, ncol = p + 1)
    CE_cv <- apply(cv_CE, 2, mean, na.rm = T)
    CE_cv_SE <- apply(cv_CE, 2, function(x) sqrt(var(x)/length(x)))

    data.frame(n_pred = 0:p,
               cv_CE = CE_cv,
               cv_CE_SE = CE_cv_SE)
  }
  #======================================================================
  # Derive cross-validation statistics
  #----------------------------------------------------------------------
  CE <- dplyr::group_by(CE, n_pred)
  mean_CE <- dplyr::summarize_each(CE, dplyr::funs(mean))

  #======================================================================
  # Fit best model and 1SE best model to full data set
  #----------------------------------------------------------------------
  best_subsets <- dplyr::left_join(best_subsets, mean_CE, by = "n_pred")
  best_subset <- which.min(best_subsets$cv_CE)
  best_form <- best_subsets$form[best_subset]
  max_CE <- min(best_subsets$cv_CE) + best_subsets$cv_CE_SE[best_subset]
  best_subsets_1SE <- best_subsets[best_subsets$cv_CE < max_CE,]
  best_form_1SE <- best_subsets_1SE$form[which.min(best_subsets_1SE$n_pred)]
  if(family == "negbin"){
    best_model <- glm_nb(best_form, mf, link = link, ...)
    if(best_form == best_form_1SE){
      best_model_1SE <- best_model
    } else {
      best_model_1SE <- glm_nb(best_form_1SE, mf, link = link, ...)
    }
  } else{
    best_model <- glm(best_form, eval(dist), mf, ...)
    if(best_form == best_form_1SE){
      best_model_1SE <- best_model
    } else {
      best_model_1SE <- glm(best_form_1SE, eval(dist), mf, ...)
    }
  }

  #======================================================================
  # Construct beset_glm object
  #----------------------------------------------------------------------
  structure(list(all_subsets = all_subsets, best_subsets = best_subsets,
                 best_model = best_model, best_model_1SE = best_model_1SE),
            class = "beset_glm")
}

#' @export
beset_lm <- function(form, train_data, test_data = NULL,
                     p_max = 10, n_folds = 10, n_repeats = 10,
                     n_cores = NULL, seed = 42){
  beset_glm(form, train_data, test_data = test_data, p_max = p_max,
            n_folds = n_folds, n_repeats = n_repeats, n_cores = n_cores,
            seed = seed)
}

glm_nb <- function (formula, data, weights, subset, na.action, start = NULL,
                    etastart, mustart, control = glm.control(...),
                    method = "glm.fit", model = TRUE, x = FALSE, y = TRUE,
                    contrasts = NULL, ..., init.theta, link = log){
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
  th <- as.vector(theta.ml(Y, mu, sum(w), w, limit = control$maxit,
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
    th <- theta.ml(Y, mu, sum(w), w, limit = control$maxit,
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

