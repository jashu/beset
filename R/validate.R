#' Cross-Validated Prediction Metrics
#'
#' \code{validate} is a generic function for cross-validating predictions from
#' the results of various model fitting functions. The function invokes
#' particular \code{\link[utils]{methods}} which depend on the
#' \code{\link[base]{class}} of the first argument.
#'
#' To obtain cross-validation statistics, first fit a model as you
#' normally would and then pass the model object to \code{validate}.
#'
#' @return a "\code{cross_valid}" object consisting of a list with the
#' following elements:
#' \describe{
#'   \item{cv_stats}{a list of cross-validated prediction metrics, each
#'   containing the mean, between-fold standard error ("btwn_fold_se"), and
#'   between-repetition min-max range ("btwn_rep_range") of each metric.}
#'   \item{predictions}{a data frame containing the
#'   hold-out predictions for each row in the training data, with a separate
#'   column for each repetition of the k-fold cross-validation}
#'   \item{fold_assignments}{a data frame of equal
#'   dimensions to \code{predictions} giving the number of the hold-out fold
#'   of the corresponding element in \code{predictions}}
#'   \item{parameters}{a list documenting how many folds and repetitions were used
#'   for cross-validation, and the seed passed to the random number generator,
#'   which will be needed to reproduce the random fold assignments}
#'   }
#'
#' @param object A model object for which a cross-validated R-squared is
#' desired.
#'
#' @param data Data frame that was used to train the model. Only needed if
#' the training data is not contained in the model \code{object}.
#'
#' @param x Model matrix that was used to train elastic net.
#'
#' @param y Response variable that was used to train elastic net.
#'
#' @param lambda \code{Numeric} value of the penalty parameter \code{lambda}
#' at which predictions are cross-validated. Default is the entire sequence
#' used to create the model, but importance scores will only be generated
#' if a single \code{lambda} value is given.
#'
#' @param offset \code{Numeric} offset vector that was used to train elastic
#' net.
#'
#' @param weights = \code{Numeric} offset vector that was used to train
#' elastic net.
#'
#' @param n_folds \code{Integer} indicating the number of folds to use for
#' cross-validation.
#'
#' @param n_reps \code{Integer} indicating the number of times cross-validation
#' should be repeated.
#'
#' @param seed \code{Integer} used to seed the random number generator.
#'
#' @param envir \code{Environment} in which to look for variables listed in
#' object call
#'
#' @inheritParams predict_metrics_
#' @inheritParams beset_glm
#'
#' @import purrr
#' @import parallel
#' @import dplyr
#' @importFrom randomForest randomForest

#' @export
validate <- function(object, ...){
  UseMethod("validate")
}

validate.function <- function(object, x, y, family = "gaussian",
                              phi = NULL, theta = NULL,
                              n_folds = 10, n_reps = 10, seed = 42, ...){
  names(y) <- names(x)
  n_obs <- length(y)
  cv_par <- set_cv_par(n_obs, n_folds, n_reps)
  n_folds <- cv_par$n_folds; n_reps <- cv_par$n_reps
  fold_ids <- create_folds(y = y, n_folds = n_folds, n_reps = n_reps,
                           seed = seed)
  train_data <- lapply(fold_ids, function(i)
    c(
      list(x = x[-i, , drop = FALSE],
           y = y[-i]),
      list(...)
      )
  )
  cv_fits <- map(train_data, ~ do.call(object, .x))
  test_data <- lapply(fold_ids, function(i)
    c(
      list(x = x[i, , drop = FALSE],
           y = y[i])
      )
  )
  y_hat <- map2(cv_fits, test_data,
                ~ predict(.x, newdata = .y$x, type = "response"))
  cv_stats <- get_cv_stats(y = y, y_hat = y_hat, family = family,
                           n_folds = n_folds, n_reps = n_reps,
                           phi = phi, theta = theta)
  fold_assignments <- get_fold_ids(fold_ids, n_reps)
  structure(
    c(cv_stats, list(
      fold_assignments = fold_assignments,
      parameters = list(family = family,
                        n_obs = n_obs,
                        n_folds = n_folds,
                        n_reps = n_reps,
                        seed = seed,
                        y = y))),
    class = "cross_valid")
}

#' @describeIn validate Cross-validation of linear models
#' @export
validate.lm <- function(object, data = NULL, n_folds = 10, n_reps = 10,
                        seed = 42, ...){
  lm_par <- set_lm_par(object, data)
  y <- lm_par$y
  n_obs <- length(y)
  cv_par <- set_cv_par(n_obs, n_folds, n_reps)
  n_folds <- cv_par$n_folds; n_reps <- cv_par$n_reps
  fold_ids <- create_folds(y = y, n_folds = n_folds, n_reps = n_reps,
                           seed = seed)
  train_data <- lapply(fold_ids, function(i)
    list(x = lm_par$x[-i, , drop = FALSE],
         y = lm_par$y[-i],
         w = lm_par$w[-i],
         offset = lm_par$offset[-i])
  )
  cv_fits <- map(train_data, ~ do.call(lm.wfit, .x))
  test_data <- lapply(fold_ids, function(i)
    list(x = lm_par$x[i, , drop = FALSE],
         y = lm_par$y[i],
         offset = lm_par$offset[i])
  )
  y_hat <- map2(cv_fits, test_data, ~ .y$x %*% coef(.x) + .y$offset)
  cv_stats <- get_cv_stats(y = y, y_hat = y_hat, family = "gaussian",
                           n_folds = n_folds, n_reps = n_reps)
  fold_assignments <- get_fold_ids(fold_ids, n_reps)
  structure(
    c(cv_stats, list(
      fold_assignments = fold_assignments,
      parameters = list(family = "gaussian",
                        n_obs = n_obs,
                        n_folds = n_folds,
                        n_reps = n_reps,
                        seed = seed,
                        y = y))),
    class = "cross_valid")
}

#' @describeIn validate Cross-validation of GLMs
#' @export
validate.glm <- function(object, data = NULL, n_folds = 10, n_reps = 10,
                         seed = 42, ...){
  glm_par <- set_glm_par(object, data)
  y <- glm_par$y
  n_obs <- length(y)
  cv_par <- set_cv_par(n_obs, n_folds, n_reps)
  n_folds <- cv_par$n_folds; n_reps <- cv_par$n_reps
  fold_ids <- create_folds(y = y, n_folds = n_folds, n_reps = n_reps,
                           seed = seed)
  if(is.factor(y)) y <- as.integer(y) - 1L
  names(y) <- rownames(glm_par$x)
  train_data <- lapply(fold_ids, function(i)
    c(list(x = glm_par$x[-i, , drop = FALSE],
           y = glm_par$y[-i],
           weights = glm_par$weights[-i],
           offset = glm_par$offset[-i]),
      glm_par$other_args)
  )
  cv_fits <- map(train_data, ~ do.call(glm_par$fitter, .x))
  test_data <- lapply(fold_ids, function(i)
    list(x = glm_par$x[i, , drop = FALSE],
         y = glm_par$y[i],
         weights = glm_par$weights[i],
         offset = glm_par$offset[i])
  )
  y_hat <- map2(cv_fits, test_data,
                ~ .x$family$linkinv(.y$x %*% coef(.x) + .y$offset))
  cv_stats <- get_cv_stats(y = y, y_hat = y_hat, family = glm_par$family_name,
                           n_folds = n_folds, n_reps = n_reps,
                           theta = object$theta)
  fold_assignments <- get_fold_ids(fold_ids, n_reps)
  structure(
    c(cv_stats, list(
      fold_assignments = fold_assignments,
      parameters = list(family = glm_par$family_name,
                        n_obs = n_obs,
                        n_folds = n_folds,
                        n_reps = n_reps,
                        seed = seed,
                        y = y))),
    class = "cross_valid")
}

#' @describeIn validate Cross-validation of GLM nets
#' @export
validate.glmnet <- function(object, x = NULL, y = NULL, lambda = NULL,
                            offset = NULL, weights = NULL,
                            n_folds = 10, n_reps = 10, seed = 42,
                            envir = .GlobalEnv, ...){
  arg_list <- as.list(object$call)
  if(is.null(x)){
    x <- eval(arg_list$x, envir = envir)
  }
  n_obs <- nrow(x)
  if(is.null(y)){
    y <- eval(arg_list$y, envir = envir)
  }
  if(is.null(x) || is.null(y))
    stop("You must supply the original `x` and `y` used to train glmnet.")
  if(!is.null(lambda) && length(lambda) > 1){
    stop("Only one value of lambda may be specified.")
  }
  offset <- if(is.null(offset) && object$offset){
    offset <- eval(arg_list$offset, envir = envir)
  }
  if(is.null(weights)){
    weights <- if(hasName(arg_list, "weights"))
      eval(arg_list$weights, envir = envir) else rep(1, nrow(x))
  }

  data <- list(x = x, y = y, offset = offset, weights = weights)
  other_args <- arg_list[setdiff(names(arg_list), c("", names(data)))]
  other_args$lambda <- object$lambda
  params <- set_cv_par(n_obs, n_folds, n_reps)
  n_folds <- params$n_folds; n_reps <- params$n_reps
  fold_ids <- create_folds(y = y, n_folds = n_folds, n_reps = n_reps,
                           seed = seed)
  train_data <- lapply(fold_ids, function(i)
    c(list(x = data$x[-i, , drop = FALSE]), purrr::map(data[-1], ~.x[-i]),
      other_args)
  )
  cv_fits <- map(train_data, ~ suppressMessages(do.call(glmnet::glmnet, .x)))
  if(is.factor(y)) y <- as.integer(y) - 1L
  names(y) <- rownames(x)
  test_data <- lapply(fold_ids, function(i)
    c(list(x = data$x[i, , drop = FALSE]),
      purrr::map(data[-1], ~.x[i]))
  )
  y_hat <- map2(cv_fits, test_data,
                ~ predict(object = .x, newx = .y$x, s = lambda,
                          type = "response", newoffset = .y$offset)
                )
  family <- if(hasName(arg_list, "family"))
    eval(arg_list$family) else "gaussian"
  y_hat <- map(y_hat, function(x){
    pad <- length(lambda) - ncol(x)
    if(pad > 0) cbind(x, matrix(mean(y), ncol = pad, nrow = nrow(x))) else
      x
  })
  fold_stats <- if(is.null(lambda)){
   map(y_hat, ~ apply(.x, 2, function(x)
      predict_metrics_(y = y[rownames(.x)], y_hat = x, family = family))) %>%
    map(transpose) %>% transpose %>% map(transpose) %>% map(simplify_all)
  } else {
    map(y_hat, ~predict_metrics_(y = y[rownames(.x)], y_hat = .x,
                                        family = family)) %>%
      transpose %>% simplify_all
  }
  btwn_fold_error <- if(is.null(lambda)){
    map(fold_stats, ~ map_dbl(., ~ if(is.numeric(.))
      sd(., na.rm = TRUE)/sqrt(n_folds) else NA_real_)
    )} else {
    map(fold_stats, ~ if(is.numeric(.)){
      sd(., na.rm = TRUE)/sqrt(n_folds)
    } else NA_real_)
  }
  repeats <- paste("Rep", 1:n_reps, "$", sep = "")
  y_hat <- map(repeats, ~ y_hat[grepl(.x, names(y_hat))]) %>%
    map(~ reduce(. , rbind)) %>% map(~.x[names(y),])
  names(y_hat) <- paste("Rep", 1:n_reps, sep = "")
  rep_stats <- if(is.null(lambda)){
    map(y_hat, ~ apply(.x, 2, function(x)
      predict_metrics_(y = y, y_hat = x, family = family))
    ) %>% map(transpose) %>% transpose %>% map(transpose) %>% map(simplify_all)
  } else map(y_hat, ~ predict_metrics_(y = y, y_hat = .x, family = family)) %>%
    transpose %>% simplify_all
  btwn_rep_range <- NULL
  if(n_reps > 1){
    btwn_rep_range <- if(is.null(lambda)){
      map(rep_stats, ~ map(., range))
      } else {
      map(rep_stats, range)
      }
  }
  if(is.null(btwn_rep_range)){
    btwn_rep_range <- map(btwn_fold_error,
                          ~ map(.x, ~ return(c(NA, NA))))
  }
  cv_means <- if(is.null(lambda)) {
    map(rep_stats, ~ map_dbl(., ~ suppressWarnings(mean(., na.rm = TRUE))))
  } else map(rep_stats, ~ mean(., na.rm = TRUE))
  cv_stats <- list(mean = cv_means, btwn_fold_se = btwn_fold_error,
                   btwn_rep_range = btwn_rep_range) %>% transpose
  if(is.null(lambda)) cv_stats <- cv_stats %>% map(transpose)
  cv_stats <- list(stats = cv_stats,
                   predictions = y_hat)
  fold_assignments <- get_fold_ids(fold_ids, n_reps)
  structure(
  c(cv_stats, list(
    fold_assignments = fold_assignments,
    parameters = list(family = family,
                      n_obs = n_obs,
                      n_folds = n_folds,
                      n_reps = n_reps,
                      seed = seed,
                      y = y))),
  class = "cross_valid")
}

#' @describeIn validate Extract test error estimates from "beset_elnet"
#' objects with nested cross-validation
#' @export
validate.nested <- function(object,
                            metric = "auto",
                            oneSE = TRUE, ...){
  metric <- tryCatch(
    match.arg(metric, c("auto", "auc", "mae", "mce", "mse", "rsq")),
    error = function(c){
      c$message <- gsub("arg", "metric", c$message)
      c$call <- NULL
      stop(c)
    }
  )
  tryCatch(
    if(
      (metric == "auc" && object$family != "binomial") ||
      (metric == "mae" && object$family == "binomial")
    ) error = function(c){
      c$message <- paste(metric, "not available for", object$family, "models")
      c$call <- NULL
      stop(c)
    }
  )
  if(metric == "auto"){
    metric <- if(object$family == "gaussian") "mse" else "mce"
  }
  family <- object$family
  fold_ids <- object$fold_assignments
  n_folds <- object$n_folds
  n_reps <- object$n_reps
  repeats <- paste("Rep", 1:n_reps, "$", sep = "")
  rep_idx <- map(repeats, ~ grepl(.x, names(object$beset)))
  best_models <- map(
    object$beset, ~
      if(inherits(.x, "elnet")){
        beset:::get_best.beset_elnet(.x,  metric = metric, oneSE = oneSE)
      } else {
        beset:::get_best.beset_glm(.x,  metric = metric, oneSE = oneSE)
      }
  )
  if(inherits(object, "elnet")){
    alphas <- map_dbl(best_models, "alpha")
    lambdas <- map_dbl(best_models, "best_lambda")
    best_idx <- pmap_int(
      list(s = map(object$beset, ~.x$stats$test), a = alphas, l = lambdas),
      function(s, a, l) with(s, which(alpha == a & lambda == l)))
  } else {
    best_form <- map_chr(best_models, "formula")
    best_idx <- map2_int(map(object$beset, ~.x$stats$test), best_form,
                         ~ which(.x$form == .y))
  }
  test_stats <- map2_df(object$beset, best_idx, ~.x$stats$test[.y,])
  tt <- terms(object)
  y0 <-  model.response(model.frame(tt, object$data))
  Terms <- delete.response(tt)
  m <- model.frame(Terms, object$data, xlev = object$xlevels)
  X <- model.matrix(Terms, m, contrasts.arg = object$contrasts)
  if(inherits(object, "elnet") && "(Intercept)" %in% colnames(X)) X <- X[,-1]
  newoffset <- object$offset
  if(is.null(newoffset)){
    newoffset <- rep(0, length(y0))
  } else {
    names(newoffset) <- rownames(object$data)
    newoffset <- newoffset[names(y)]
  }

  y <- map(fold_ids, ~ y0[.x])
  y <- map(rep_idx, ~ y[.x]) %>% map(as_vector)
  y_hat <-
    if(inherits(object, "elnet")){
      map2(best_models, fold_ids,
                ~ predict(.x, newx = X[.y,], s = .x$best_lambda,
                          type = "response", newoffset = newoffset[.y]))
    } else {
      map2(
        best_models, fold_ids,
        ~ .x$family$linkinv(X[.y, names(coef(.x))] %*% coef(.x) + newoffset[.y]))
    }
  y_hat <- map(rep_idx, ~ y_hat[.x]) %>% map(simplify_all) %>% map(as_vector)
  theta <- NULL
  if(family == "negbin"){
    theta <- map_dbl(best_models, "theta") %>% mean
  }
  rep_stats <- map2(
    y, y_hat,
    ~ predict_metrics_(y = .x, y_hat = .y, family = family, theta = theta)
  ) %>% transpose %>% simplify_all %>% as_data_frame
  stat_names <- intersect(names(test_stats),
                          c("auc", "mae", "mce", "mse", "rsq"))
  test_stats <- list(
    mean = map(rep_stats, ~ mean(.x, na.rm = TRUE)),
    btwn_fold_se = map(test_stats[stat_names],
                       ~ sd(.x, na.rm = TRUE) / sqrt(n_folds)),
    btwn_rep_range =  map(rep_stats, ~ range(.x, na.rm = TRUE))
  ) %>% transpose

  fold_assignments <- get_fold_ids(fold_ids, n_reps)
  fold_ids <- map(rep_idx, ~ fold_ids[.x]) %>% map(as_vector)
  names(fold_ids) <- names(y_hat) <- names(fold_assignments)
  structure(
    list(
      stats = test_stats,
      predictions = map2_df(y_hat, fold_ids, ~ .x[order(.y)]),
      fold_assignments = fold_assignments,
      parameters = list(family = family,
                        n_obs = nrow(fold_assignments),
                        n_folds = n_folds,
                        n_reps = n_reps,
                        seed = object[[1]]$parameters$cv$seed,
                        y = y0)
      ), class = "cross_valid"
  )
}

#' @export
#' @describeIn validate Cross-validation of random forests
validate.randomForest <- function(object, data = NULL, x = NULL,
                                  n_folds = 10, n_reps = 10,
                                  seed = 42, ..., parallel_type = NULL,
                                  n_cores = NULL, cl = NULL){
  rf_par <- as.list(object$call)
  rf_par$keep.forest = TRUE
  other_args <- setdiff(names(rf_par), c("", "formula", "data"))
  if(is.null(data)){
    data <- eval(rf_par$data)
  }
  if(is.null(x)){
    x <- eval(rf_par$x)
  }
  if(is.null(data) && is.null(x)){
    stop("You must supply the original `data` or `x` used to train the random
         forest.")
  }
  if(is.null(n_cores) || n_cores > 1){
    parallel_control <- setup_parallel(
      parallel_type = parallel_type, n_cores = n_cores, cl = cl)
    have_mc <- parallel_control$have_mc
    n_cores <- parallel_control$n_cores
    cl <- parallel_control$cl
  }
  if(is.null(x)) x <- data[labels(object$terms)]
  y <- object$y
  n_obs <- length(y)
  cv_par <- set_cv_par(n_obs, n_folds, n_reps)
  n_folds <- cv_par$n_folds; n_reps <- cv_par$n_reps
  fold_ids <- create_folds(y = y, n_folds = n_folds, n_reps = n_reps,
                           seed = seed)
  if(is.factor(y)) y <- as.integer(y) - 1L
  train_data <- lapply(fold_ids, function(i)
    c(list(x = x[-i, , drop = FALSE],
           y = y[-i],
           xtest = x[i, , drop = FALSE],
           ytest = y[i]),
      rf_par[other_args])
  )
  cv_fits <- if(n_cores > 1L){
    if(have_mc){
      parallel::mclapply(train_data, function(x){
        set.seed(seed); do.call(randomForest, x)
      }, mc.cores = n_cores)
    } else {
      parallel::parLapply(cl, train_data, function(x){
        set.seed(seed); do.call(randomForest, x)
      })
    }
  } else {
    lapply(train_data, function(x){
      set.seed(seed); do.call(randomForest, x)
    })
  }

  y_hat <- map2(cv_fits, train_data, ~ as.matrix(predict(.x, .y$xtest)))
  for(i in seq_along(y_hat)) rownames(y_hat[[i]]) <- fold_ids[[i]]
  family <- ifelse(n_distinct(y) == 2, "binomial", "gaussian")
  cv_stats <- get_cv_stats(y = y, y_hat = y_hat, family = family,
                           n_folds = n_folds, n_reps = n_reps)
  fold_assignments <- get_fold_ids(fold_ids, n_reps)
  structure(
    c(cv_stats, list(
      fold_assignments = fold_assignments,
      parameters = list(family = family,
                        n_obs = n_obs,
                        n_folds = n_folds,
                        n_reps = n_reps,
                        seed = seed,
                        y = y))),
    class = "cross_valid")
}

#' @export
#' @describeIn validate Retrieve cross-validation metrics from "beset_rf" object
validate.beset_rf <- function(object, ...){
  object$stats
}
