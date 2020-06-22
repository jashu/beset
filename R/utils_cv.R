# Utility functions for cross-validation methods
#' @import dplyr
#' @import purrr

set_lm_par <- function(object, data){
  x <- object[["x"]]; y <- object[["y"]]; data <- object[["model"]];
  contrasts.arg <- object$contrasts; family_name <- object$family$family
  if(is.null(x) && !is.null(data)){
    x <- model.matrix(object$terms, data, contrasts.arg = contrasts.arg)
  }
  if(is.null(y) && !is.null(data)){
    y <-  purrr::as_vector(data[all.vars(object$terms)][1])
  }
  if(is.null(x) || is.null(y)){
    stop(paste(
      "Model data not found in model object.",
      "Use `data` argument to supply data frame to which model was fit"))
  }
  names(y) <- rownames(x)
  w <- object$weights
  if(is.null(w)) w <- rep(1, nrow(x))
  offset <- object$offset
  if(is.null(offset)) offset <- rep(0, nrow(x))
  list(x = x, y = y, w = w, offset = offset)
}

set_glm_par <- function(object, data){
  x <- object[["x"]]; y <- object[["y"]]; data <- object[["model"]];
  contrasts.arg <- object$contrasts; family_name <- object$family$family
  fitter <- "glm.fit"
  if(grepl("Negative Binomial", family_name)){
    family_name <- "negbin"
    fitter <- "glm_nb"
  }
  if(is.null(x) && !is.null(data)){
    x <- model.matrix(object$terms, data, contrasts.arg = contrasts.arg)
  }
  if(is.null(y) && !is.null(data)){
    y <-  purrr::as_vector(data[all.vars(object$terms)][1])
  }
  if(is.null(x) || is.null(y)){
    stop(paste(
      "Model data not found in model object.",
      "Use `data` argument to supply data frame to which model was fit"))
  }
  if(family_name == "binomial" && !is.factor(y)) y <- factor(y)
  names(y) <- rownames(x)
  weights <- object$prior.weights
  if(is.null(weights)) weights <- rep(1, nrow(x))
  offset <- object$offset
  if(is.null(offset)) offset <- rep(0, nrow(x))
  other_args <- list(family = object$family,
                     control = object$control,
                     intercept = hasName(object$coefficients, "(Intercept)"))
  list(x = x, y = y, weights = weights, offset = offset,
       other_args = other_args, family_name = family_name, fitter = fitter)
}

set_zi_par <- function(object){
  x <- model.matrix(object, "count")
  z <- model.matrix(object, "zero")
  data <- object$model
  y <- if(!is.null(object$y)){
    object$y
  } else if(!is.null(data)){
    model.response(data)
  } else {
    stop(
      paste(
        "Model data not found in model object.",
        "Refit `zeroinfl` and set either `model` or `y` argument to `TRUE`.")
    )
  }
  names(y) <- rownames(x)
  weights <- object$weights
  offset <- object$offset
  control <- object$control
  control$start <- object$start
  list(x = x, y = y, z = z, weights = weights, offset = offset,
       dist = object$dist, link = object$link, control = control)
}

set_cv_par <- function(n_obs, n_folds, n_reps, silent = FALSE,
                       nest_cv = FALSE){
  if(n_folds > n_obs) n_folds <- n_obs
  if(n_folds < n_obs){
    fold_size <- n_obs / n_folds
    if(nest_cv) fold_size <- fold_size * (n_folds - 1) / n_folds
    if(fold_size < 2){
      if(!silent){
        warning(
          paste(
            "Your sample size is too small to perform ", n_folds,
            "-fold cross-validation.\n",
            "  Performing leave-one-out cross-validation instead.", sep = ""
          ), immediate. = TRUE
        )
      }
      n_folds <- n_obs
    }
  }
  if(n_folds == n_obs && n_reps > 1){
    if(!silent){
      message(
        paste(
          "\nNOTE: Repetitions of leave-one-out cross-validation are pointless",
          "  and will not be performed.", sep = ""
        )
      )
    }
    n_reps <- 1
  }
  list(n_folds = n_folds, n_reps = n_reps)
}

get_cv_stats <-  function(y, y_hat, family, n_folds, n_reps, phi = NULL,
                          theta = NULL){
  fold_stats <- map(
    y_hat, ~ predict_metrics_(
      y = if(is.list(.x)) y[names(.x$mu)] else y[rownames(.x)],
      y_hat = if(is.list(.x)) .x$y_hat else .x,
      family = family,
      mu = if(is.list(.x)) .x$mu else NULL,
      phi = if(is.list(.x)) .x$phi else NULL,
      theta = theta
    )
  ) %>% transpose %>% simplify_all
  btwn_fold_error <- map(
    fold_stats, function(x){
      na_adjustment <- sum(is.na(x)) / n_reps
      sd(x, na.rm = TRUE) / sqrt(n_folds - na_adjustment)
    }
  )
  repeats <- paste("Rep", 1:n_reps, "$", sep = "")
  hold_out_pred <- map(repeats, ~ y_hat[grepl(.x, names(y_hat))])
  if(is.list(y_hat[[1]])){
    hold_out_mu <- hold_out_pred %>% map(~ map(.x, "mu")) %>%
      map(~ reduce(. , c)) %>% map(~.x[names(y)])
    hold_out_phi <- hold_out_pred %>% map(~ map(.x, "phi")) %>%
      map(~ reduce(. , c)) %>% map(~.x[names(y)])
    hold_out_yhat <- hold_out_pred %>% map(~ map(.x, "y_hat")) %>%
      map(~ reduce(. , c)) %>% map(~.x[names(y)])
    names(hold_out_mu) <- names(hold_out_phi) <- names(hold_out_yhat) <-
      paste("Rep", 1:n_reps, sep = "")
    hold_out_mu <- as_tibble(hold_out_mu)
    hold_out_phi <- as_tibble(hold_out_phi)
    hold_out_yhat <- as_tibble(hold_out_yhat)
    hold_out_pred <- list(
      mu = hold_out_mu,
      phi = hold_out_phi,
      yhat = hold_out_yhat
    )
    rep_stats <- map(
      1:n_reps,
      ~ predict_metrics_(
        y = y, y_hat = hold_out_yhat[[.x]], family = family,
        mu = hold_out_mu[[.x]], phi = hold_out_phi[[.x]], theta = theta
      )
    ) %>% transpose %>% simplify_all %>% as_tibble
  } else {
    hold_out_pred <- map(hold_out_pred, ~ reduce(. , rbind)) %>%
      map(~.x[names(y),])
    names(hold_out_pred) <- paste("Rep", 1:n_reps, sep = "")
    hold_out_pred <- as_tibble(hold_out_pred)
    rep_stats <- map(
      hold_out_pred, ~ predict_metrics_(
        y = y, y_hat = .x, family = family, phi = phi, theta = theta
      )
    ) %>% transpose %>% simplify_all %>% as_tibble
  }
  btwn_rep_range <- if(n_reps > 1){
    map(rep_stats, ~ range(.x, na.rm = TRUE))
  } else map(rep_stats, ~c(NA, NA))
  cv_means <- map(rep_stats, ~mean(.x, na.rm = TRUE))
  cv_stats <- list(mean = cv_means,
                   btwn_fold_se = btwn_fold_error,
                   btwn_rep_range = btwn_rep_range) %>% transpose()
  list(stats = cv_stats, predictions = hold_out_pred)
}

get_fold_ids <- function(fold_ids, n_reps){
  repeats <- paste("Rep", 1:n_reps, "$", sep = "")
  holdout_ids <- map(repeats, ~ fold_ids[grepl(.x, names(fold_ids))])
  holdout_fold <- map(
    holdout_ids, ~ reduce(
      map2(.x, seq_along(.x), ~rep(.y, length(.x))), c)) %>%
  map2(holdout_ids, ~ .x[order(reduce(.y, c))])
  names(holdout_fold) <- paste("Rep", 1:n_reps, sep = "")
  as_tibble(holdout_fold)
}

create_folds <- function(y, n_folds, n_reps, seed = 42){
  if(!is.factor(y) && length(unique(y)) <= 5) y <- factor(y)
  set.seed(seed)
  folds <- purrr::map(1:n_reps, ~ stratify_folds(y, n_folds))
  fold_ids <- expand.grid(Fold = 1:n_folds, Rep = 1:n_reps)
  fold_names <- paste("Fold", fold_ids$Fold, ".Rep", fold_ids$Rep, sep = "")
  fold_ids <- if(length(y) == n_folds){
    as.list(1:length(y))
  } else {
    purrr::map2(fold_ids$Fold, fold_ids$Rep, ~ which(folds[[.y]] == .x))
  }
  names(fold_ids) <- fold_names
  fold_ids
}

assign_folds <- function(n, n_folds = 10){
  fold_ids_1 <- rep(1:n_folds, n %/% n_folds)
  fold_ids_1 <- sample(fold_ids_1)
  modulus <- n - length(fold_ids_1)
  if(modulus > 0){
    fold_ids_2 <- sample.int(n_folds, size = modulus)
    c(fold_ids_1, fold_ids_2)
  } else fold_ids_1
}

stratify_folds <- function(y, n_folds = 10){
  folds <- vector("integer", length(y))
  if(is.numeric(y)){
    if(min(y) == 0 && sum(y == 0) >= n_folds){
      # make separate fold assignments for zero vs. non-zero values to insure
      # both 0 and count process is represented in every fold
      folds[y == 0] <- assign_folds(sum(y == 0), n_folds)
      folds[y != 0] <- stratify_folds(y[y != 0], n_folds)
      return(folds)
    } else {
      cuts <- length(y) %/% n_folds
      if(cuts < 2) cuts <- 2
      if(cuts > 5) cuts <- 5
      breaks <- unique(quantile(y, probs = seq(0, 1, length = cuts)))
      y <- as.integer(cut(y, breaks, include.lowest = TRUE))
    }
  }
  # Verify that each strata can appear in every fold
  if(all(table(y) >= n_folds)){
    purrr::walk(unique(y), function(x){
      folds[y == x] <<- assign_folds(sum(y == x), n_folds)
    })
  # Otherwise abort stratification and use simple randomization
  } else {
    folds <- assign_folds(length(y), n_folds)
  }
  folds
}
