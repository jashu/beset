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

set_cv_par <- function(n_obs, n_folds, n_reps){
  out <- list(n_folds = n_folds, n_reps = n_reps)
  if(n_obs / n_folds < 2){
    message("Performing leave-one-out cross-validation")
    out$n_folds <- n_obs
    out$n_reps <- 1
  }
  if(n_reps > out$n_reps) message(
    paste("NOTE: Repetitions of leave-one-out cross-validation\n",
          "are pointless and will not be performed.", sep = ""))
  out
}

get_cv_stats <-  function(y, y_hat, family, n_folds, n_reps,
                          phi = NULL, theta = NULL){
  fold_stats <- map(
    y_hat, ~ predict_metrics_(y = y[rownames(.x)], y_hat = .x, family = family,
                             phi = phi, theta = theta)) %>%
    transpose %>% simplify_all
  btwn_fold_error <- map(
    fold_stats, function(x){
      na_adjustment <- sum(is.na(x)) / n_reps
      sd(x, na.rm = TRUE) / sqrt(n_folds - na_adjustment)
    })
  repeats <- paste("Rep", 1:n_reps, "$", sep = "")
  hold_out_pred <- map(
    repeats, ~ y_hat[grepl(.x, names(y_hat))]) %>%
    map(~ reduce(. , rbind)) %>% map(~.x[names(y),])
  names(hold_out_pred) <- paste("Rep", 1:n_reps, sep = "")
  hold_out_pred <- as_data_frame(hold_out_pred)
  rep_stats <- map(hold_out_pred,
                   ~ predict_metrics_(y = y, y_hat = .x, family = family,
                                      phi = phi, theta = theta)) %>%
    transpose %>% simplify_all %>% as_data_frame
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
  as_data_frame(holdout_fold)
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

assign_folds <- function(y, n_folds = 10){
  fold_ids_1 <- rep(1:n_folds, length(y) %/% n_folds)
  fold_ids_1 <- sample(fold_ids_1)
  modulus <- length(y) - length(fold_ids_1)
  if(modulus > 0){
    fold_ids_2 <- sample.int(n_folds, size = modulus)
    c(fold_ids_1, fold_ids_2)
  } else fold_ids_1
}

stratify_folds <- function(y, n_folds = 10){
  folds <- vector("integer", length(y))
  if(is.numeric(y)){
    if(min(y) == 0){
      # make separate fold assignments for zero vs. non-zero values to insure
      # both 0 and count process is represented in every fold
      folds[y == 0] <- assign_folds(y[y == 0], n_folds)
      folds[y != 0] <- stratify_folds(y[y != 0], n_folds)
    } else {
      cuts <- length(y) %/% n_folds
      if(cuts < 2) cuts <- 2
      if(cuts > 5) cuts <- 5
      breaks <- unique(quantile(y, probs = seq(0, 1, length = cuts)))
      y <- as.integer(cut(y, breaks, include.lowest = TRUE))
    }
  }
  purrr::walk(unique(y), function(x){
    folds[y == x] <<- assign_folds(y[y == x], n_folds)
  })
  folds
}
