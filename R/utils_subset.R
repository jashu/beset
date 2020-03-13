#' @import purrr
#' @import dplyr

#============================================================================
# Create covariate matrices and response vectors for use with beset functions
#----------------------------------------------------------------------------
make_args <- function(form, data, family, link, contrasts,
                      weights = NULL, offset = NULL,
                      start = NULL, etastart = NULL, mustart = NULL,
                      epsilon, maxit, check = TRUE, ...){
  if(inherits(data, "data_partition")){
    return(
      list(
        train = make_args(form, data$train, family, link, contrasts, weights,
                          offset, start, etastart, mustart, epsilon, maxit,
                          check, ...)$train,
        test = make_args(form, data$test, family, link, contrasts, weights,
                         offset, start, etastart, mustart, epsilon, maxit,
                         check = FALSE, ...)$train
      )
    )
  }
  link <- if(!is.null(link)){
    check_link(family, link)
  } else {
    switch(family,
           binomial = "logit",
           gaussian = "identity",
           poisson = "log",
           negbin = "log")
  }
  mf <- if(is.null(attr(data, "terms"))){
    check_names(names(data))
    environment(form) <- environment()
    if(is.null(offset)) offset <- rep(0, nrow(data))
    if(is.null(weights)) weights <- rep(1, nrow(data))
    model.frame(form, data, na.action = na.omit, offset = offset,
                weights = weights)
  } else {
    data
  }
  x <- stats::model.matrix(form, mf, contrasts)
  todrop <- grep("`\\(offset\\)`|`\\(weights\\)`", colnames(x))
  if(length(todrop)) x <- x[, -todrop]
  y <- model.response(mf)
  # insure y is a factor if family is binomial
  if(family == "binomial" && !is.factor(y)) y <- factor(y)
  family <- if(family == "negbin"){
    stats::poisson(link = link)
  } else {
    do.call(family, list(link = link))
  }
  control <- glm.control(epsilon = epsilon, maxit = maxit)
  offset <- model.offset(mf)
  weights <- model.weights(mf)
  if("alpha" %in% names(list(...))){
    list(
      train = c(
        list(x = x, y = y, family = family$family, weights = weights,
             offset = offset, intercept = "(Intercept)" %in% colnames(x)),
        list(...)
      ),
      test = NULL
    )
  } else {
    list(
      train = list(
        x = x, y = y, weights = weights, start = start, etastart = etastart,
        mustart = mustart, offset = offset, family = family, control = control,
        intercept = "(Intercept)" %in% colnames(x)
      ),
      test = NULL
    )
  }
}

#======================================================================
# Make list of all possible formulas with number of predictors <= p_max
#----------------------------------------------------------------------

get_subsets <- function(m, force_in, p_max){
  factors <- attr(m$terms, "dataClasses") == "factor"
  factors <- names(factors[factors])
  factors <- factors[factors %in% force_in]
  if(length(factors) > 0){
    dummies <- map(
      factors, ~ paste(.x, m$xlevels[[.x]][-1], sep = "")
    ) %>% as_vector
    force_in <- c(setdiff(force_in, factors), dummies)
  }
  force_in <- c("(Intercept)", force_in)
  forced <- map(force_in, ~ colnames(m$train$x) == .) %>% transpose %>%
    map_lgl(~ reduce(., `|`))
  p <- min(ncol(m$train$x) - sum(forced), p_max)
  pred <- map(1:p, ~ combn(colnames(m$train$x)[!forced], ., simplify = FALSE)
  ) %>% reduce(c)
  form_list <- map_chr(pred, ~ paste0(.x, collapse = " + "))
  pred_idx <- map(pred, ~ which(colnames(m$train$x) %in% .))
  if(any(forced)){
    pred_idx <- map(pred, ~ c(which(forced), which(colnames(m$train$x) %in% .)))
    pred_idx <- c(list(which(forced)), pred_idx)
    pred_forced <- colnames(m$train$x)[forced]
    if(m$train$intercept) pred_forced <- pred_forced[-1]
    if(length(pred_forced)){
      pred_forced <- paste0(pred_forced, collapse = " + ")
      form_list <- paste(pred_forced, form_list, sep = " + ")
    } else {
      pred_forced <- "1"
    }
    form_list <- c(pred_forced, form_list)
  }
  form_list <- paste("~", form_list)
  #======================================================================
  # Determine number of predictors for each model in list
  #----------------------------------------------------------------------
  n_pred <- pred_idx %>% map_int(~length(.))
  if(m$train$intercept) n_pred <- n_pred - 1L
  #======================================================================
  # Eliminate formulae that exceed p_max
  #----------------------------------------------------------------------
  form_list <- form_list[n_pred <= p_max]
  pred_idx <- pred_idx[n_pred <= p_max]
  n_pred <- n_pred[n_pred <= p_max]
  tibble(pred_idx = pred_idx, form = form_list, n_pred = n_pred)
}

get_subset_stats <- function(j, m, fitter){
  fit_args <- c(
    list(x = m$train$x[, j, drop = FALSE]), m$train[-1])
  fit <- try(suppressMessages(do.call(fitter, fit_args)), silent = TRUE)
  if(inherits(fit, "try-error")) {
    return(NULL)
  }
  fam <- if(grepl("Negative Binomial", fit$family$family))
    "negbin" else fit$family$family
  p <- fit$rank
  if (fam %in% c("gaussian", "Gamma", "inverse.gaussian"))
    p <- p + 1L
  y_hat <- fit$family$linkinv(
    m$train$x[, j, drop = FALSE] %*% fit$coefficients + m$train$offset)
  y_obs <- m$train$y
  if(is.factor(y_obs)) y_obs <- as.integer(y_obs) - 1L
  fit_stats <- c(
    list(aic = fit$aic),
    suppressWarnings(
      predict_metrics_(y_obs, y_hat, family = fam, theta = fit$theta)
    )
  )
  test_stats <- NULL
  test_preds <- NULL
  if(!is.null(m$test)){
    y_hat <- fit$family$linkinv(
      m$test$x[, j, drop = FALSE] %*% fit$coefficients + m$test$offset
    )
    y_obs <- m$test$y
    if(is.factor(y_obs)) y_obs <- as.integer(y_obs) - 1L
    test_stats <- suppressWarnings(
      predict_metrics_(y_obs, y_hat, fam, theta = fit$theta)
    )
    test_preds <- y_hat
  }
  list(fit_stats = fit_stats, test_stats = test_stats, test_preds = test_preds)
}

