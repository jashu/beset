
#======================================================================
# Make list of all possible formulas with number of predictors <= p_max
#----------------------------------------------------------------------

get_subsets <- function(m, p_max){
  forced <- attr(m$train, "forced")
  p <- min(ncol(m$train$x) - length(forced), p_max)
  pred <- map(
    1:p, ~ combn(colnames(m$train$x)[-forced], ., simplify = FALSE)
  ) %>% reduce(c)
  form_list <- map_chr(pred, ~ paste0(.x, collapse = " + "))
  pred_idx <- map(pred, ~ which(colnames(m$train$x) %in% .))
  if(length(forced)){
    pred_idx <- map(pred, ~ c(forced, which(colnames(m$train$x) %in% .)))
    pred_idx <- c(list(forced), pred_idx)
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

