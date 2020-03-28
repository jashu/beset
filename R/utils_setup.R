#' @import purrr
#' @import dplyr

#============================================================================
# Create covariate matrices and response vectors for use with beset functions
#----------------------------------------------------------------------------
make_args <- function(form, data, family, link, contrasts, force_in = NULL,
                      weights = NULL, offset = NULL,
                      start = NULL, etastart = NULL, mustart = NULL,
                      epsilon, maxit, check = TRUE, ...){
  if(inherits(data, "data_partition")){
    return(
      list(
        train = make_args(
          form, data$train, family, link, contrasts, force_in, weights, offset,
          start, etastart, mustart, epsilon, maxit, check, ...
          )$train,
        test = make_args(
          form, data$test, family, link, contrasts, force_in, weights, offset,
          start, etastart, mustart, epsilon, maxit, check = FALSE, ...
          )$train
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
    model.frame(form, data, na.action = na.omit)
  } else {
    offset <- model.offset(data)
    weights <- model.weights(data)
    data$`(offset)` <- data$`(weights)` <- NULL
    data
  }
  x <- stats::model.matrix(form, mf, contrasts)
  intercept <- "(Intercept)" %in% colnames(x)
  frc_cols <- forced_cols(mf, x, force_in)
  y <- model.response(mf)
  # insure y is a factor if family is binomial
  if(family == "binomial" && !is.factor(y)) y <- factor(y)
  family <- if(family == "negbin"){
    stats::poisson(link = link)
  } else {
    do.call(family, list(link = link))
  }
  if(is.null(offset)) offset <- rep(0, nrow(data))
  if(is.null(weights)) weights <- rep(1, nrow(data))
  # arguments passed to glmnet
  if("alpha" %in% names(list(...))){
    if(intercept){
      x <- x[, -1]
      if(!is.null(force_in)){
        frc_cols <- frc_cols[-1] - 1L
      }
    }
    penalty <- rep(1, ncol(x))
    penalty[frc_cols] <- 0
    list(
      train = c(
        list(x = x, y = y, family = family$family, weights = weights,
             offset = offset, intercept = intercept, penalty.factor = penalty,
             thresh = epsilon, maxit = maxit),
        list(...)
      ),
      test = NULL
    )
  # arguments passed to glm
  } else {
    list(
      train = structure(
        list(
          x = x, y = y, weights = weights, start = start, etastart = etastart,
          mustart = mustart, offset = offset, family = family,
          control = glm.control(epsilon = epsilon, maxit = maxit),
          intercept = intercept
        ),
        forced = frc_cols
      ),
      test = NULL
    )
  }
}

forced_cols <- function(mf, x, force_in = NULL){
  factors <- attr(terms(mf), "dataClasses") == "factor"
  factors <- names(factors[factors])
  factors <- factors[factors %in% force_in]
  intercept <- "(Intercept)" %in% colnames(x)
  if(length(factors) > 0){
    dummies <- map(
      factors, ~ {
        x <- paste(.x, .getXlevels(terms(mf), mf)[[.x]], sep = "")
        if(intercept) x[-1] else x
      }
    ) %>% as_vector
    force_in <- c(setdiff(force_in, factors), dummies)
  }
  if(intercept){
    force_in <- c("(Intercept)", force_in)
  }
  match(force_in, colnames(x))
}

