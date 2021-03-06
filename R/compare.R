#' Compare Predictive Performance of Two Models
#'
#' @param yhat1 A data frame consisting of cross-validated predictions from a
#' benchmark model, an object containing such a data frame, e.g., a
#' "cross_valid" object returned by \code{\link{validate}}, or an object that
#' can be passed to \code{\link{validate}}.
#'
#' @param yhat2 An object of the same type as \code{yhat1} to be compared
#'
#' @param n_rep \code{Integer} giving the number of bootstrap replicates to
#' perform for each repetition of cross-validated predictions. For example,
#' if \code{yhat1} and \code{yhat2} contain 10 columns of predictions, the
#' default value of \code{n_rep = 1000} will result in
#' \eqn{1000 \times 10 = 10,000} replicates total.
#'
#' @param conf Confidence level for the difference between model performance.
#'
#' @inheritParams beset_glm
#'
#' @import parallel
#' @import purrr
#' @export
compare <- function(yhat1, yhat2, n_rep = 1000, conf = 0.95,
                    parallel_type = NULL, n_cores = NULL, cl = NULL,...){
  UseMethod("compare")
}
#' @export
#' @describeIn compare S3 method for class 'cross_valid'
compare.cross_valid <- function(
  yhat1, yhat2, n_rep = 1000, conf = 0.95, parallel_type = NULL, n_cores = NULL,
  cl = NULL, ...
){
  y <- yhat1$parameters$y
  if(!all(y == yhat2$parameters$y)){
    stop("Observed responses for `yhat1` and `yhat2` do not match")
  }
  out <- list(Model1 = map_dbl(yhat1$stats, "mean"),
              Model2 = map_dbl(yhat2$stats, "mean"))
  out$Delta <- out$Model2 - out$Model1
  compare.default(
    yhat1 = yhat1$predictions,
    yhat2 = yhat2$predictions,
    y = y, n_rep = n_rep, conf = conf,
    family = c(yhat1$parameters$family, yhat2$parameters$family),
    theta = c(yhat1$parameters$theta, yhat2$parameters$theta),
    mu = list(yhat1$predictions$mu, yhat2$predictions$mu),
    phi = list(yhat1$predictions$phi,yhat2$predictions$phi),
    parallel_type = parallel_type, n_cores = n_cores, cl = cl,
    out = out
  )
}
compare.default <- function(
  yhat1, yhat2, n_rep = 1000, conf = 0.95, parallel_type = NULL, n_cores = NULL,
  cl = NULL, y, family = "gaussian", mu = NULL, phi = NULL, theta = NULL, ...
){
  extra_args <- list(...)
  out <- extra_args$out
  dim_yhat1 <- dim(yhat1); dim_yhat2 <- dim(yhat2)
  if(!identical(dim_yhat1, dim_yhat2)){
    stop("`yhat1` and `yhat2` must have the same dimensions")
  }
  if(!(NROW(y) == NROW(yhat1) && NROW(y) == NROW(yhat2) &&
       NROW(yhat1) == NROW(yhat2))){
    stop("`y`, `yhat1`, and `yhat2` must have the same number of observations")
  }
  if(length(family) == 1) family <- rep(family, 2)
  # implement hierarchy of families for comparing models fit under different
  # distributions; the logic is to compare predictions under the distribution
  # with more parameters
  if(family[1] != family[2]){
    if(any(family == "zinb" | family == "negbin") & any(family == "zip")){
      family[family == "zip"] <- "zinb"
    } else if(any(family == "zinb" | family == "negbin")){
      family[] <- "negbin"
    }
  }
  if(length(theta) == 1) theta <- rep(theta, 2)
  if(is.null(n_cores) || n_cores > 1){
    parallel_control <- setup_parallel(
      parallel_type = parallel_type, n_cores = n_cores, cl = cl)
    have_mc <- parallel_control$have_mc
    n_cores <- parallel_control$n_cores
    cl <- parallel_control$cl
  }
  if(is.null(out)){
    out <- list(
      Model1 = predict_metrics_(
        y = y, y_hat = yhat1, family = family[1], theta = theta[1],
        mu = mu[[1]], phi = phi[[1]]
      ) %>% as_vector,
      Model2 = predict_metrics_(
        y = y, y_hat = yhat2, family = family[2], theta = theta[2],
        mu = mu[[2]], phi = phi[[2]]
      ) %>% as_vector
    )
    out$Delta <- out$Model2 - out$Model1
  }
  predictive_gain <- if(n_cores > 1L){
    if(have_mc){
      mclapply(1:n_rep, resample_pred_diff, y = y, yhat1 = yhat1, yhat2 = yhat2,
               family = family, theta = theta, mu = mu, phi = phi,
               mc.cores = n_cores)
    } else {
      parLapply(cl, resample_pred_diff, y = y, yhat1 = yhat1, yhat2 = yhat2,
                family = family, theta = theta, mu = mu, phi = phi)
    }
  } else {
    lapply(1:n_rep, resample_pred_diff, y = y, yhat1 = yhat1, yhat2 = yhat2,
           family = family, theta = theta, mu = mu, phi = phi)
  }
  predictive_gain <- transpose(predictive_gain) %>% simplify_all
  a <- (1 - conf)/2
  a <- c(a, 1 - a)
  ci_gain <- map(predictive_gain, ~ quantile(.x, probs = a, na.rm = TRUE))
  out$`95% CI` <- ci_gain
  out$predictive_gain <- predictive_gain
  names(out)[4] <- paste(format(conf * 100, digits = 2), "%", " CI", sep = "")
  structure(out, class = "predictive_gain", family = family[1])
}

resample_pred_diff <- function(
  seed, y, yhat1, yhat2, family, theta = NULL, mu = NULL, phi = NULL
){
  set.seed(seed)
  i <- sample(seq_along(y), replace = TRUE)
  # needs to be modified to resample mu and phi also, if they are present
  pred1 <- if(is.vector(yhat1)){
    if(grepl("^zi", family[1])){
      predict_metrics_(
        y[i], yhat1[i], family[1], theta[1], mu[[1]][i], phi[[1]][i]
      )
    } else predict_metrics_(y[i], yhat1[i], family[1], theta[1])
  } else {
    if(grepl("^zi", family[1])){
      pmap(list(yhat1, mu[[1]], phi[[1]]), function(yhat, m, p){
        predict_metrics_(y[i], yhat[i], family[1], theta[1], m[i], p[i])
      }) %>% transpose %>% simplify_all
    } else {
      map(yhat1, ~ predict_metrics_(y[i], .x[i], family[1], theta[1])) %>%
        transpose %>% simplify_all
    }
  }
  pred2 <- if(is.vector(yhat2)){
    if(grepl("^zi", family[2])){
      predict_metrics_(
        y[i], yhat2[i], family[2], theta[2], mu[[2]][i], phi[[2]][i]
      )
    } else predict_metrics_(y[i], yhat2[i], family[2], theta[2])
  } else {
    if(grepl("^zi", family[2])){
      pmap(list(yhat2, mu[[2]], phi[[2]]), function(yhat, m, p){
        predict_metrics_(y[i], yhat[i], family[2], theta[2], m[i], p[i])
      }) %>% transpose %>% simplify_all
    } else {
      map(yhat2, ~ predict_metrics_(y[i], .x[i], family[2], theta[2])) %>%
        transpose %>% simplify_all
    }
  }
  map2(pred1, pred2, ~ .y - .x)
}

#' @export
#' @describeIn compare S3 method for class 'numeric'
compare.numeric <- function(
  yhat1, yhat2, n_rep = 1000, conf = 0.95, parallel_type = NULL, n_cores = NULL,
  cl = NULL, y, ...){
  compare.default(
    yhat1, yhat2, n_rep, conf, parallel_type, n_cores, cl, y, ...
  )
}

#' @export
#' @describeIn compare two 'beset' class models
compare.beset <- function(
  yhat1, yhat2, n_rep = 1000, conf = 0.95, parallel_type = NULL, n_cores = NULL,
  cl = NULL, ...
){
  compare.cross_valid(
    validate(yhat1, ...), validate(yhat2, ...), n_rep, conf, parallel_type,
    n_cores, cl, ...
  )
}

#' @export
#' @describeIn compare two 'glm' or 'lm' class models
compare.lm <- function(
  yhat1, yhat2, n_rep = 1000, conf = 0.95, parallel_type = NULL, n_cores = NULL,
  cl = NULL, ...
){
  compare.cross_valid(
    validate(yhat1, ...), validate(yhat2, ...), n_rep, conf, parallel_type,
    n_cores, cl, ...
  )
}

#' @export
#' @describeIn compare two 'zeroinfl' class models
compare.zeroinfl <- function(yhat1, yhat2, ...){
  compare.cross_valid(
    validate(yhat1, ...), validate(yhat2, ...), n_rep, conf, parallel_type,
    n_cores, cl, ...
  )
}
