#' Bootstrapped Confidence Intervals
#'
#' Functions for calculating the bias corrected accelerated confidence intervals
#' from either univariate or multivariate data.
#'
#' @seealso \code{\link[boot]{boot}}, \code{\link[boot]{boot.ci}},
#'  \code{\link{cor_loo}}, \code{\link{permute}}
#'
#' @param data The data as a vector, matrix or data frame. If it is a matrix or
#' data frame then each row is considered as one multivariate observation.
#'
#' @param statistic A function which, when applied to data, returns a vector
#' containing the statistic(s) of interest. \code{statistic} whose first two
#' arguments must be as follows: the first argument must be the original data,
#' and the second argument must be a vector of indices, which the function must
#' reference when calculating the statistic. Any additional function arguments
#' must be assigned default values.
#'
#' @param n_rep The number of bootstrap replicates.
#'
#' @param conf A scalar or vector containing the confidence level(s) of the
#' required interval(s).
#'
#' @param ... Other named arguments for \code{statistic} which are passed
#' unchanged each time it is called.  For \code{boot_cor}, this could include
#' a \code{use} and/or \code{method} argument to override the defaults
#' documented in \code{\link{cor_list}}.
#'
#' @return Object of class "bootstrap", a list containing an object of class
#' "boot" and a stats data frame or list with the following vectors:
#' \describe{
#'  \item{statistic}{observed statistic for each measurement}
#'  \item{ci_lower}{lower endpoint of the confidence interval}
#'  \item{ci_upper}{upper endpoint of the confidence interval}
#' }
#' @export

boot_ci <- function(data, statistic, n_rep = 1000, conf = 0.95, ...){
  boot_out <- boot::boot(data = data,
                         statistic = statistic,
                         R = n_rep,
                         ...)
  .get_ci <- function(i){
    boot_ci <-  boot::boot.ci(boot_out, conf = conf, type = "bca", index = i)
    boot_ci$bca[4:5]
  }
  index <- seq_along(boot_out$t0)
  ci <- t(parallel::mcmapply(.get_ci, index))
  stats <- cbind(boot_out$t0, ci)
  colnames(stats) <- c("statistic", "lower", "upper")
  stats <- as.data.frame(stats)
  return(structure(list(boot = boot_out, stats = stats), class = "bootstrap"))
}

#' @describeIn boot_ci Bootstrapped confidence interval for the mean
#' @export
boot_mean <- function(data, n_rep = 1000, conf = 0.95){
  .boot_mean <- function(data, indices){
    apply(data[indices,], 2, mean, na.rm = T)
  }
  out <- boot_ci(data, .boot_mean, n_rep, conf)
  colnames(out$stats)[1] <- "mean"
  return(out)
}

#' @describeIn boot_ci Bootstrapped confidence interval for correlations
#' @export

boot_cor <- function(data, n_rep = 1000, conf = 0.95, ...){
  .boot_cor <- function(data, indices, ...){
      cor(data[indices,], ...)$coef
  }
  out <- boot_ci(data, .boot_cor, n_rep, conf, ...)
  cors <- cor(data, ...)
  colnames(out$stats)[1] <- "coef"
  out$stats <- structure(c(cors, out$stats[,2:3]), class = "cor_list")
  attr(out$stats, "coef") <- attr(cors, "coef")
  return(out)
}

