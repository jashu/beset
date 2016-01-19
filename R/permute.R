#' Permutation Test for Multiple Comparisons
#'
#' Computes the empirical distribution of the maximum of a test statistic
#' repeated across multiple measurements under the null hypothesis.
#'
#' @section Warning:
#' If you are writing your own \code{statistic} function, be sure that either
#' your variables are centered and scaled or, better yet, that your
#' \code{statistic} is itself standardized.
#'
#' @seealso \code{\link{boot_ci}}, \code{\link{cor_dif}}
#'
#' @param form A formula of the form \code{group ~ var1 + var2 + ...}. That is,
#'  the left-hand side specifies the grouping factor, and the right-hand side
#'  specifies the (numeric) variables that are to be evaluated for group
#'  differences.
#'
#' @param data Data frame from which variables specified in \code{form} are to
#'  be taken.
#'
#' @param statistic A function which when applied to data returns a vector
#' containing the effect-size statistic(s) of interest. \code{statistic} must
#' take exactly two arguments. The first argument passed will be a vector of
#' group labels. The second will be a data frame of measurements whose rows
#' correspond to the group labels. \strong{If the statistic is signed (i.e., it
#' reflects the direction of the group difference) and you desire a two-tailed
#' test, be sure that \code{statistic} returns the \emph{absolute value} of the
#' effect-size statistic.}
#'
#' @param n_perm The number of permutations to perform.
#'
#' @param loo For \code{permute_cor_dif}, logical indicating whether to use
#' leave-one-out estimates of the correlation coefficients (\code{TRUE} by
#' default) or the original coefficients (\code{FALSE}).
#'
#' @return Object of class "permutation", a list containing the following items:
#' \describe{
#'  \item{sampling_distribution}{Numeric vector containing the maximum group
#'  difference in terms of absolute value (largest group difference regardless
#'  of direction) obtained for each random permutation of the group labels. This
#'  corresponds to the distribution of the maximum effect sizes that would be
#'  obtained from your multivariate data set under the null hypothesis.}
#' \item{effect_size}{Numeric vector containing the effect size for each of the
#'  dependent variables. For \code{permute_mean_dif} and \code{permute_cor_dif}
#'  this will be Cohen's \eqn{d}}
#' \item{p_value}{Numeric vector containing the empirical probability of
#' obtaining the observed test statistic for each of the dependent variables,
#' given that the null hypothesis is true and conditioned on the number of
#' variables assessed and their covariance.}
#' }
#'
#' @export

permute <- function(form, data, statistic, n_perm = 1000){
  mf <- model.frame(form, data, na.action = na.pass)
  x <- unlist(mf[,1])
  y <- mf[,2:ncol(mf)]
  observed <- as.array(statistic(x, y))
  .permute <- function(seed){
    set.seed(seed)
    x <- sample(x)
    permutation <- statistic(x, y)
    max(permutation, na.rm = T)
  }
  seed <- 1:n_perm
  distribution <- parallel::mcmapply(.permute, seed)
  get_p <- function(observed, distribution){
    (sum(distribution >= observed) + 1) / (n_perm + 1)
  }
  p <- apply(observed, 1, get_p, distribution)
  return(structure(
    list(sampling_distribution = distribution,
         effect_size = as.vector(observed),
         p_value = p),
    class = "permutation"))
}

#' @describeIn permute Permutation test for mean differences between groups.
#' @export

permute_mean_dif <- function(form, data, n_perm = 1000){
  .mean_dif <- function(x, y){
    group <- unique(x)
    y <- apply(y, 2, scale)
    mean1 <- apply(y[x == group[1],], 2, mean, na.rm = T)
    mean2 <- apply(y[x == group[2],], 2, mean, na.rm = T)
    d <- abs(mean1 - mean2)
    return(d)
  }
  permute(form, data, .mean_dif, n_perm)
}

#' @describeIn permute Permutation test for correlation differences between
#' groups. \code{permute_cor_dif} performs a permutation test on the
#' between-group z-score difference for all bivariate corrleations obtained
#' from two independent samples.
#' @export

permute_cor_dif <- function(form, data, loo = TRUE, n_perm = 1000){
  if(loo){
    .cor_dif <- function(x, y){
      group <- unique(x)
      cor1 <- cor_loo(y[x == group[1],])
      cor2 <- cor_loo(y[x == group[2],])
      d <- abs(atanh(cor1) - atanh(cor2))
      return(d)
    }
  } else {
    .cor_dif <- function(x, y){
      group <- unique(x)
      cor1 <- cor_list(y[x == group[1],])
      cor2 <- cor_list(y[x == group[2],])
      d <- abs(atanh(cor1) - atanh(cor2))
      return(d)
    }
  }
  permute(form, data, .cor_dif, n_perm)
}
