#' Robust Significance Testing for Between-Group Differences in Correlations
#'
#' Given the same set of variables from two independent samples, \code{cor_dif}
#' calculates bivariate correlations for each sample and tests whether each
#' correlation is signficantly different between samples.
#'
#' \code{cor_dif} was designed to test multiple hypotheses of group differences
#' in correlations between elements of a network, such as whether a treatment
#' alters the correlation between the activity of two brain regions. Fisher's z
#' test is traditionally used for this purpose, but may yield inaccurate
#' probability estiamtes when a) one or both variables being correlated is not
#' normally distributed, b) the sample size for each group is small (n <= 10),
#' or c) there are multiple pairwise correlations from the same data set.
#' Although the latter is a multiple-comparisons problem, conventional
#' statistical corrections assume that the hypothesis tests are independent,
#' which is assuredly not the case for multiple hypothesis tests on correlations
#' arising from a networked system.
#'
#' \code{cor_dif} can be used with small sample sizes and, by default, uses
#' leave-one-out cross validation to produce more realistic estimates of the
#' correlation coefficients themselves, provides bootstrapped estimates of the
#' 95\% confidence intervals for the correlation coefficients, and conducts
#' permutation tests to calculate p-values conditioned on the number of
#' correlations being assessed. (For example, if 15 correlation coefficients are
#' being computed, then the resulting p value corresponds to the probability of
#' calulcating 15 correlation differences and obtaining at least one as extreme
#' as the observed difference.)
#'
#' @param form A formula of the form \code{groups ~ x1 + x2 + ...}. That is, the
#'  left-hand side indicates the grouping factor and the right-hand side
#'  indicates the (numeric) measurements to correlate within each group.
#'
#' @param data Data frame from which variables specified in \code{form} are to
#'  be taken.
#'
#' @param n_iter Number of iterations to use for the bootstrap and permutation
#' tests.
#'
#' @param loo Logical indicating whether permutation tests for group differences
#' should be performed on the original correlation coefficients (\code{FALSE})
#' or on the leave-one-out coefficients (\code{TRUE}).
#'
#' @seealso \code{\link{cor_loo}}, \code{\link{boot_cor}},
#' \code{\link{permute_cor_dif}}
#'
#' @import dplyr
#' @export
#'
cor_dif <- function(form, data, n_iter = 1000, loo = TRUE){
  mf <- model.frame(form, data, na.action = na.pass)
  x <- unlist(mf[,1])
  if(!is.factor(x)) x <- factor(x)
  y <- mf[,2:ncol(mf)]
  group <- levels(x)
  if(length(group) > 2){
    warning(paste("You have specified more than 2 groups.",
                  "I will only compute correlation differences",
                  "for the first 2 groups.", "If this is not what you want,",
                  "abort now or forever hold your peace."), immediate. = TRUE)
  }
  numeric_vars <- sapply(y, is.numeric)
  if(!all(numeric_vars)){
    warning(paste("The following variables are not numeric",
                  "and will be excluded from the analysis:",
                  names(numeric_vars)[!numeric_vars]),
            immediate. = TRUE)
    y <- y[, numeric_vars]
  }
  group1 <- y[x == group[1], ]
  group2 <- y[x == group[2], ]
  group1_cor <- boot_cor(group1, n_iter)
  group2_cor <- boot_cor(group2, n_iter)
  group1_cor$stats$loo_r <- cor_loo(group1)
  group2_cor$stats$loo_r <- cor_loo(group2)
  data <- cbind(x, y); names(data)[1] <- names(mf)[1]
  group_dif <- permute_cor_dif(form, data, loo, n_iter)
  return(structure(
    list(group1 = group1_cor, group2 = group2_cor,
         diff = group_dif),
    class = "cor_dif"))
}
