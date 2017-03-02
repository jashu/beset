#' beset: Best Subset Predictive Modeling
#'
#' \code{beset} is a portmanteau of BEst subSET, which references the overall
#' objective of this package: to improve prediction accuracy and
#' interpretability of generalized linear models by automatically performing
#' variable selection, either explicitly (e.g., with subset selection) or
#' implicitly (e.g., with coefficient shrinkage), based on cross-validation.
#'
#' It has three main goals:
#' \itemize{
#'   \item Establish a common, user-friendly interface for best subset selection
#'    that works with several model fitting functions, including
#'    \code{\link[stats]{lm}}, \code{\link[stats]{glm}},
#'    \code{\link[MASS]{glm.nb}}, and \code{\link[pscl]{zeroinfl}}.
#'   \item Make it fast and easy to cross-validate generalized linear models,
#'    with a streamlined, parallelized implementation that can take advantage of
#'    multicore computing architectures.
#'   \item Make elastic-net regression more accessible by providing a simple
#'    \code{\link[stats]{lm}}-like wrapper and summary method for
#'    \code{\link[glmnet]{glmnet}}.
#' }
#'
#' @docType package
#' @name beset

NULL
