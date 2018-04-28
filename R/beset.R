#' beset: Best Subset Predictive Modeling
#'
#' \code{beset} is a portmanteau of BEst subSET, which references the overall
#' objective of this package: to improve the generalizability of data-driven
#' model selection through cross-validation. To learn more about \code{beset},
#' start with the vignettes: \code{browseVignettes(package = "beset")}
#'
#' @section Overarching goals:
#' \enumerate{
#'   \item Provide a fast and easy way to cross-validate GLMs.
#'   \item Establish a common, user-friendly interface for best subset selection
#'     that works with several model fitting functions (\code{\link[stats]{lm}},
#'     \code{\link[stats]{glm}}, and \code{\link[MASS]{glm.nb}}.
#'   \item Make elastic-net regression more accessible and interpretable by
#'     providing a wrapper to \code{\link[glmnet]{glmnet}} that maintains the
#'     same user interface as \code{\link[stats]{lm}} and
#'     \code{\link[stats]{glm}} provides informative summary and plot methods.
#' }
#'
#' @section Overview of principal functions:
#' \describe{
#'  \item{\code{\link{partition}}}{Makes it quick and easy to split a data
#'    frame into a train/test partition and keep track of which train set goes
#'    with which test set.}
#'  \item{\code{\link{beset_glm}}}{Performs best subset selection using repeated
#'    cross-validation to find the optimal number of predictors for several
#'    families of generalized linear models.}
#'  \item{\code{\link{beset_elnet}}}{Enhances elastic-net regression with
#'    \code{\link[glmnet]{glmnet}} by 1) allowing the user to specify a model
#'    using R's formula syntax, 2) allowing the user to simultaneously tune both
#'    alpha and lambda using cross-validation, and 3) providing a summary output
#'    that ranks the relative importance of the predictors that survived
#'    shrinkage and gives some statistics indicating how well the model fits the
#'    training data and predicts new data.}
#'  \item{\code{\link{validate}}}{An S3 object system that makes it very easy
#'   to obtain cross-validated prediction stats for a previously fit model.}
#'  }
#' @docType package
#' @name beset

NULL
