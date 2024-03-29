% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/importance.R
\name{importance}
\alias{importance}
\alias{importance.nested}
\alias{importance.beset_rf}
\alias{importance.part_depend}
\alias{plot.variable_importance}
\title{Relative Variable Importance}
\usage{
importance(object, ...)

\method{importance}{nested}(object, ...)

\method{importance}{beset_rf}(object, ...)

\method{importance}{part_depend}(object, ...)

\method{plot}{variable_importance}(x, p_max = 20, labels = NULL, max_import = NULL, ...)
}
\arguments{
\item{object}{A model object for which variable importance scores are
desired.}

\item{...}{Additional named arguments that define the model selection rules.
See \code{\link{summary.beset}}.}

\item{p_max}{Maximum number of predictors to include in variable importance
plot}

\item{labels}{(Optional) two-column \code{data.frame} where column 1
gives the variable names used in the model and column 2 gives a corresponding
descriptive label. If \code{labels} are defined, the variable importance plot
will replace the model variable names with their descriptive labels.}

\item{max_import}{(Optional) numeric value specifying the maximum value of
the x axis for the relative importance plot. Default is the upper limit of
the largest importance score}
}
\description{
\code{importance} is a generic function for obtaining resampled variable
importance scores. The function invokes particular
\code{\link[utils]{methods}} which depend on the \code{\link[base]{class}}
of the first argument.
}
\section{Methods (by class)}{
\itemize{
\item \code{importance(nested)}: Determine variable importance for "nested"
\code{\link{beset_glm}} and \code{\link{beset_elnet}} objects

\item \code{importance(beset_rf)}: Relative importance method for "beset_rf" objects

\item \code{importance(part_depend)}: Relative importance based on partial dependence

}}
\examples{
rf <- beset_rf(Fertility ~ ., data = swiss)
importance(rf)

}
