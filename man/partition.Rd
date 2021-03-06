% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/partition.R
\name{partition}
\alias{partition}
\title{Partition Data into Training and Test Model Frames}
\usage{
partition(
  data,
  y,
  frac = 0.5,
  x = NULL,
  offset = NULL,
  weights = NULL,
  na_action = na.omit,
  seed = 42
)
}
\arguments{
\item{data}{Data frame to be partitioned.}

\item{y}{\code{Character} string giving the name of the column containing the
response variable to be predicted.}

\item{frac}{The fraction of data that should be included in the training set.
Default is \code{0.5}.}

\item{x}{(Optional) \code{character} vector giving the names of the columns
containing the predictor variables. If omitted, defaults to all columns
other than those named as \code{y}, \code{offset}, or \code{weights}.}

\item{offset}{(Optional) \code{character} string giving the name of the
column containing a model offset. An offset is a known predictor that is
added to a linear model \emph{as is} (with a beta coefficient of 1) rather
than having its beta coefficient optimized. If given, an offset must be
included for both the \code{train} and \code{test} data frames.}

\item{weights}{(Optional) \code{character} string giving the name of the
column containing observation weights. Use these if you want some rows of the
data frame to exert more or less influence than others on a model fit. If
given, the \code{weights} column is only applied during model training;
a \code{weights} column in the \code{test} data will be ignored.}

\item{na_action}{\code{Function} defining how \code{NA}s shoud be treated.
Options include \code{\link[stats]{na.omit}} (default),
\code{\link[stats]{na.fail}}, \code{\link[stats]{na.exclude}}, and
\code{\link[stats]{na.pass}}.}

\item{seed}{\code{Integer} value for seeding the random number generator. See
\code{\link[base]{set.seed}}.}
}
\value{
An object of class "data_partition": a list containing two model
frames named \code{train} and \code{test}, containing the training and
testing sets, respectively.
}
\description{
\code{partition} randomly splits a data frame into two model frames,
\code{train} and \code{test}, which are returned as a
"data_partition" structure.
}
\details{
\code{partition} creates a train/test split among the rows of a data frame
based on stratified random sampling within the factor levels of a
classification outcome or the quartiles of a numeric outcome. This insures
that the training and test samples will be closely matched in terms of class
incidence or frequency distribution of the outcome measure. \code{partition}
includes a \code{seed} argument so that the randomized partitioning is
reproducible. The \code{train} and \code{test} data frames are returned
bound together in a \code{\link{data_partition}} structure so that their
common ancestry is maintained and self-documented. For example, if you name
your \code{\link{data_partition}} "\code{data}", you can intuitively access
the training set with \code{data$train} and its corresponding test set with
\code{data$test}.
}
\examples{
data <- mtcars
factor_names <- c("cyl", "vs", "am", "gear", "carb")
data[factor_names] <- purrr::map_dfc(data[factor_names], factor)
data <- partition(data, y = "mpg")

}
\seealso{
\code{\link[base]{set.seed}}, \code{\link{data_partition}}
}
