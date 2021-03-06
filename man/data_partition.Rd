% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_partition.R
\name{data_partition}
\alias{data_partition}
\title{"data_partition" Constructor}
\usage{
data_partition(
  train,
  test,
  y,
  x = NULL,
  offset = NULL,
  weights = NULL,
  na_action = na.omit
)
}
\arguments{
\item{train}{A \code{\link[base]{data.frame}} containing the training data to
be used for model fitting.}

\item{test}{A \code{\link[base]{data.frame}} containing the test data to be
used for model validation.}

\item{y}{\code{Character} string giving the name of the column containing the
response variable to be predicted.}

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
}
\value{
A \code{data_partition} object containing a \code{train} data frame
and a \code{test} data frame.
}
\description{
Constructs an object of class \code{data_partition}.
}
\details{
A \code{data_partition} object is a list containing exactly two data frames
(\code{train} and \code{test}). This object will normally be constructed by
passing a single data frame to \code{\link{partition}}. Use this constructor
function in the event that you wish to manually link two independent data
sets: one to be used for model training and the other to be used for model
testing.

\code{data_partition} objects can be passed as the \code{data} argument to
the \code{beset} modeling functions (\code{\link{beset_lm}},
\code{\link{beset_glm}}, and \code{\link{beset_elnet}}), in which case these
functions will train and cross-validate models using the \code{train} data
and append additional evaluation metrics using the \code{test} data. Note
that in earlier development versions, these functions provided an
optional \code{test_data} argument for this purpose. This has been removed
and you are now required to construct a \code{data_partition} object
beforehand because the \code{data_partition} constructor performs a number
of important checks to insure that your \code{test} data are compatible with
your \code{train} data: 1) all predictor and response variables are present
in both data sets, 2) the levels of all factor variables are the same for
both data sets, 3) if an offset variable is used for model training, an
offset variable is provided for predicting the test data, and 4) unless
\code{na_action} is set to \code{na.pass}, both data frames contain
complete cases with no missing data. The \code{data_partition} constructor
will alert you to potential issues, attempt to resolve them, and return an
error if it can't.
}
\examples{
train <- mtcars[1:16,]
test <- mtcars[17:32,]
factor_names <- c("cyl", "vs", "am", "gear", "carb")
train[factor_names] <- purrr::map_dfc(train[factor_names], factor)
test[factor_names] <- purrr::map_dfc(test[factor_names], factor)
data <- data_partition(train, test, "mpg")

}
