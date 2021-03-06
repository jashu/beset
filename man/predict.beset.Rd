% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods_predict.R
\name{predict.beset}
\alias{predict.beset}
\title{Predict Methods for \code{beset} Objects}
\usage{
\method{predict}{beset}(
  object,
  newdata,
  type = "response",
  newoffset = NULL,
  alpha = NULL,
  lambda = NULL,
  n_pred = NULL,
  metric = "auto",
  oneSE = TRUE,
  na.action = na.pass,
  tt = NULL,
  ...
)
}
\arguments{
\item{object}{A \code{beset_elnet} object.}

\item{newdata}{An optional data frame in which to look for variables with
    which to predict.  If omitted, the fitted values are used.}

\item{type}{Type of prediction required. Type "link" gives the linear
predictors for "binomial" and "poisson" models; for "gaussian" models it
gives the fitted values. Type "response" gives the fitted probabilities for
"binomial", fitted mean for "poisson", and fitted values for  "gaussian".
Type "coefficients" computes the coefficients. Note that for "binomial"
models, results are returned only for the class corresponding to the second
level of the factor response. Type "class" applies only to "binomial" models,
and produces the class label corresponding to the maximum probability. Type
"nonzero" returns a list of the indices of the nonzero coefficients.}

\item{newoffset}{If an offset is used in the fit, then one must be supplied
for making predictions (except for \code{type="coefficients"} or
\code{type="nonzero"})}

\item{alpha}{(Optional) \code{numeric} value to force selection of
elastic-net model with the given \code{alpha} parameter. If left \code{NULL},
the best value of \code{alpha} will be chosen using the cross-validation
\code{metric} and \code{oneSE} rule.}

\item{lambda}{(Optional) \code{numeric} value to force selection of elastic-
net model with the given \code{lambda} parameter. If left \code{NULL},
the best value of \code{lambda} will be chosen using the cross-validation
\code{metric} and \code{oneSE} rule.}

\item{n_pred}{(Optional) \code{integer} number of predictors that the best
model should contain. If specified, all other arguments are ignored.}

\item{metric}{\code{Character} string giving prediction metric on which to
base model selection. Can be one of \code{"mce"} (mean cross entropy--the
default), \code{"mse"} (mean squared error), \code{"aic"} (Akaike informatio
criterion, not applicable for elastic net), \code{"auc"} (area under the ROC
curve--only applicable if response is binomial), or \code{"mae"} (mean
absolute error--only applicable if response is numeric).}

\item{oneSE}{\code{Logical} indicating whether or not to use the "one
standard error" rule. If \code{TRUE} (default) the simplest model within one
standard error of the optimal model is returned. If \code{FALSE} the model
with the optimal cross-validation performance is returned.}

\item{na.action}{function determining what should be done with missing
    values in \code{newdata}.  The default is to predict \code{NA}.}

\item{...}{further arguments passed to or from other methods.}
}
\description{
Predict Methods for \code{beset} Objects
}
