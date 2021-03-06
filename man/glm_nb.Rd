% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/glm_nb.R
\name{glm_nb}
\alias{glm_nb}
\title{Fit a Negative Binomial Generalized Linear Model}
\usage{
glm_nb(
  x,
  y,
  weights = rep(1, nobs),
  start = NULL,
  etastart = NULL,
  mustart = NULL,
  offset = rep(0, nobs),
  family = poisson(),
  control = list(),
  intercept = TRUE
)
}
\arguments{
\item{x}{A design matrix of dimension \code{n * p}}

\item{y}{A vector of observations of length \code{n}.}

\item{weights}{(Optional) \code{numeric vector} of prior weights placed on
the observations during model fitting. Default is \code{NULL}.}

\item{start}{(Optional) starting values for the parameters in the linear
predictor.}

\item{etastart}{(Optional) starting values for the linear predictor.}

\item{mustart}{(Optional) starting values for the vector of means.}

\item{offset}{A vector of length \code{nobs} specifying an \emph{a priori}
known component that will be added to the linear predictor
before applying the link function. Default is NULL.}

\item{family}{The result of a call to either \code{\link[stats]{poisson}} or
\code{\link[MASS]{negative.binomial}}.}

\item{control}{A list of parameters for controlling the fitting process to be
passed to \code{\link[stats]{glm.control}}}
}
\description{
A modification of \code{\link[MASS]{glm.nb}} to provide a more efficient
workhorse function analagous to \code{\link[stats]{glm.fit}} where the
response vector, design matrix, and family have already been calculated.
}
