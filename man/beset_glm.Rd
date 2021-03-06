% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/beset_glm.R
\name{beset_glm}
\alias{beset_glm}
\alias{beset_lm}
\title{Best Subset Selection for Generalized Linear Models}
\usage{
beset_glm(
  form,
  data,
  family = "gaussian",
  link = NULL,
  p_max = 10,
  force_in = NULL,
  nest_cv = FALSE,
  n_folds = 10,
  n_reps = 10,
  seed = 42,
  contrasts = NULL,
  offset = NULL,
  weights = NULL,
  start = NULL,
  etastart = NULL,
  mustart = NULL,
  epsilon = 1e-08,
  maxit = 25,
  skinny = FALSE,
  n_cores = NULL,
  parallel_type = NULL,
  cl = NULL
)

beset_lm(
  form,
  data,
  p_max = 10,
  force_in = NULL,
  weights = NULL,
  contrasts = NULL,
  offset = NULL,
  nest_cv = FALSE,
  n_folds = 10,
  n_reps = 10,
  seed = 42,
  n_cores = NULL,
  parallel_type = NULL,
  cl = NULL
)
}
\arguments{
\item{form}{A model \code{\link[stats]{formula}}.}

\item{data}{Either a \code{\link{data_partition}} object containing data sets
to be used for both model training and testing, or a single data frame that
will be used for model training and cross-validation.}

\item{family}{Character string naming the error distribution to be used in
the model. Available families are listed under 'List of available families
and link functions'.}

\item{link}{(Optional) character string naming the link function to be used in
the model. Available links and their defaults differ by \code{family} and are
listed under 'List of available families and link functions'.}

\item{p_max}{Maximum number of predictors to attempt to fit. Default is 10.}

\item{force_in}{(Optional) character vector containing the names of any
predictor variables that should be included in every model. (Note that if
there is an intercept, it is forced into every model by default.)}

\item{nest_cv}{\code{Logical} value indicating whether to perform nested
cross-validation. If \code{nest_cv = TRUE}, the cross-validation used to
select the best model is nested within a cross-validation used to estimate
prediction error on a new sample, thus providing as estimate of test error
that is free from potential selection bias. Because this multiplicatively
increases compute times by a factor equal to the number of folds, the
default is \code{FALSE}. Note that setting this parameter to \code{TRUE} will
provide more informative summary output regarding the uncertatinty in the
selection procedure itself, i.e., how often a given model is chosen as
"best" according to the given criteria, and is necessary in order for the
returned objects to work with certain \code{beset} methods, such as
\code{\link{compare}} and \code{\link{importance}}.}

\item{n_folds}{\code{Integer} indicating the number of folds to use for
cross-validation.}

\item{n_reps}{\code{Integer} indicating the number of times cross-validation should
be repeated (with different randomized fold assignments).}

\item{seed}{\code{Integer} used to seed the random number generator when
assigning observations to folds.}

\item{contrasts}{(Optional) \code{list}. See the \code{contrasts.arg} of
\code{\link[stats]{model.matrix.default}}.}

\item{offset}{(Optional) \code{vector} of length \code{nobs} specifying an
\emph{a priori} known component that will be added to the
linear predictor before applying the link function. Useful for
the "\code{poisson}" family (e.g. log of exposure time), or for
refining a model by starting at a current fit. Default is
\code{NULL}.}

\item{weights}{(Optional) \code{numeric vector} of prior weights placed on
the observations during model fitting. Default is \code{NULL}.}

\item{start}{(Optional) starting values for the parameters in the linear
predictor.}

\item{etastart}{(Optional) starting values for the linear predictor.}

\item{mustart}{(Optional) starting values for the vector of means.}

\item{epsilon}{\code{Numeric} value of positive convergence tolerance ε; the
iterations converge when \eqn{|dev - dev_{old}|/(|dev| + 0.1) < ε}. Default
is \code{1e-8}.}

\item{maxit}{\code{Integer} giving the maximal number of IWLS iterations.
Default is 25.}

\item{skinny}{\code{Logical} value indicating whether or not to return a
"skinny" model. If \code{FALSE} (the default), the return object will include
a copy of the model \code{\link[stats]{terms}}, \code{data},
\code{contrasts}, and a record of the \code{xlevels} of the factors used in
fitting. If these features are not needed, setting \code{skinny = TRUE} will
prevent these copies from being made.}

\item{n_cores}{Integer value indicating the number of workers to run in
parallel during subset search and cross-validation. By default, this will
be set to one fewer than the maximum number of physical cores you have
available, as indicated by \code{\link[parallel]{detectCores}}. Set to 1 to
disable parallel processing.}

\item{parallel_type}{(Optional) character string indicating the type of
parallel operation to be used, either \code{"fork"} or \code{"sock"}. If
omitted and \code{n_cores > 1}, the default is \code{"sock"} for Windows and
otherwise either \code{"fork"} or \code{"sock"} depending on which process is
being run.}

\item{cl}{(Optional) \code{\link[parallel]{parallel}} or
\code{\link[snow]{snow}} cluster for use if \code{parallel_type = "sock"}.
If not supplied, a cluster on the local machine is automatically created.}
}
\value{
A "beset_glm" object with the following components:
\describe{
   \item{stats}{a list with three data frames:
     \describe{
       \item{fit}{statistics for every possible combination of predictors:
         \describe{
           \item{n_pred}{the total number of predictors in model; note that
              the number of predictors for a factor variable corresponds to the
              number of factor levels minus 1}
           \item{form}{formula for model}
           \item{aic}{\eqn{-2*log-likelihood + k*npar}, where \eqn{npar}
             represents the number of parameters in the fitted model, and
             \eqn{k = 2}}
           \item{dev}{twice the difference between the log-likelihoods of the
             saturated and fitted models, multiplied by the scale parameter}
           \item{mae}{mean absolute error}
           \item{mce}{mean cross entropy, estimated as
             \eqn{-log-likelihood/N}, where \eqn{N} is the number of
             observations}
           \item{mse}{mean squared error}
           \item{r2}{R-squared, calculated as
             \eqn{1 - deviance/null deviance}}
           }
         }
     \item{cv}{a data frame containing cross-validation statistics
     for the best model for each \code{n_pred} listed in \code{fit_stats}.
     Each metric is computed using \code{\link{predict_metrics}}, with
     models fit to \eqn{n-1} folds and predictions made on the left-out fold.
     Each metric is followed by its standard error. The data frame
     is otherwise the same as that documented for \code{fit}, except
     AIC is omitted.}
     \item{test}{if \code{test_data} is provided, a data frame
     containing prediction metrics for the best model for each \code{n_pred}
     listed in \code{fit} as applied to the \code{test_data}.}
     }
   }
  \item{fold_assignments}{list giving the row indices for the holdout
   observations for each fold and/or repetition of cross-validation}
  \item{n_folds}{number of folds used in cross-validation}
  \item{n_reps}{number of repetitions used in cross-validation}
  \item{family}{name of error distribution used in the model}
  \item{link}{name of link function used in the model}
  \item{terms}{the \code{\link[stats]{terms}} object used}
  \item{data}{the \code{data} argument}
  \item{offset}{the offset vector used}
  \item{contrasts}{(where relevant) the contrasts used}
  \item{xlevels}{(where relevant) a record of the levels of the factors used
       in fitting}
    }
}
\description{
\code{beset_glm} performs best subset selection using repeated
cross-validation to find the optimal number of predictors for several
families of generalized linear models.
}
\details{
\code{beset_glm} performs best subset selection for generalized linear
models, fitting a separate model for each possible combination of predictors
(all models that contain exactly 1 predictor, all models that contain
exactly 2 predictors, and so forth). For each number of predictors,
\code{beset_glm} first picks the model with the best fit and then
estimates how well this model predicts new data using \code{k}-fold
cross-validation (how well, on average, a model trained using \eqn{k - 1}
folds predicts the left-out fold).
}
\section{Cross-validation details}{

\code{beset_glm} randomly partitions the data set into \code{n_folds} *
\code{n_repeats} folds within strata (factor levels for factor outcomes,
percentile-based groups for numeric outcomes). This insures that the folds
will be matched in terms of the outcome's frequency distribution.
\code{beset_glm} also insures the reproducibility of your analysis by
requiring a \code{seed} to the random number generator as one of its
arguments.
}

\section{List of available families and link functions}{

\describe{
 \item{\code{"gaussian"}}{The Gaussian family accepts the links
   \code{"identity"} (default), \code{"log"}, and \code{"inverse"}.}
 \item{\code{"binomial"}}{The binomial family accepts the links
   \code{"logit"} (default), \code{"probit"}, \code{"cauchit"}, \code{"log"}, and
   \code{"cloglog"} (complementary log-log).}
 \item{\code{"poisson"}}{The Poisson family accepts the links \code{"log"}
   (default), \code{"sqrt"}, and \code{"identity"}.}
 \item{\code{"negbin"}}{The negative binomial family accepts the links
   \code{"log"} (default), \code{"sqrt"}, and \code{"identity"}.}
 }
}

\section{Warnings}{

\enumerate{
 \item \code{beset_glm} handles missing data by performing listwise deletion.
  No other options for handling missing data are provided at this time. The
  user is encouraged to deal with missing values prior to running this
  function.
 \item \code{beset_glm} is intended for use with additive models only.  An
 exhaustive search over the space of possible interactions and/or non-linear
 effects is computationally prohibitive, but I hope to offer a greedy search
 option in the future. In the meantime and in general, I would recommend the
   \href{https://cran.r-project.org/web/packages/earth/index.html}{MARS}
   technique if you require this feature.
 \item \code{beset_glm} is best suited for searching over a small number of
 predictors (less than 10). For a large number of predictors (more than 20),
 \code{\link{beset_elnet}} is recommended instead. However, note that
 \code{\link{beset_elnet}} only works with a more restricted set of
 distributions.
}
}

\examples{
subset1 <- beset_glm(Fertility ~ ., data = swiss)
summary(subset1)

# Force variables to be included in model
subset2 <- beset_glm(Fertility ~ ., data = swiss,
                     force_in = c("Agriculture", "Examination"))
summary(subset2)

# Use nested cross-validation to evaluate error in selection
subset3 <- beset_glm(Fertility ~ ., data = swiss, nest_cv = TRUE)
summary(subset3)
}
\seealso{
\code{\link[stats]{glm}},
\code{\link[base]{set.seed}}, \code{\link[MASS]{glm.nb}}
}
