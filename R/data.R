#' Self-schema and depression in adolescents
#'
#' A dataset containing diffusion model metrics from the self-referent encoding
#' task (SRET) and total counts from the short form of the Children's Depression
#' Inventory (CDI:S).
#'
#' This is an abridged data set from a study by Dainer-Best et al. (2018) used
#' to illustrate best subset selection with negative binomial models. Complete
#' data from this study can be obtained from the
#' \href{http://dx.doi.org/10.18738/T8/XK5PXX}{Texas Data Repository}.
#'
#' @format A data frame with 408 rows and 11 variables:
#' \describe{
#'   \item{dep}{total depression score from the CDI:S}
#'   \item{v.positive}{drift rate to positive words}
#'   \item{v.negative}{drift rate to negative words}
#'   \item{zr.positive}{relative starting point for positive words}
#'   \item{zr.negative}{relative starting point for negative words}
#'   \item{a}{threshold separation}
#'   \item{t0}{response time constant}
#'   \item{d}{differences in speed of response execution}
#'   \item{szr}{intertrial variability of starting point}
#'   \item{sv}{intertrial variability of drift}
#'   \item{st0}{intertrial variability of non-decisional components}
#' }
#' @references Dainer-Best, J., Lee, H. Y., Shumake, J. D., Yeager, D. S., &
#'   Beevers, C. G. (2018). Determining optimal parameters of the self-referent
#'   encoding task: A large-scale examination of self-referent cognition and
#'   depression. Psychological Assessment. Advance online publication.
#'   \url{http://dx.doi.org/10.1037/pas0000602}
#' @source \url{https://dataverse.tdl.org/api/access/datafile/942}
"adolescents"

#' Prostate Cancer Study
#'
#' Baseline exam results on prostate cancer patients from Dr. Donn Young at The
#' Ohio State University Comprehensive Cancer Center.
#'
#' @format A data frame with 380 rows and 8 variables:
#' \describe{
#'   \item{tumor}{tumor penetration of prostatic capsule, 0 = no penetration}
#'   \item{age}{in years}
#'   \item{race}{white or black}
#'   \item{dpros}{results of the digital rectal exam, 4 levels}
#'   \item{dcaps}{detection of capsular involvement in rectal exam (yes, no)}
#'   \item{psa}{antigen (mg/ml)}
#'   \item{vol}{tumor volume obtained from ultrasound (cm3)}
#'   \item{gleason}{total gleason score (0 - 10)}
#' }
#' @references Hosmer and Lemeshow (2000) Applied Logistic Regression: Second
#' Edition.
#' @source \url{https://raw.github.com/0xdata/h2o/master/smalldata/logreg/prostate.csv}
"prostate"

