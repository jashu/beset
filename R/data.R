#' Self-schema and depression in adolescents
#'
#' A dataset containing measures from the self-referent encoding task (SRET) and
#' total counts from the short form of the Children's Depression Inventory
#' (CDI:S).
#'
#' This is an abridged data set from a study by Dainer-Best et al. (2018) used
#' to illustrate best subset selection with negative binomial models. Complete
#' data from this study can be obtained from the
#' \href{http://dx.doi.org/10.18738/T8/XK5PXX}{Texas Data Repository}.
#'
#' @format A data frame with 270 rows and 9 variables:
#' \describe{
#'   \item{dep}{total depression score from the CDI:S}
#'   \item{num.pos.endorsed}{number of positive words endorsed as
#'         self-descriptive}
#'   \item{num.neg.endorsed}{number of negative words endorsed as
#'         self-descriptive}
#'   \item{numnegrecalled}{number of negative words from task that were later
#'         recalled}
#'   \item{numposrecalled}{number of positive words from task that were later
#'         recalled}
#'   \item{numSRnegrecalled}{number of negative words that were both endorsed
#'         and later recalled}
#'   \item{numSRposrecalled}{number of positive words that were both endorsed
#'         and later recalled}
#'   \item{negRT}{mean latency to endorse negative words in seconds}
#'   \item{posRT}{mean latency to endorse positive words in seconds}
#' }
#' @references Dainer-Best, J., Lee, H. Y., Shumake, J. D., Yeager, D. S., &
#'   Beevers, C. G. (2018). Determining optimal parameters of the self-referent
#'   encoding task: A large-scale examination of self-referent cognition and
#'   depression. Psychological Assessment. Advance online publication.
#'   \url{http://dx.doi.org/10.1037/pas0000602}
#' @source \url{https://dataverse.tdl.org/api/access/datafile/935}
"adolescents"
