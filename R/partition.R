#' Partition Data into Training and Test Sets
#'
#' Wrapper for \code{\link[caret]{createDataPartition}}
#'
#' \code{\link[caret]{createDataPartition}} returns row position integers
#' corresponding to the training data, as derived from a random sampling within
#' the factor levels of a classification outcome or the quintiles of a numeric
#' outcome. This insures that the training and test samples will be closely
#' matched in terms of class incidence or frequency distribution of the outcome
#' measure. The typical use case for \code{\link[caret]{createDataPartition}}
#' is to create a single partition, which will include a preceding call to
#' \code{\link[base]{set.seed}}, so that the randomization is reproducible, and
#' a subsequent construction of separate training and testing data frames using
#' the row positions that are returned.
#'
#' \code{partition} consolidates the above steps into a single function call.
#' In addition, it insures that binary response variables are converted to
#' factors using \code{\link{binary_to_factor}} and that the partition can be
#' reproduced by including a randomization seed as a function argument, and it
#' creates analysis-ready \code{train} and \code{test} data frames, which are
#' bound together in a "\code{\link{data_partition}}" structure so that their
#' common ancestry is maintained and self-documented. For example, if you name
#' your \code{\link{data_partition}} "\code{data}", you can intutively access
#' the training set with \code{data$train} and its corresponding test set with
#' \code{data$test}. Moreover, you can lazily pass a
#' \code{\link{data_partition}} object to any predictive modeling function in
#' the \code{\link{beset}} package, and it will automatically choose the
#' appropriate data frame for the task at hand.
#'
#' @param data Data frame to be partitioned.
#'
#' @param y Name of the response variable, without quotes.
#'
#' @param p The fraction of data that should be included in the training set.
#'
#' @param seed A single value, interpreted as an integer, for seeding the random
#' number generator. See \code{\link[base]{set.seed}}.
#'
#' @return An object of class "data_partition": a list containing two data
#' frames named \code{train} and \code{test}, containing the training and
#' testing sets, respectively.
#'
#' @seealso \code{\link[caret]{createDataPartition}},
#' \code{\link{binary_to_factor}}, \code{\link[base]{set.seed}},
#' \code{\link{data_partition}}
#'
#' @export

partition <- function(data, y, p = .75, seed = 42){
  a <- as.list(match.call())
  response <- eval(a$y, data)
  if(length(unique(response)) == 2 && !is.factor(response)){
    response <- binary_to_factor(response)
  }
  set.seed(seed, kind = "default")
  inTrain <- caret::createDataPartition(y = response, p = p, list = F)
  test <- data[-inTrain,]
  train <- data[inTrain,]
  data_partition(train, test, as.character(a$y))
}
