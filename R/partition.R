#' Partition Data into Training and Test Sets
#'
#' \code{partition} randomly splits a data frame into two data frames,
#' \code{train} and \code{test}, which are returned as a
#' \code{\link{data_partition}} structure.
#'
#' \code{partition} creates a train/test split among the rows of a data frame
#' based on stratified random sampling within the factor levels of a
#' classification outcome or the quintiles of a numeric outcome. This insures
#' that the training and test samples will be closely matched in terms of class
#' incidence or frequency distribution of the outcome measure. \code{partition}
#' includes a \code{seed} argument so that the randomized partitioning is
#' reproducible. The \code{train} and \code{test} data frames are returned
#' bound together in a \code{\link{data_partition}} structure so that their
#' common ancestry is maintained and self-documented. For example, if you name
#' your \code{\link{data_partition}} "\code{data}", you can intuitively access
#' the training set with \code{data$train} and its corresponding test set with
#' \code{data$test}.
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
#' @seealso \code{\link[base]{set.seed}}, \code{\link{data_partition}}
#'
#' @export

partition <- function(data, y, p = .75, seed = 42){
  a <- as.list(match.call())
  y <- eval(a$y, data)
  if(length(unique(y)) == 2 && !is.factor(y)){
    y <- factor(y)
  }
  in_train <- vector("logical", length(y))
  y_orig <- y
  if(is.numeric(y)){
    if(min(y_orig) == 0) y <- y[y != 0]
    cuts <- floor(length(y)/2)
    if(cuts < 2) cuts <- 2
    if(cuts > 5) cuts <- 5
    breaks <- unique(quantile(y, probs = seq(0, 1, length = cuts)))
    y <- cut(y, breaks, include.lowest = TRUE)
    y <- as.integer(y)
    if(min(y_orig) == 0){
      y_orig[y_orig != 0] <- y
      y <- y_orig
    }
  }
  purrr::walk(unique(y), function(x){
    n <- sum(y == x)
    assignments <- vector("logical", n)
    assignments[1:round(p*n)] <- TRUE
    set.seed(seed, kind = "default")
    in_train[y == x] <<- sample(assignments)
  })
  train <- data[in_train,]
  test <- data[!in_train,]
  data_partition(train, test, as.character(a$y))
}
