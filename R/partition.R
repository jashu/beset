#' Partition Data into Training and Test Model Frames
#'
#' \code{partition} randomly splits a data frame into two model frames,
#' \code{train} and \code{test}, which are returned as a
#' "data_partition" structure.
#'
#' \code{partition} creates a train/test split among the rows of a data frame
#' based on stratified random sampling within the factor levels of a
#' classification outcome or the quartiles of a numeric outcome. This insures
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
#' @param frac The fraction of data that should be included in the training set.
#'   Default is \code{0.5}.
#'
#' @param seed \code{Integer} value for seeding the random number generator. See
#'   \code{\link[base]{set.seed}}.
#'
#' @inheritParams data_partition
#'
#' @return An object of class "data_partition": a list containing two model
#' frames named \code{train} and \code{test}, containing the training and
#' testing sets, respectively.
#'
#' @examples
#' data <- mtcars
#' factor_names <- c("cyl", "vs", "am", "gear", "carb")
#' data[factor_names] <- purrr::map_dfc(data[factor_names], factor)
#' data <- partition(data, y = "mpg")
#'
#' @seealso \code{\link[base]{set.seed}}, \code{\link{data_partition}}
#'
#' @export

partition <- function(data, y, frac = 0.5, x = NULL, offset = NULL,
                      weights = NULL, na_action = na.omit, seed = 42){
  check_names(names(data))
  if(is.null(x)) x <- setdiff(names(data), c(y, offset, weights))
  missing_y <- is.na(data[[y]])
  if(any(missing_y))
    warning(paste(sum(missing_y), "rows missing response data were dropped."))
  vars <- c(y, x, offset, weights)
  data <- data[!missing_y, vars]
  y_obs <- data[[y]]
  if(length(unique(y_obs)) == 2 && !is.factor(y_obs)){
    y_obs <- factor(y_obs)
  }
  in_train <- vector("logical", length(y_obs))
  strata <- y_obs
  if(is.numeric(y_obs)){
    if(min(y_obs) == 0) y_obs <- y_obs[y_obs != 0]
    cuts <- floor(length(y_obs)/2)
    if(cuts < 2) cuts <- 2
    if(cuts > 5) cuts <- 5
    breaks <- unique(quantile(y_obs, probs = seq(0, 1, length = cuts)))
    y_obs <- cut(y_obs, breaks, include.lowest = TRUE)
    y_obs <- as.integer(y_obs)
    if(min(strata) == 0){
      strata[strata != 0] <- y_obs
    } else {
      strata <- y_obs
    }
  }
  set.seed(seed, kind = "default")
  purrr::walk(unique(strata), function(x){
    n <- sum(strata == x)
    assignments <- vector("logical", n)
    assignments[1:round(frac * n)] <- TRUE
    in_train[strata == x] <<- sample(assignments)
  })
  data_partition(train = data[in_train,], test = data[!in_train,], y = y,
                 x = x, offset = offset, weights = weights,
                 na_action = na_action)
}

