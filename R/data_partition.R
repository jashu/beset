#' \code{data_partition} Constructor
#'
#' Constructs an object of class {data_partition}.
#'
#' A \code{data_partition} object is a list containing exactly two data frames
#' (\code{train} and \code{test}). This object will normally be constructed by
#' passing a single data frame to \code{\link{partition}}. Use this constructor
#' function in the event that you wish to manually link two independent data
#' sets: one to be used for model training and the other to be used for model
#' testing. For your protection, construction will fail if any variables in
#' the training set are not found in the test set.
#'
#' @param train A data frame containing the training set.
#'
#' @param test A data frame containing the test set.
#'
#' @export
#'
data_partition <- function(train, test){
  if(!("data.frame" %in% class(train) && "data.frame" %in% class(test))){
    stop("Both training and test data must be stored as data frames.")
  }
  if(!all(names(train) %in% names(test))){
    stop(paste("All of the variables in the training data must also be",
               "included in the test data."))
  }
  structure(list(train = train, test = test), class = "data_partition")
}

