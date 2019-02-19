#' "data_partition" Constructor
#'
#' Constructs an object of class \code{data_partition}.
#'
#' A \code{data_partition} object is a list containing exactly two data frames
#' (\code{train} and \code{test}). This object will normally be constructed by
#' passing a single data frame to \code{\link{partition}}. Use this constructor
#' function in the event that you wish to manually link two independent data
#' sets: one to be used for model training and the other to be used for model
#' testing.
#'
#' \code{data_partition} objects can be passed as the \code{data} argument to
#' the \code{beset} modeling functions (\code{\link{beset_lm}},
#' \code{\link{beset_glm}}, and \code{\link{beset_elnet}}), in which case these
#' functions will train and cross-validate models using the \code{train} data
#' and append additional evaluation metrics using the \code{test} data. Note
#' that in earlier development versions, these functions provided an
#' optional \code{test_data} argument for this purpose. This has been removed
#' and you are now required to construct a \code{data_partition} object
#' beforehand because the \code{data_partition} constructor performs a number
#' of important checks to insure that your \code{test} data are compatible with
#' your \code{train} data: 1) all predictor and response variables are present
#' in both data sets, 2) the levels of all factor variables are the same for
#' both data sets, 3) if an offset variable is used for model training, an
#' offset variable is provided for predicting the test data, and 4) unless
#' \code{na_action} is set to \code{na.pass}, both data frames contain
#' complete cases with no missing data. The \code{data_partition} constructor
#' will alert you to potential issues, attempt to resolve them, and return an
#' error if it can't.
#'
#' @param train A \code{\link[base]{data.frame}} containing the training data to
#' be used for model fitting.
#'
#' @param test A \code{\link[base]{data.frame}} containing the test data to be
#' used for model validation.
#'
#' @param x (Optional) \code{character} vector giving the names of the columns
#' containing the predictor variables. If omitted, defaults to all columns
#' other than those named as \code{y}, \code{offset}, or \code{weights}.
#'
#' @param y \code{Character} string giving the name of the column containing the
#' response variable to be predicted.
#'
#' @param offset (Optional) \code{character} string giving the name of the
#' column containing a model offset. An offset is a known predictor that is
#' added to a linear model \emph{as is} (with a beta coefficient of 1) rather
#' than having its beta coefficient optimized. If given, an offset must be
#' included for both the \code{train} and \code{test} data frames.
#'
#' @param weights (Optional) \code{character} string giving the name of the
#' column containing observation weights. Use these if you want some rows of the
#' data frame to exert more or less influence than others on a model fit. If
#' given, the \code{weights} column is only applied during model training;
#' a \code{weights} column in the \code{test} data will be ignored.
#'
#' @param na_action \code{Function} defining how \code{NA}s shoud be treated.
#' Options include \code{\link[stats]{na.omit}} (default),
#' \code{\link[stats]{na.fail}}, \code{\link[stats]{na.exclude}}, and
#' \code{\link[stats]{na.pass}}.
#'
#' @return A \code{data_partition} object containing a \code{train} data frame
#' and a \code{test} data frame.
#'
#' @examples
#' train <- mtcars[1:16,]
#' test <- mtcars[17:32,]
#' factor_names <- c("cyl", "vs", "am", "gear", "carb")
#' train[factor_names] <- purrr::map_dfc(train[factor_names], factor)
#' test[factor_names] <- purrr::map_dfc(test[factor_names], factor)
#' data <- data_partition(train, test, "mpg")
#'
#' @import purrr
#' @export

data_partition <- function(train, test, y, x = NULL, offset = NULL,
                           weights = NULL, na_action = na.omit){
  if(!("data.frame" %in% class(train) && "data.frame" %in% class(test))){
    stop("Both training and test data must be stored as data frames.")
  }
  check_names(names(train)); check_names(names(test))
  if(is.null(x)) x <- setdiff(names(train), c(y, offset, weights))
  all_names <- c(y, x, offset, weights)
  not_in_train <- all_names[!hasName(train, all_names)]
  if(length(not_in_train)){
    stop(paste("The following variables were not found in `train`:\n\t",
               paste0(not_in_train, collapse = "\n\t"), sep = ""))
  }
  not_in_test <- all_names[!hasName(test, all_names)]
  not_in_test <- setdiff(not_in_test, weights)
  if(length(not_in_test)){
    stop(paste("The following variables were not found in `test`:\n\t",
               paste0(not_in_test, collapse = "\n\t"), sep = ""))
  }
  factors_in_train <- map(train[map_lgl(train, is.factor)], levels)
  factors_in_test <- map(test[map_lgl(test, is.factor)], levels)
  in_train_not_test <- map2(factors_in_train, factors_in_test,
                            ~ setdiff(.x, .y))
  if(any(map_int(in_train_not_test, length))){
    warn_data <- tibble(
      factor = names(in_train_not_test),
      `levels unobserved in test` = map_chr(in_train_not_test,
                                            ~ paste0(.x, collapse = ", "))
    ) %>% filter(nchar(`levels unobserved in test`) > 0)
    warning(paste(
      "The following factor levels were observed in the train set\n",
      "but not the test set:\n\t",
      paste0(
        apply(warn_data, 1, function(x) paste0(x, collapse = "\t")),
        collapse = "\n\t"),
      sep = ""))
  }
  in_test_not_train <- map2(factors_in_train, factors_in_test,
                            ~ setdiff(.y, .x))
  if(any(map_int(in_test_not_train, length))){
    warn_data <- tibble(
      factor = names(in_test_not_train),
      `levels dropped from test` = map_chr(in_test_not_train,
                                            ~ paste0(.x, collapse = ", "))
    ) %>% filter(nchar(`levels dropped from test`) > 0)
    warning(paste(
      "The following factor levels were dropped from the test set\n",
      "because they were not observed in the train set:\n\t",
      paste0(
        apply(warn_data, 1, function(x) paste0(x, collapse = "\t")),
        collapse = "\n\t"),
      sep = ""))
    test[names(factors_in_train)] <- imap_dfc(
      factors_in_train, ~ factor(test[[.y]], levels = .x))
  }
  form <- formula(paste(y, paste0(x, collapse = " + "), sep = " ~ "))
  train_offset <- if(is.null(offset))
    rep(0, nrow(train)) else train[[offset]]
  train_weights <- if(is.null(weights))
    rep(1, nrow(train)) else train[[weights]]
  train_mf <- model.frame(formula = form, data = train, na.action = na_action,
                          drop.unused.levels = FALSE, offset = train_offset,
                          weights = train_weights)
  # Warn user if any rows were dropped
  n_drop <- nrow(train) - nrow(train_mf)
  if(n_drop > 0){
    warning(paste(
      n_drop, "rows with missing data were dropped from `train`.\n"))
  }
  test_offset <- if(is.null(offset))
    rep(0, nrow(test)) else test[[offset]]
  test_mf <- model.frame(formula = form, data = test, na.action = na_action,
                         drop.unused.levels = FALSE, xlev = factors_in_train,
                         offset = test_offset)
  n_drop <- nrow(test) - nrow(test_mf)
  if(n_drop > 0){
    warning(paste(
      n_drop, "rows with missing data were dropped from `test`.\n"))
  }
  structure(list(train = train_mf, test = test_mf), class = "data_partition")
}
