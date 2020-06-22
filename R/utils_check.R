#' CHECKING UTILS

# create a dictionary that maps names of model frame to column indices of
# model matrix
mf_to_mm <- function(mf){
  idx_start <- 2L
  lapply(mf[-1], function(x){
    idx <- idx_start
    if(is.factor(x)){
      idx_stop <- idx_start + length(base::levels(x)) - 2L
      idx <- idx_start:idx_stop
      idx_start <<- idx_stop + 1L
    } else {
      idx_start <<- idx_start + 1L
    }
    idx
  })
}


## Error argument
check_error <- function(
  error = c("auto", "btwn_fold_se", "btwn_rep_range", "none")
){
  tryCatch(
    match.arg(error, c("auto", "btwn_fold_se", "btwn_rep_range", "none")),
    error = function(c){
      c$message <- gsub("arg", "error", c$message)
      c$call <- NULL
      stop(c)
    })
}

## Metric argument
check_metric <- function(
  metric = c("auto", "auc", "mae", "mce", "mse", "rsq")
){
  tryCatch(
    match.arg(metric, c("auto", "auc", "mae", "mce", "mse", "rsq")),
    error = function(c){
      c$message <- gsub("arg", "metric", c$message)
      c$call <- NULL
      stop(c)
    })
}

## Variable Names
check_names <- function(var_names){
  new_names <- make.names(var_names)
  bad_names <- var_names[!var_names %in% new_names]
  if(length(bad_names))
    stop(paste(
      "Variable names may contain only letters, numbers, '.', and '_',\n",
      "  and must start with a letter or '.' followed by a letter.\n",
      "  Rename the following variables before proceeding:\n\t",
      paste0(bad_names, collapse = "\n\t"), sep = ""))
  return(NULL)
}

## Check for valid family arg
check_family <- function(family){
  tryCatch(
    match.arg(family, c("binomial", "gaussian", "poisson", "negbin")),
    error = function(c){
      c$message <- gsub("arg", "family", c$message)
      c$call <- NULL
      stop(c)
    }
  )
}
## Check for valid link arg
check_link <- function(family, link){
  tryCatch(
    if(family == "binomial"){
      match.arg(link, c("logit", "probit", "cauchit", "log", "cloglog"))
    } else if(family == "gaussian"){
      match.arg(link, c("identity", "log", "inverse"))
    } else if(family %in% c("negbin", "poisson")){
      match.arg(link, c("log", "sqrt", "identity"))
    },
    error = function(c){
      c$message <- gsub("'arg'", paste("'link' for", family, "family"),
                        c$message)
      c$call <- NULL
      stop(c)
    })
}
## Check for linear dependencies and remove them
check_lindep <- function(mf){
  X <- stats::model.matrix(terms(mf), mf)
  # remove intercept if present
  valid_names <- setdiff(colnames(X), "(Intercept)")
  X <- X[, valid_names, drop = FALSE]
  new_X <- rm_lindep(X)
  mm_keep <- which(colnames(X) %in% colnames(new_X)) + 1
  mm_dict <- mf_to_mm(mf)
  mf_keep <- c(TRUE, purrr::map_lgl(mm_dict, ~ all(.x %in% mm_keep)))
  new_mf <- mf[mf_keep]
  attr(new_mf, "terms") <- NULL
  if(ncol(new_mf) < ncol(mf)){
    lindep_vars <- setdiff(names(mf), colnames(new_mf))
    dependx <- "dependency"; predx <- "predictor"
    if(length(lindep_vars) > 1){
      dependx <- "dependencies"; predx <- "predictors"
    }
    warning(
      paste("Found ", length(lindep_vars), " linear ", dependx, ". ",
            " Removed the following ", predx, ":\n\t",
            paste0(lindep_vars, collapse = "\n\t"),
            sep = ""),
      immediate. = TRUE
    )
  }
  new_mf
}

rm_lindep <- function(X){
  # factor the matrix using QR decomposition
  qr_ob <- qr(X)
  # extract R matrix
  R <- qr.R(qr_ob)
  if (is.null(dim(R)[2]) || qr_ob$rank == dim(R)[2]){
    # there are no linear combinations; return original X
    X
  } else {
    # extract the independent columns and remove
    p1 <- 1:qr_ob$rank
    X <- X[, colnames(R[p1, p1])]
    rm_lindep(X)
  }
}
