# CHECKING UTILS
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
check_lindep <- function(form, mf){
  attr(mf, "terms") <- NULL
  new_mf <- rm_lindep(form, mf)
  if(!identical(mf, new_mf)){
    lindep_vars <- setdiff(names(mf), names(new_mf))
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
rm_lindep <- function(form, mf){
  y <- names(mf[1])
  form <- formula(paste(y, "~ ."))
  # Correct non-standard column names
  # names(mf) <- make.names(names(mf))
  mm <- stats::model.matrix(form, mf)
  # factor the matrix using QR decomposition
  qr_ob <- qr(mm)
  # extract R matrix
  R <- qr.R(qr_ob)
  if (is.null(dim(R)[2]) || qr_ob$rank == dim(R)[2]){
    # there are no linear combinations; return original mf (-terms attribute)
    attr(mf, "terms") <- NULL
    mf
  } else {
    # extract the independent columns and remove
    p1 <- 1:qr_ob$rank
    X <- R[p1, p1]
    mm_keep <- which(colnames(mm) %in% colnames(X))
    mm_dict <- mf_to_mm(mf)
    mf_keep <- purrr::map_lgl(mm_dict, ~ all(.x %in% mm_keep))
    rm_lindep(form, mf[c(T, mf_keep)])
  }
}

## Check for problems with p_max and n_folds
check_cv <- function(n, p, binom, n_folds){
  max_folds <- floor(n/5)
  if(max_folds < 2){
    if(binom){
      stop("Your sample size for the minority class is too small.")
    } else {
      stop("Your sample size is too small.")
    }
  }
  alt_folds <- min(n_folds, max_folds)
  alt_p <- floor(n/alt_folds * (alt_folds - 1))
  alt_p <- min(p, alt_p)
  while(alt_folds < max_folds && alt_p < p){
    alt_folds <- alt_folds + 1L
    alt_p <- floor(n/alt_folds * (alt_folds - 1))
    alt_p <- min(p, alt_p)
  }
  list(p = alt_p, folds = alt_folds)
}
