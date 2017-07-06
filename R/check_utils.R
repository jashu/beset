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
## Check for linear dependence
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
## Check for problems with p_max and n_folds args

check_cv <- function(n, p, binom, n_folds){
  alt_folds <- n / (n - p * 2)
  alt_p <- p
  while(!dplyr::between(alt_folds, 1, 10)){
    alt_p <- alt_p - 1L
    alt_folds <- n / (n - alt_p * 2L)
  }
  if(alt_p < 1){
    if(binom){
      stop("Your sample size for the minority class is too small.")
    } else {
      stop("Your sample size is too small.")
    }
  }
  list(p = as.integer(alt_p), folds = as.integer(alt_folds))
}
