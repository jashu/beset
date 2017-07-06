fit_glm <- function(mf, form, family, link){
  if(family == "negbin"){
    suppressWarnings(glm_nb(formula = form, data = mf, link = link))
  } else {
    suppressWarnings(stats::glm(formula = form,
                                family = do.call(family, list(link = link)),
                                data = mf, model = FALSE, y = FALSE))
  }
}

# create a dictionary that maps names of model frame to column indices of
# model matrix
mf_to_mm <- function(mf){
  idx_start <- 2L
  lapply(mf[-1], function(x){
    idx <- idx_start
    if(is.factor(x)){
      idx_stop <- idx_start + length(levels(x)) - 2L
      idx <- idx_start:idx_stop
      idx_start <<- idx_stop + 1L
    } else {
      idx_start <<- idx_start + 1L
    }
    idx
  })
}

# simplified version of model.frame for internal use that won't break under
# non-syntatic variable names
model_frame <- function (form, data){
  vars <- rownames(attr(terms(form, data = data), "factors"))
  vars <- gsub("`\\\\", "", vars)
  vars <- gsub("\\\\`", "", vars)
  data <- data[vars]
  data[complete.cases(data),]
}

# check for linear dependencies and remove them
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
