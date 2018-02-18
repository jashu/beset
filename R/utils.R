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
  vars <- gsub("`", "", vars)
  data <- data[vars]
  data[complete.cases(data),]
}
