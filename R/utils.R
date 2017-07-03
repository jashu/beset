fit_glm <- function(mf, form, family, link){
  if(family == "negbin"){
    suppressWarnings(glm_nb(formula = form,
                            data = mf,
                            link = link))
  } else {
    suppressWarnings(stats::glm(formula = form,
                                family = do.call(family, list(link = link)),
                                data = mf, model = FALSE, y = FALSE))
  }
}
