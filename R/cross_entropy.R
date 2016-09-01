

#' @export
cross_entropy <- function(object, test_data){
  y_obs <- unlist(test_data[, names(object$model)[1]])
  y_bar <- mean(y_obs, na.rm = TRUE)
  y_hat <- predict(object, test_data, type="response")
  if(class(object) == "zeroinfl"){
    if(object$dist == "poisson") family <- "zip" else family <- "zinb"
  } else family <- object$family$family
  if(grepl("Negative Binomial", family)) family <- "negbin"
  log_prob <- switch(family,
                 gaussian = dnorm(y_obs, mean = y_hat,
                                  sd = sigma(object), log = TRUE),
                 binomial = y_obs * log(y_hat) + (1 - y_obs) * log(1 - y_hat),
                 poisson = dpois(y_obs, lambda = y_hat, log = TRUE),
                 negbin = dnbinom(y_obs, mu = y_hat, size=object$theta,
                                  log = TRUE),
                 zip = VGAM::dzipois(y_obs,
                                     lambda = predict(object,
                                                      test_data,
                                                      type = "count"),
                                     pstr0 = predict(object,
                                                     test_data,
                                                     type = "zero"),
                                     log = TRUE),
                 zinb = VGAM::dzinegbin(y_obs,
                                        size = object$theta,
                                        munb = predict(object,
                                                       test_data,
                                                       type = "count"),
                                        pstr0 = predict(object,
                                                        test_data,
                                                        type = "zero"),
                                        log = TRUE))
  -1 * mean(log_prob)
}
