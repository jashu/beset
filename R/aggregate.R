#' @importFrom rlang !!
#' @import purrr
aggregate.nested <- function(x, metric, ...){
  stats <- map(x$beset, "stats")
  cv <- map(stats, "cv")
  parameters <- intersect(
    c("lambda", "alpha", "n_pred", "form", metric), names(cv[[1]])
  )
  cv <- map(cv, ~ .x[parameters])
  fit <- map(stats, "fit") %>% map(~.x[parameters])
  test <- map(stats, "test") %>% map(~.x[parameters])
  if("form" %in% names(cv[[1]])){
    fit <- map2(fit, cv, ~.x[.x$form %in% .y$form,] %>% dplyr::select(-form))
    test <- map2(test, cv, ~.x[.x$form %in% .y$form,] %>% dplyr::select(-form))
    cv <- map(cv, ~.x %>% dplyr::select(-form))
  } else {
    for(i in 1:length(cv)){
      n <- nrow(cv[[i]])
      cv[[i]]$n_pred <- fit[[i]]$n_pred <- test[[i]]$n_pred <- 1:n
    }
  }
  metric_sym <- rlang::sym(metric)
  se <- rlang::expr(sd(!!metric_sym) / sqrt(x$n_folds))
  fit <- reduce(fit, rbind)
  test <- reduce(test, rbind)
  cv <- reduce(cv, rbind)
  if(inherits(x, "elnet")){
    fit_means <- fit %>% group_by(n_pred, alpha)
    fit_ses <-  fit %>% group_by(n_pred, alpha)
    test_means <- test %>% group_by(n_pred, alpha)
    test_ses <-  test %>% group_by(n_pred, alpha)
    cv <- cv %>% group_by(n_pred, alpha)
  } else {
    fit_means <- fit %>% group_by(n_pred)
    fit_ses <- fit %>% group_by(n_pred)
    test_means <- test %>% group_by(n_pred)
    test_ses <-  test %>% group_by(n_pred)
    cv <- cv %>% group_by(n_pred)
  }
  fit_means <-  fit_means %>% summarize_all(mean)
  fit_ses <-  fit_ses %>% summarize(se = !!se)
  fit <- fit_means
  fit$train <- fit_means[[metric]]
  fit$train_lower <- fit$train - fit_ses$se
  fit$train_upper <- fit$train + fit_ses$se
  test_means <- test_means %>% summarize_all(mean)
  test_ses <- test_ses %>% summarize(se = !!se)
  test <- test_means
  test$test <- test_means[[metric]]
  test$test_lower <- test$test - test_ses$se
  test$test_upper <- test$test + test_ses$se
  cv$cv <- map_dbl(cv[[metric]], "mean")
  cv$se <- map_dbl(cv[[metric]], "btwn_fold_se")
  fit[[metric]] <- cv[[metric]] <- test[[metric]] <- NULL
  cv <- cv %>% summarize_all(mean) %>%
    mutate(cv_lower = cv - se, cv_upper = cv + se) %>% select(-se)
  join_vars <- intersect(names(fit), names(test))
  data <- inner_join(fit, cv, by = join_vars)
  data <- inner_join(data, test, by = join_vars) %>% ungroup
  if(inherits(x, "elnet")){
    data %>% select(lambda, alpha, train, train_lower, train_upper,
                    cv, cv_lower, cv_upper, test, test_lower, test_upper)
  } else {
    data
  }
}
