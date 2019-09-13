#' Partial-Dependence Effects
#'
#' Generate data for partial dependence plots.
#'
#' @param object A "beset" model object for which partial dependence plots
#' are desired.
#'
#' @param x (Optional) \code{character} vector giving name(s) of predictors
#' for which partial dependence effects should be estimated. If omitted,
#' partial dependence effects will be estimated for all predictors in the
#' model.
#'
#' @param cond (Optional) named \code{list} of the values to use for one or more
#' predictors when estimating the partial effect of the predictor(s) named in
#' \code{x}. Any predictor omitted from this list will be assigned the mean
#' observed value for continuous variables, or the modal level for factors.
#'
#' @param x_lab (Optional) \code{character} vector giving replacement labels
#' for \code{x}. Labels should be listed in the same order as the predictor
#' names in \code{x}, or, if \code{x} is omitted, the same column order as the
#' predictors in the model data frame. If omitted, plots will be labeled with
#' the names of the predictors as they appear in the model object.
#'
#' @param y_lab (Optional) \code{character} string giving replacement label
#' for the response variable in the model. If omitted, plots will be labeled
#' with the name of the response as it appears in the model object.
#'
#' @param ... Additional arguments affecting the plots produced.
#'
#' @param order If plotting a grid, order in which partial dependence plots
#' should appear. Options are \code{"delta"} (default), which arranges plots
#' based on the magnitude of change in \code{y} over the range of each \code{x},
#' and \code{"import"}, which arranges plots based on the variable importance
#' scores returned by the object's \code{\link{importance}} method.
#'
#' @param p_max Maximum number of partial dependence effects to plot.
#' Default is 16, resulting in a 4 X 4 grid.
#'
#' @inheritParams summary.beset
#' @inheritParams beset_glm

#' @export
dependence <- function(object, ...){
  UseMethod("dependence")
}

#' @rdname dependence
#' @export
dependence.nested <- function(object, x = NULL, cond = NULL,
                              alpha = NULL, lambda = NULL, n_pred = NULL,
                              metric = "auto", oneSE = TRUE,
                              x_lab = NULL, y_lab = NULL, ...,
                              parallel_type = NULL, n_cores = NULL, cl = NULL){
  if(inherits(object, "elnet"))
    do.call(dependence.nested_elnet, as.list(match.call())[-1])
  else
    tryCatch(
      error = function(c){
        c$message <- paste("Partial dependence method not yet implemented,",
                           "for `beset_lm` and `beset_glm`.")
        c$call <- NULL
        stop(c)
      }
    )
}

#' @rdname dependence
#' @export
dependence.beset <- function(
  object, x = NULL, cond = NULL, alpha = NULL, lambda = NULL, n_pred = NULL,
  metric = "auto", oneSE = TRUE, x_lab = NULL, y_lab = NULL, ...,
  parallel_type = NULL, n_cores = NULL, cl = NULL
){
  if(inherits(object, "rf")){
    return(dependence.beset_rf(object, x, cond, x_lab, y_lab, ...,
      parallel_type = parallel_type, n_cores = n_cores, cl = cl
    ))
  }
  model_data <- list(...)
  if(is.null(object$data)){
    object$data <- model_data$data
    object$terms <- model_data$terms
    object$contrasts <- model_data$contrasts
    object$xlevels <- model_data$xlevels
  }
  data <- object$data
  if(is.null(data)){
      stop(
        paste("Partial dependence cannot be calculated from skinny objects.",
              "Please rerun `beset_elnet` and set `skinny` to `FALSE`.")
      )
  }
  data <- mutate_if(data, is.logical, factor)
  if(is.null(x)) x <- labels(object$terms)
  if(is.null(x_lab)) x_lab <- x
  y <- setdiff(all.vars(object$terms), labels(object$terms))
  if(is.null(y_lab)) y_lab <- y
  if(length(x) > 1){
    return(map2(x, x_lab,
                ~ dependence.beset(object, x = .x, cond = cond,
                                   alpha = alpha, lambda = lambda,
                                   metric = metric, oneSE = oneSE,
                                   n_pred = n_pred, x_lab = .y, y_lab = y_lab)))
  }
  xes <- setdiff(labels(object$terms), x)
  covars <- data[xes]
  covars <- map_df(covars, function(x){
    if(is.factor(x))
      factor(levels(forcats::fct_infreq(x))[1], levels = levels(x)) else
        mean(x, na.rm = TRUE)
  })
  if(!is.null(cond)){
    for(varname in names(cond)) covars[, varname] <- cond[[varname]]
  }
  x_obs <- data[[x]]
  x_plot <- if(is.factor(x_obs)){
    factor(levels(x_obs), levels = levels(x_obs))
  } else if(is.character(x_obs) ||
            (is.integer(x_obs) && n_distinct(x_obs) < 100)){
    sort(unique(x_obs))
  } else {
    seq(min(x_obs, na.rm = T), max(x_obs, na.rm = T), length.out = 100)
  }
  plot_data <- vector("list", length(x_plot))
  plot_data <- map_df(plot_data, ~ return(covars))
  plot_data <- bind_cols(x = x_plot, plot_data)
  names(plot_data)[1] <- x
  y_hat <- predict(object, newdata = plot_data,
                   newoffset = object$parameters$offset,
                   alpha = alpha, lambda = lambda,
                   metric = metric, oneSE = oneSE)
  y_obs <- object$parameters$y
  if(is.factor(y_obs)) y_obs <- as.integer(y_obs) - 1
  impact <- (max(y_hat) - min(y_hat)) / (max(y_obs) - min(y_obs))
  plot_data <- data_frame(
    x = x_plot,
    y = y_hat
  )
  partial_effects <- data_frame(
    variable = x_lab,
    delta = impact
  )
  structure(
    list(plot_data = plot_data, partial_effects = partial_effects),
    class = "part_depend"
  )
}

dependence.nested_elnet <- function(
  object, x = NULL, cond = NULL, alpha = NULL, lambda = NULL,
  metric = "auto", oneSE = TRUE, x_lab = NULL, y_lab = NULL, ...,
  parallel_type = NULL, n_cores = NULL, cl = NULL
){
  if(is.null(object$data)){
    stop(
      paste("Partial dependence cannot be calculated from skinny objects.",
            "Please rerun `beset_elnet` and set `skinny` to `FALSE`.")
    )
  }
  if(is.null(x)) x <- labels(object$terms)
  y <- setdiff(all.vars(object$terms), labels(object$terms))
  if(is.null(x_lab)) x_lab <- x
  if(is.null(y_lab)) y_lab <- y
  variable_importance <- importance(object, alpha = alpha, lambda = lambda,
                                    metric = metric, oneSE = oneSE)
  keep <- map(x, ~ grep(.x, variable_importance$variable)) %>% as_vector
  variable_importance <- variable_importance[keep,]
  matching_labels <- map(x, ~ grep(.x, variable_importance$variable))
  walk2(x_lab, matching_labels,
       function(x, i) variable_importance$variable[i] <<- x)
  variable_importance <- variable_importance %>% group_by(variable) %>%
    summarize_all(max)
  if(is.null(n_cores) || n_cores > 1){
    parallel_control <- beset:::setup_parallel(
      parallel_type = parallel_type, n_cores = n_cores, cl = cl)
    have_mc <- parallel_control$have_mc
    n_cores <- parallel_control$n_cores
    cl <- parallel_control$cl
  }
  pd <- if(n_cores > 1L){
    if(is.null(cl)){
      parallel::mclapply(object$beset, dependence, x = x, cond = cond,
                         alpha = alpha, lambda = lambda, metric = metric,
                         oneSE = oneSE, x_lab = x_lab, y_lab = y_lab,
                         data = object$data, terms = object$terms,
                         contrasts = object$contrasts,
                         xlevels = object$xlevels, mc.cores = n_cores)
        } else {
          parallel::parLapply(cl, object$beset, dependence, x = x,
                              cond = cond, alpha = alpha, lambda = lambda,
                              metric = metric, oneSE = oneSE,
                              x_lab = x_lab, y_lab = y_lab,
                              data = object$data, terms = object$terms,
                              contrasts = object$contrasts,
                              xlevels = object$xlevels)
        }
      } else {
      lapply(object$beset, dependence, x = x, cond = cond,
             alpha = alpha, lambda = lambda, metric = metric,
             oneSE = oneSE,  x_lab = x_lab, y_lab = y_lab,
             data = object$data, terms = object$terms,
             contrasts = object$contrasts,
             xlevels = object$xlevels)
      }
  if(!is.null(cl)) parallel::stopCluster(cl)
  if(length(x) > 1){
    pd <- transpose(pd)
    names(pd) <- x_lab
    impact <- map(pd, ~ map_dbl(.x, ~.x$partial_effects$delta))
    plot_data <- map(pd, ~ imap_dfr(.x, function(x, i){
      df <- x$plot_data; df$model <- i; df}))
    partial_dependence <- imap(
      plot_data, ~ ggplot(data = .x, aes(x = x, y = y)) +
      geom_line(aes(group = model), alpha = .1) +
      (if(is.factor(.x$x)) {
        stat_summary(aes(group=1), fun.y=mean, geom="line",
                     colour="tomato2", size = 1.5)
      } else {
        geom_smooth(method = "lm", size = 1.5, color = "tomato2")
      }) + xlab(.y) + ylab(y_lab) + theme_bw()
    )
  } else {
    impact <- map_dbl(pd, ~.x$partial_effects$delta)
    plot_data <- imap_dfr(pd, function(x, i){
      df <- x$plot_data; df$model <- i; df})
    partial_dependence <- ggplot(data = plot_data, aes(x = x, y = y)) +
      geom_line(aes(group = model), alpha = .1) +
        (if(is.factor(plot_data$x)) {
          stat_summary(aes(group=1), fun.y=mean, geom="line",
                       colour="tomato2", size = 1.5)
        } else {
          geom_smooth(method = "lm", size = 1.5, color = "tomato2")
        }) + xlab(x_lab) + ylab(y_lab) + theme_bw()
  }
  variable_importance <- variable_importance %>%
    inner_join(
      data_frame(
        variable = names(impact),
        delta = map_dbl(impact, median),
        delta_low = map_dbl(impact, ~ quantile(.x, .025)),
        delta_high = map_dbl(impact, ~ quantile(.x, .975))
      ), by = "variable")
  structure(
    list(partial_dependence = partial_dependence,
         variable_importance = variable_importance),
    class = "part_depend"
  )
}

#' @export
dependence.randomForest <- function(
  object, data = NULL, x = NULL, y = NULL, cond = NULL, x_lab = NULL,
  y_lab = NULL, make_plot = TRUE, ...
){
  if(is.null(x)) x <- labels(object$terms)
  if(is.null(y)) y <- setdiff(all.vars(object$terms), labels(object$terms))
  if(is.null(x_lab)) x_lab <- x
  if(is.null(y_lab)) y_lab <- y
  if(is.null(data)){
    rf_par <- as.list(object$call)
    data <- eval(rf_par$data)
  }
  if(length(x) > 1){
    out <- map2(x, x_lab, ~ dependence.randomForest(
      object, data = data, x = .x, y = y, cond = cond, x_lab = .y,
      y_lab = y_lab, make_plot = make_plot))
    names(out) <- x_lab
    out <- transpose(out)
    out$variable_importance <- reduce(out$variable_importance, rbind)
    class(out) <- "part_depend"
    return(out)
  }
  xes <- setdiff(names(data), x)
  x_obs <- as_vector(data[x])
  covars <- data[xes]
  covars <- map_df(covars, function(x){
    if(is.factor(x))
      factor(levels(forcats::fct_infreq(x))[1], levels = levels(x)) else
        mean(x, na.rm = TRUE)
  })
  if(!is.null(cond)){
    for(varname in names(cond)) covars[, varname] <- cond[[varname]]
  }
  x_plot <- if(is.factor(x_obs)){
    factor(levels(x_obs), levels = levels(x_obs))
  } else if(is.character(x_obs) ||
            (is.integer(x_obs) && n_distinct(x_obs) < 100)){
    sort(unique(x_obs))
  } else {
    seq(min(x_obs, na.rm = T), max(x_obs, na.rm = T), length.out = 100)
  }
  plot_data <- vector("list", length(x_plot))
  plot_data <- map_df(plot_data, ~ return(covars))
  plot_data <- bind_cols(x = x_plot, plot_data)
  names(plot_data)[1] <- x
  type <- if(object$type == "classification") "prob" else "response"
  y_hat <- predict(object, newdata = plot_data, type = type)
  if(inherits(y_hat, "matrix")) y_hat <- y_hat[, "1"]
  y_obs <- object$y
  if(is.factor(y_obs)) y_obs <- as.integer(y_obs) - 1
  impact <- (max(y_hat) - min(y_hat)) / (max(y_obs) - min(y_obs))
  plot_data <- data_frame(
    x = x_plot,
    y = y_hat
  )
  partial_dependence <- if(make_plot){
    ggplot(data = plot_data, aes(x = x, y = y)) +
    geom_line(aes(group = 1)) + theme_classic() + xlab(x_lab) + ylab(y_lab)
  } else plot_data
  variable_importance <- data_frame(
    variable = x_lab,
    delta = impact
  )
  structure(
    list(partial_dependence = partial_dependence,
         variable_importance = variable_importance),
    class = "part_depend"
  )
}

#' @export
#' @rdname dependence
dependence.beset_rf <- function(
  object, x = NULL, cond = NULL, x_lab = NULL, y_lab = NULL, ...,
  parallel_type = NULL, n_cores = NULL, cl = NULL
){
  if(is.null(n_cores) || n_cores > 1){
    parallel_control <- beset:::setup_parallel(
      parallel_type = parallel_type, n_cores = n_cores, cl = cl)
    have_mc <- parallel_control$have_mc
    n_cores <- parallel_control$n_cores
    cl <- parallel_control$cl
  }
  if(is.null(x)) x <- names(object$forests[[1]]$call$x)
  y <- setdiff(names(object$data), names(object$forests[[1]]$call$x))
  if(is.null(x_lab)) x_lab <- x; if(is.null(y_lab)) y_lab <- y
  variable_importance <- importance(object)
  keep <- map(x, ~ grep(.x, variable_importance$variable)) %>% as_vector
  variable_importance <- variable_importance[keep,]
  matching_labels <- map(x, ~ grep(.x, variable_importance$variable))
  walk2(x_lab, matching_labels,
        function(x, i) variable_importance$variable[i] <<- x)
  data <- object$data
  attr(data, "terms") <- NULL
  pd <- if(n_cores > 1L){
    if(have_mc){
      parallel::mclapply(object$forests, beset:::dependence.randomForest,
                         data = data, x = x, y = y, cond = cond,
                         x_lab = x_lab, y_lab = y_lab, make_plot = FALSE,
                         mc.cores = n_cores)
    } else {
      parallel::parLapply(cl, object$forests, beset:::dependence.randomForest,
                          data = data, x = x, y = y, cond = cond,
                          x_lab = x_lab, y_lab = y_lab, make_plot = FALSE)
    }
  } else {
    lapply(object$forests, beset:::dependence.randomForest, data = data,
           x = x, y = y, cond = cond, x_lab = x_lab, y_lab = y_lab,
           make_plot = FALSE)
  }
  if(!is.null(cl)) parallel::stopCluster(cl)
  pd <- transpose(pd)
  impact <- pd$variable_importance %>% map("delta") %>% transpose %>%
    simplify_all
  names(impact) <- x_lab
  plot_data <- pd$partial_dependence
  if(length(x) > 1){
    plot_data <- transpose(plot_data)
    names(plot_data) <- x_lab
    plot_data <- map(plot_data, ~ imap_dfr(.x, function(x, i){
      df <- x; df$model <- i; df}))
  } else {
    plot_data <- list(imap_dfr(plot_data, function(x, i){
      df <- x; df$model <- i; df}))
    names(plot_data) <- x_lab
  }
  partial_dependence <- imap(
    plot_data, ~ ggplot(data = .x, aes(x = x, y = y)) +
      geom_line(aes(group = model), alpha = .1) + (
        if(is.factor(.x$x)) {
          stat_summary(aes(group=1), fun.y=mean, geom="line",
                       colour="tomato2", size = 1.5)
        } else {
          geom_smooth(method = "loess", size = 1.5, color = "tomato2")
        }
      ) + xlab(.y) + ylab(y_lab) + theme_bw()
    )
  variable_importance <- variable_importance %>%
    mutate(delta = map_dbl(impact, median),
           delta_low = map_dbl(impact, ~ quantile(.x, .025)),
           delta_high = map_dbl(impact, ~ quantile(.x, .975)))
  structure(
    list(partial_dependence = partial_dependence,
         variable_importance = variable_importance),
    class = "part_depend"
  )
}

#' @export
#' @rdname dependence
plot.part_depend <- function(x, order = "delta", p_max = 16, ...){
  plot_data <- x$partial_dependence
  n_vars <- length(plot_data)
  if(n_vars == 1) return(plot(plot_data[[1]]))
  varimp <- x$variable_importance
  varimp <- if(order == "import"){
    arrange(varimp, desc(importance))
  } else {
    arrange(varimp, desc(delta))
  }
  min_y <- plot_data %>% map(~.x$data$y) %>% map_dbl(min) %>% min
  max_y <- plot_data %>% map(~.x$data$y) %>% map_dbl(max) %>% max
  p_max <- min(length(plot_data), p_max)
  plots <- map(plot_data[varimp$variable[1:p_max]],
               ~.x + ylab("") +
                 coord_cartesian(ylim = c(min_y, max_y)) +
                 theme(panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(),
                       axis.line = element_line(size = .5),
                       axis.text.x = element_text(size = 10),
                       axis.text.y = element_text(size = 10),
                       axis.title.x = element_text(size = 12),
                       axis.title.y = element_text(size = 12))
  )
  gout <- suppressWarnings(
    gridExtra::arrangeGrob(grobs = plots,
                           left = x$partial_dependence[[1]]$labels$y)
  )
  plot(gout)
}

#' @export
print.part_depend <- function(x, ...){
  plot(x, ...)
}
