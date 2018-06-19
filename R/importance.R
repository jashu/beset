#' Relative Variable Importance
#'
#' \code{importance} is a generic function for obtaining resampled variable
#' importance scores. The function invokes particular
#' \code{\link[utils]{methods}} which depend on the \code{\link[base]{class}}
#' of the first argument.
#'
#' @param object A model object for which variable importance scores are
#' desired.
#'
#' @param p_max Maximum number of predictors to include in variable importance
#' plot
#'
#' @param labels (Optional) two-column \code{data.frame} where column 1
#' gives the variable names used in the model and column 2 gives a corresponding
#' descriptive label. If \code{labels} are defined, the variable importance plot
#' will replace the model variable names with their descriptive labels.
#'
#' @param ... Additional named arguments that define the model selection rules.
#' See \code{\link{summary.beset}}.
#'
#' @examples
#' rf <- beset_rf(Fertility ~ ., data = swiss)
#' importance(rf)
#'
#' @export

importance <- function(object, ...){
  UseMethod("importance")
}

#' @export
importance.beset <- function(object, ...){
  tryCatch(
    error = function(c){
      c$message <- paste("To obtain variable importance scores,",
                         "set `nest_cv` to `TRUE` when running `beset_`.")
      c$call <- NULL
      stop(c)
    }
  )
}

#' @describeIn importance Determine variable importance for "nested"
#' \code{\link{beset_glm}} and \code{\link{beset_elnet}} objects
#' @export
importance.nested <- function(object, ...){
  imp <- summary(object, ...)$coefs$standardized
  import <- map_dbl(imp, ~ abs(.x$mean))
  min_import <- map_dbl(import - map_dbl(imp, "btwn_fold_se"), ~ max(.x, 0))
  max_import <- map_dbl(import + map_dbl(imp, "btwn_fold_se"), ~ max(.x, 0))
  scale_by <- sum(import)
  varimp <- data_frame(
    variable = names(imp),
    importance = import / scale_by,
    min_import = min_import / scale_by,
    max_import = max_import / scale_by
  ) %>% filter(variable != "(Intercept)")
  structure(varimp, class = c("variable_importance", class(varimp)))
}

#' @export
#' @describeIn importance Relative importance method for "beset_rf" objects
importance.beset_rf <- function(object, ...){
  import <- map(object$forests, ~ .x$importance[,1]) %>% transpose %>%
    simplify_all %>% map_dbl(mean)
  import_sd <- map(object$forests, ~ .x$importanceSD) %>% transpose %>%
    simplify_all %>% map_dbl(mean)
  max_idx <- which.max(import)
  scale_by <- sum(import)
  import <- import / scale_by
  import_sd <- import_sd / scale_by
  varimp <- data_frame(
    variable = names(import),
    importance = import,
    min_import = import - import_sd,
    max_import = import + import_sd
  )
  structure(varimp, class = c("variable_importance", class(varimp)))
}

#' @export
#' @describeIn importance Relative importance based on partial dependence
importance.part_depend <- function(object, ...){
import <- object$variable_importance$delta
min_import <- object$variable_importance$delta_low
max_import <- object$variable_importance$delta_high
scale_by <- sum(import)
varimp <- data_frame(
  variable = object$variable_importance$variable,
  importance = import / scale_by,
  min_import = min_import / scale_by,
  max_import = max_import / scale_by
)
structure(varimp, class = c("variable_importance", class(varimp)))
}

#' @export
#' @describeIn importance Plot method for "variable_importance" objects
plot.variable_importance <- function(x, p_max = 20, labels = NULL){
  if(!is.null(labels)){
    matching_labels <- map(labels[[1]], ~ grep(., x$variable))
    not_empty <- which(map(matching_labels, length) > 0)
    labels <- labels[not_empty,]
    matching_labels <- matching_labels[not_empty]
    walk2(labels[[2]], matching_labels, function(x, i) x$variable[i] <<- x)
  }
  x <- x %>% group_by(variable) %>% summarize_all(max) %>%
    arrange(desc(importance), desc(min_import), desc(max_import))
  impvar <- x$variable
  x$variable <- factor(x$variable, levels = rev(impvar))
  theme_set(theme_classic())
  p_max <- min(nrow(x), p_max)
  p <- ggplot(data = x[1:p_max,],
              aes(x = variable, y = importance)) +
    geom_point(col = "tomato2", size = 3) +
    geom_segment(aes(x = variable, xend = variable,
                     y = 0, yend = max(x$max_import)),
                 linetype = "dashed",
                 size = 0.1) +
    geom_segment(aes(x = variable, xend = variable,
                     y = min_import, yend = max_import),
                 linetype = "solid",
                 size = 1, col = "tomato2") +
    labs(x = "", y = "Relative variable importance") +
    coord_flip() + theme(plot.title = element_text(size = 20))
  p
}

