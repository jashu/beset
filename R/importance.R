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
#' @import ggplot2
#' @import purrr
#' @import dplyr
#' @export

importance <- function(object, ...){
  UseMethod("importance")
}

#' @export
importance.beset <- function(object, ...){
  if(!inherits(object, "rf")){
    tryCatch(
      error = function(c){
        c$message <- paste("To obtain variable importance scores,",
                           "set `nest_cv` to `TRUE` when running `beset_`.")
        c$call <- NULL
        stop(c)
      }
    )
  } else {
    importance.beset_rf(object, ...)
  }
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
  varimp <- tibble(
    variable = names(imp),
    importance = import / scale_by,
    min_import = min_import / scale_by,
    max_import = max_import / scale_by
  ) %>% filter(variable != "(Intercept)")
  if(!is.null(object$xlevels)){
    for(i in seq_along(object$xlevels)){
      root_name <- names(object$xlevels)[i]
      for(x in object$xlevels[[i]]){
        varimp$variable <- gsub(
          paste(root_name, x, sep = ""), root_name, varimp$variable, fixed = T
        )
      }
    }
  }
  varimp <- varimp %>% group_by(variable) %>% summarize_all(max) %>% ungroup
  structure(varimp, class = c("variable_importance", class(varimp)))
}

#' @export
#' @describeIn importance Relative importance method for "beset_rf" objects
importance.beset_rf <- function(object, ...){
  import_name <- c("MeanDecreaseAccuracy", "%IncMSE")
  import_name <- intersect(import_name, colnames(object$forests[[1]]$importance))
  import <- map(object$forests, ~ .x$importance[, import_name]) %>%
    transpose %>% simplify_all %>% map_dbl(mean)
  import_sd <- map(object$forests, ~.x$importanceSD)
  if(inherits(import_sd[[1]], "matrix")){
    import_sd <- map(import_sd, ~.x[, import_name])
  }
  import_sd <- import_sd %>% transpose %>% simplify_all %>% map_dbl(mean)
  scale_by <- sum(import)
  import <- import / scale_by
  import_sd <- import_sd / scale_by
  varimp <- tibble(
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
varimp <- tibble(
  variable = object$variable_importance$variable,
  importance = import / scale_by,
  min_import = min_import / scale_by,
  max_import = max_import / scale_by
)
structure(varimp, class = c("variable_importance", class(varimp)))
}

#' @export
#' @rdname importance
plot.variable_importance <- function(x, p_max = 20, labels = NULL, ...){
  if(!is.null(labels)){
    if(n_distinct(labels[[2]]) < nrow(labels)){
      duplicate_labels <- table(labels[[2]])
      duplicate_labels <- duplicate_labels[duplicate_labels > 1]
      stop(
        paste0("\nYou have provided the same label for more than one variable.\n",
            "Please fix these duplicate labels:\n\t",
            paste0(names(duplicate_labels), collapse = "\n\t"),
            sep = ""),
        call. = FALSE
      )
    }
    matching_labels <- map(labels[[1]], ~ match(., x$variable))
    not_empty <- which(!map_lgl(matching_labels, is.na))
    labels <- labels[not_empty,]
    matching_labels <- matching_labels[not_empty]
    walk2(labels[[2]], matching_labels, function(a, i) x$variable[i] <<- a)
  }
  x <- x %>% arrange(desc(importance), desc(min_import), desc(max_import))
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

