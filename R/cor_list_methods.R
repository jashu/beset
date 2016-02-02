#' Summarizing Correlation Lists
#'
#' Methods for the \code{cor_list} class.
#'
#' @param object An object of class \code{\link{cor_list}}.
#' @export

summary.cor_list <- function (object,
                              x1 = everything(),
                              x2 = everything(),
                              sort = TRUE){
  bad_args <- try(is.character(x1), silent = T)
  if(class(bad_args) == "logical"){
    stop("Remove the quotation marks from your x1 argument.")
  }
  bad_args <- try(is.character(x2), silent = T)
  if(class(bad_args) == "logical"){
    stop("Remove the quotation marks from your x2 argument.")
  }
  output <- as.data.frame(object)
  var_names <- unique(output$x1)
  x1 <- lazyeval::lazy(x1); x2 <- lazyeval::lazy(x2)
  x1 <- dplyr::select_vars_(var_names, x1)
  x2 <- dplyr::select_vars_(var_names, x2)
  output <- output[output$x1 %in% x1 & output$x2 %in% x2,]
  if (sort) output <- dplyr::arrange(output, x1, desc(output[,3]))
  output[,3:ncol(output)] <- round(output[,3:ncol(output)], 2)
  attr(output, "row.names") <- NULL
  structure(output, class = "cor_list_summary")
}

#' @export
as.data.frame.cor_list <- function (object){
  # j is assigned to first column of data frame because it is the clustered
  # index; i is assigned to second column because it permutes within j
  output <- data.frame(x1 = object$x1,
                       x2 = object$x2,
                       coef = object$coef,
                       stringsAsFactors = FALSE)
  if(!is.null(object$lower)){
    output$lower <- object$lower
    output$upper <- object$upper
  }
  coef_name <- attr(object, "coef")
  if(!is.null(coef_name)) names(output)[3] <- coef_name
  output
}

#' @export
print.cor_list <- function(object){
  temp <- summary(object)
  paths <- paste(temp$x1, temp$x2, sep = " <-> ")
  output <-  data.frame(correlates = paths,
                        coef = temp[[3]],
                        stringsAsFactors = FALSE)
  names(output)[2] <- names(temp)[3]
  rownames(output) <- NULL
  if(nrow(output) > 50){
    print(head(output, 25, addrownums = FALSE))
    if(nrow(output) > 25){
      cat("............................................\n",
          "\t(", nrow(output) - 50, " correlations omitted)\n",
          "............................................\n",sep="")
      print(tail(output, 25, addrownums = FALSE))
    }
  } else print(output, row.names = FALSE)
}

#' @export
print.cor_list_summary <- function(object){
  paths <- paste(object$x1, object$x2, sep = " <-> ")
  x1s <- unique(object$x1)
  for(i in seq_along(x1s)){
    temp <-  data.frame(correlates = paths[object$x1 == x1s[i]],
                        coef = object[[3]][object$x1 == x1s[i]],
                        stringsAsFactors = FALSE)
    if(!is.null(object$lower)){
      temp$lower <- object$lower[object$x1 == x1s[i]]
      temp$upper <- object$upper[object$x1 == x1s[i]]
    }
    names(temp)[2] <- names(object)[3]
    print(temp, row.names = FALSE)
    cat("\n")
  }
}

