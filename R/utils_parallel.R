# Parallel Utilities

setup_parallel <- function(parallel_type = NULL, n_cores = NULL, cl = NULL,
                           ...){
  available_cores <- parallel::detectCores(logical = FALSE) - 1L
  if(is.null(n_cores)){
    n_cores <- available_cores
  } else {
    n_cores <- min(n_cores, available_cores)
  }
  have_mc <- have_snow <- FALSE
  if (n_cores > 1L) {

    if(is.null(parallel_type)){
      parallel_type <- if(.Platform$OS.type == "windows") "sock" else "fork"
    }
    if (parallel_type == "fork"){
      have_mc <- .Platform$OS.type != "windows"
    } else if (parallel_type == "sock"){
      have_snow <- TRUE
    } else stop("parallel type must be either 'fork' or 'sock'")
    if (!have_mc && !have_snow) n_cores <- 1L
  }
  if (have_snow && is.null(cl)) {
    cl <- parallel::makePSOCKcluster(rep("localhost", n_cores))
  }
  list(have_mc = have_mc, n_cores = n_cores, cl = cl)
}
