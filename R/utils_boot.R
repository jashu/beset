boot_elnet <- function(data, lambda, n_reps = 1000, seed = 42,
                       parallel_type = NULL, n_cores = NULL, cl = NULL){
  parallel_control <- setup_parallel(parallel_type = parallel_type,
                                     n_cores = n_cores,
                                     cl = cl,
                                     data = data, lambda = lambda)
  have_mc <- parallel_control$have_mc
  n_cores <- parallel_control$n_cores
  cl <- parallel_control$cl
  n <- nrow(data$x)
  fn <- function(i){
    fit_arg <- list(x = data$x[i,], y = data$y[i],
                    weights = data$weights[i], offset = data$offset[i])
    fit_arg <- c(fit_arg,
                 data[setdiff(names(data), c("x", "y", "weights", "offset"))])
    m <- do.call(glmnet::glmnet, fit_arg)
    as.vector(coef(m, s = lambda))
  }
  idx_array <-  ordinary_array(n, n_reps, seed)
  boot_ids <- as_data_frame(idx_array) %>% transpose %>% simplify_all
  res <- if (n_cores > 1L) {
    if (have_mc) {
      parallel::mclapply(boot_ids, fn, mc.cores = n_cores)
    } else {
      parallel::parLapply(cl, boot_ids, fn)
    }
  } else lapply(boot_ids, fn)
  res %>% transpose %>% simplify_all
}

boot_ci <- function(t, t0, n, seed, conf){
  if(t0 == 0) return(c(NA_real_, NA_real_))
  t_o <- t
  t <- t_o[is.finite(t_o)]
  w <- qnorm(sum(t < t0)/length(t))
  if (!is.finite(w))
    stop("estimated adjustment 'w' is infinite")
  alpha <- (1 + c(-conf, conf))/2
  zalpha <- qnorm(alpha)
  L <- empinf_reg(t_o, n, seed)
  a <- sum(L^3)/(6 * sum(L^2)^1.5)
  if (!is.finite(a))
    stop("estimated adjustment 'a' is NA")
  adj_alpha <- pnorm(w + (w + zalpha)/(1 - a * (w + zalpha)))
  qq <- norm_inter(t, adj_alpha)
  qq[, 2L]
}

ordinary_array <- function(n, R, seed){
  # R x n array of bootstrap indices, resampled within strata.
  # This is the function which generates a regular bootstrap array
  # using equal weights within each stratum.
  set.seed(seed)
  output <- sample.int(n, n*R, replace=TRUE)
  dim(output) <- c(R, n)
  output
}

freq_array <- function(i_array){
  # converts R x n array of bootstrap indices into
  # R X n array of bootstrap frequencies
  n <- ncol(i_array)
  result <- t(apply(i_array, 1, tabulate, n))
  result
}

empinf_reg <- function(t, n, seed){
  #  Function to estimate empirical influence values using regression.
  #  This method regresses the observed bootstrap values on the bootstrap
  #  frequencies to estimate the empirical influence values
  t_o <- t
  fins <- seq_along(t)[is.finite(t)]
  t <- t[fins]
  R <- length(t)
  f <- freq_array(ordinary_array(n, length(t_o), seed))[fins, -1L]
  X <- cbind(rep(1, R), f/n)
  m <- glm.fit(x = X, y = t)
  beta <- coefficients(m)[-1L]
  l <- rep(0, n)
  l[-1] <- beta
  l - mean(l)
}

norm_inter <- function(t, alpha){
  #  Interpolation on the normal quantile scale.  For a non-integer
  #  order statistic this function interpolates between the surrounding
  #  order statistics using the normal quantile scale.  See equation
  #  5.8 of Davison and Hinkley (1997)
  t <- t[is.finite(t)]
  R <- length(t)
  rk <- (R+1)*alpha
  k <- trunc(rk)
  inds <- seq_along(k)
  out <- inds
  kvs <- k[k>0 & k<R]
  tstar <- sort(t, partial = sort(union(c(1, R), c(kvs, kvs+1))))
  ints <- (k == rk)
  if (any(ints)) out[inds[ints]] <- tstar[k[inds[ints]]]
  out[k == 0] <- tstar[1L]
  out[k == R] <- tstar[R]
  not <- function(v) xor(rep(TRUE,length(v)),v)
  temp <- inds[not(ints) & k != 0 & k != R]
  temp1 <- qnorm(alpha[temp])
  temp2 <- qnorm(k[temp]/(R+1))
  temp3 <- qnorm((k[temp]+1)/(R+1))
  tk <- tstar[k[temp]]
  tk1 <- tstar[k[temp]+1L]
  out[temp] <- tk + (temp1-temp2)/(temp3-temp2)*(tk1 - tk)
  cbind(round(rk, 2), out)
}

