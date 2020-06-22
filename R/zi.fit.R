zi.fit <- function(x, y, z, weights, offset, dist, link, control){
  ziPoisson <- function(parms, trunc.start = FALSE) {
    mu <- as.vector(exp(X %*% parms[1:kx] + offsetx))
    if (trunc.start)
      phi <- rep(0, length(mu))
    else phi <- as.vector(linkinv(Z %*% parms[(kx + 1):(kx + kz)] + offsetz))
    loglik0 <- log(phi + exp(log(1 - phi) - mu))
    loglik1 <- log(1 - phi) + dpois(Y, lambda = mu, log = TRUE)
    if (trunc.start){
      sum(weights[Y1] * loglik1[Y1]) -
        sum(weights[Y1] * log(1 - exp(loglik0[Y1])))
    } else sum(weights[Y0] * loglik0[Y0]) + sum(weights[Y1] * loglik1[Y1])
  }
  ziNegBin <- function(parms, trunc.start = FALSE) {
    mu <- as.vector(exp(X %*% parms[1:kx] + offsetx))
    if (trunc.start)
      phi <- rep(0, length(mu))
    else phi <- as.vector(linkinv(Z %*% parms[(kx + 1):(kx + kz)] + offsetz))
    theta <- exp(parms[(kx + kz) + 1])
    loglik0 <- log(
      phi + exp(
        log(1 - phi) +
          suppressWarnings(dnbinom(0, size = theta, mu = mu, log = TRUE))
        )
      )
    loglik1 <- log(1 - phi) +
      suppressWarnings(dnbinom(Y, size = theta, mu = mu, log = TRUE))
    if (trunc.start){
      sum(weights[Y1] * loglik1[Y1]) -
        sum(weights[Y1] * log(1 - exp(loglik0[Y1])))
    } else {
      sum(weights[Y0] * loglik0[Y0]) + sum(weights[Y1] * loglik1[Y1])
    }
  }
  ziGeom <- function(parms, trunc.start = FALSE){
    ziNegBin(c(parms, 0), trunc.start)
  }
  countGradPoisson <- function(parms) {
    eta <- as.vector(X %*% parms[1:kx] + offsetx)[Y1]
    mu <- exp(eta)
    colSums(
      (
        (Y[Y1] - mu) -
          exp(
            ppois(0, lambda = mu, log.p = TRUE) -
              ppois(0, lambda = mu, lower.tail = FALSE, log.p = TRUE) + eta
            )
      ) * weights[Y1] * X[Y1, , drop = FALSE]
    )
  }
  countGradGeom <- function(parms) {
    eta <- as.vector(X %*% parms[1:kx] + offsetx)[Y1]
    mu <- exp(eta)
    colSums(
      (
        (Y[Y1] - mu * (Y[Y1] + 1)/(mu + 1)) -
          exp(
            pnbinom(0, mu = mu, size = 1, log.p = TRUE) -
              pnbinom(0, mu = mu, size = 1, lower.tail = FALSE, log.p = TRUE) -
              log(mu + 1) + eta
            )
      ) * weights[Y1] * X[Y1, , drop = FALSE]
    )
  }
  countGradNegBin <- function(parms) {
    eta <- as.vector(X %*% parms[1:kx] + offsetx)[Y1]
    mu <- exp(eta)
    theta <- exp(parms[kx + 1])
    logratio <- pnbinom(0, mu = mu, size = theta, log.p = TRUE) -
      pnbinom(0, mu = mu, size = theta, lower.tail = FALSE,
              log.p = TRUE)
    rval <- colSums(
      (
        (Y[Y1] - mu * (Y[Y1] + theta)/(mu + theta)) -
          exp(logratio + log(theta) - log(mu + theta) + eta)
      ) * weights[Y1] * X[Y1, , drop = FALSE]
    )
    rval2 <- sum(
      (
        digamma(Y[Y1] + theta) - digamma(theta) + log(theta) - log(mu + theta) +
          1 - (Y[Y1] + theta)/(mu + theta) + exp(logratio) *
          (log(theta) - log(mu + theta) + 1 - theta/(mu + theta))
      ) * weights[Y1]) * theta
    c(rval, rval2)
  }
  gradPoisson <- function(parms) {
    eta <- as.vector(X %*% parms[1:kx] + offsetx)
    mu <- exp(eta)
    etaz <- as.vector(Z %*% parms[(kx + 1):(kx + kz)] + offsetz)
    muz <- linkinv(etaz)
    clogdens0 <- -mu
    dens0 <- muz * (1 - as.numeric(Y1)) + exp(log(1 - muz) + clogdens0)
    wres_count <- ifelse(
      Y1, Y - mu, -exp(-log(dens0) + log(1 - muz) + clogdens0 + log(mu))
    )
    wres_zero <- ifelse(
      Y1, -1/(1 - muz) * linkobj$mu.eta(etaz),
      (linkobj$mu.eta(etaz) - exp(clogdens0) * linkobj$mu.eta(etaz))/dens0
    )
    colSums(cbind(wres_count * weights * X, wres_zero * weights * Z))
  }
  gradGeom <- function(parms) {
    eta <- as.vector(X %*% parms[1:kx] + offsetx)
    mu <- exp(eta)
    etaz <- as.vector(Z %*% parms[(kx + 1):(kx + kz)] + offsetz)
    muz <- linkinv(etaz)
    clogdens0 <- dnbinom(0, size = 1, mu = mu, log = TRUE)
    dens0 <- muz * (1 - as.numeric(Y1)) + exp(log(1 - muz) + clogdens0)
    wres_count <- ifelse(
      Y1, Y - mu * (Y + 1)/(mu + 1),
      -exp(-log(dens0) + log(1 - muz) + clogdens0 - log(mu + 1) + log(mu))
    )
    wres_zero <- ifelse(
      Y1, -1/(1 - muz) * linkobj$mu.eta(etaz),
      (linkobj$mu.eta(etaz) - exp(clogdens0) * linkobj$mu.eta(etaz))/dens0
    )
    colSums(cbind(wres_count * weights * X, wres_zero * weights * Z))
  }
  gradNegBin <- function(parms) {
    eta <- as.vector(X %*% parms[1:kx] + offsetx)
    mu <- exp(eta)
    etaz <- as.vector(Z %*% parms[(kx + 1):(kx + kz)] + offsetz)
    muz <- linkinv(etaz)
    theta <- exp(parms[(kx + kz) + 1])
    clogdens0 <- dnbinom(0, size = theta, mu = mu, log = TRUE)
    dens0 <- muz * (1 - as.numeric(Y1)) + exp(log(1 - muz) + clogdens0)
    wres_count <- ifelse(
      Y1, Y - mu * (Y + theta)/(mu + theta),
      -exp(-log(dens0) + log(1 - muz) + clogdens0 + log(theta) -
             log(mu + theta) + log(mu))
    )
    wres_zero <- ifelse(
      Y1, -1/(1 - muz) * linkobj$mu.eta(etaz),
      (linkobj$mu.eta(etaz) - exp(clogdens0) * linkobj$mu.eta(etaz))/dens0
    )
    wres_theta <- theta * ifelse(
      Y1,
      digamma(Y + theta) - digamma(theta) + log(theta) - log(mu + theta) + 1 -
        (Y + theta)/(mu + theta),
      exp(-log(dens0) + log(1 - muz) + clogdens0) *
        (log(theta) - log(mu + theta) + 1 - theta/(mu + theta))
    )
    colSums(
      cbind(wres_count * weights * X, wres_zero * weights * Z, wres_theta)
    )
  }
  loglikfun <- switch(
    dist, poisson = ziPoisson, geometric = ziGeom, negbin = ziNegBin
  )
  gradfun <- switch(
    dist, poisson = gradPoisson, geometric = gradGeom, negbin = gradNegBin
  )
  linkobj <- make.link(link)
  linkinv <- linkobj$linkinv
  X <- x
  Z <- z
  Y <- y
  Y <- as.integer(round(Y + 0.001))
  n <- length(Y)
  kx <- NCOL(X)
  kz <- NCOL(Z)
  Y0 <- Y <= 0
  Y1 <- Y > 0
  offsetx <- offset$count
  if(is.null(offsetx)) offsetx <- 0
  if(length(offsetx) == 1) offsetx <- rep.int(offsetx, n)
  offsetx <- as.vector(offsetx)
  offsetz <- offset$zero
  if(is.null(offsetz)) offsetz <- 0
  if(length(offsetz) == 1) offsetz <- rep.int(offsetz, n)
  offsetz <- as.vector(offsetz)
  start <- control$start
  method <- control$method
  hessian <- control$hessian
  ocontrol <- control
  control$method <- control$hessian <- control$EM <- control$start <- NULL
  if(is.null(start)){
    model_zero <- glm.fit(
      Z, as.integer(Y0), weights = weights, family = binomial(link = link),
      offset = offsetz
    )
    countloglikfun <- function(parms, ...){
      loglikfun(
        c(parms[1:kx], rep(0, kz), parms[-(1:kx)]), trunc.start = TRUE
      )
    }
    countgradfun <- switch(
      dist, poisson = countGradPoisson, geometric = countGradGeom,
      negbin = countGradNegBin
    )
    lmstart <- lm.wfit(
      X[Y1, , drop = FALSE], log(Y[Y1]) - offsetx[Y1], weights[Y1]
    )$coefficients
    lmstart <- ifelse(is.na(lmstart), 0, lmstart)
    fit <- tryCatch(
      optim(
        par = c(lmstart, if (dist == "negbin") 0 else NULL),
        fn = countloglikfun, gr = countgradfun,
        method = method, control = control, hessian = FALSE,
      ),
      error = function(e) list(convergence = 1)
    )
    if(fit$convergence > 0){
      model_count <- glm.fit(
        X, Y, family = poisson(), weights = weights, offset = offsetx
      )
      start <- list(
        count = model_count$coefficients, zero = model_zero$coefficients
      )
      if(dist == "negbin") start$theta <- 1
    } else {
      start <- list(count = fit$par[1:kx], zero = model_zero$coefficients)
      if (length(fit$par) > kx) start$theta <- exp(fit$par[-(1:kx)])
    }
    if(ocontrol$EM & dist == "poisson"){
      mui <- model_count$fitted
      probi <- model_zero$fitted
      probi <- probi/(probi + (1 - probi) * dpois(0, mui))
      probi[Y1] <- 0
      ll_new <- loglikfun(c(start$count, start$zero))
      ll_old <- 2 * ll_new
      while (abs((ll_old - ll_new)/ll_old) > control$reltol){
        ll_old <- ll_new
        model_count <- glm.fit(
          X, Y, weights = weights * (1 - probi), offset = offsetx,
          family = poisson(), start = start$count
        )
        model_zero <- suppressWarnings(
          glm.fit(
            Z, probi, weights = weights, offset = offsetz,
            family = binomial(link = link), start = start$zero
          )
        )
        mui <- model_count$fitted
        probi <- model_zero$fitted
        probi <- probi/(probi + (1 - probi) * dpois(0, mui))
        probi[Y1] <- 0
        start <- list(count = model_count$coefficients,
                      zero = model_zero$coefficients)
        ll_new <- loglikfun(c(start$count, start$zero))
      }
    }
    if(ocontrol$EM & dist == "geometric"){
      mui <- model_count$fitted
      probi <- model_zero$fitted
      probi <- probi/(probi + (1 - probi) * dnbinom(0, size = 1, mu = mui))
      probi[Y1] <- 0
      ll_new <- loglikfun(c(start$count, start$zero))
      ll_old <- 2 * ll_new
      while (abs((ll_old - ll_new)/ll_old) > control$reltol) {
        ll_old <- ll_new
        model_count <- suppressWarnings(
          glm.fit(
            X, Y, weights = weights * (1 - probi), offset = offsetx,
            family = MASS::negative.binomial(1), start = start$count
          )
        )
        model_zero <- suppressWarnings(
          glm.fit(
            Z, probi, weights = weights, offset = offsetz,
            family = binomial(link = link), start = start$zero
          )
        )
        start <- list(count = model_count$coefficients,
                      zero = model_zero$coefficients)
        mui <- model_count$fitted
        probi <- model_zero$fitted
        probi <- probi/(probi + (1 - probi) * dnbinom(0, size = 1, mu = mui))
        probi[Y1] <- 0
        ll_new <- loglikfun(c(start$count, start$zero))
      }
    }
    if (ocontrol$EM & dist == "negbin") {
      mui <- model_count$fitted
      probi <- model_zero$fitted
      probi <- probi /
        (probi + (1 - probi) * dnbinom(0, size = start$theta, mu = mui))
      probi[Y1] <- 0
      ll_new <- loglikfun(c(start$count, start$zero, log(start$theta)))
      ll_old <- 2 * ll_new
      offset <- offsetx
      while (abs((ll_old - ll_new)/ll_old) > control$reltol) {
        ll_old <- ll_new
        model_count <- suppressWarnings(
          MASS::glm.nb(
            Y ~ 0 + X + offset(offset), weights = weights *
              (1 - probi), start = start$count, init.theta = start$theta
          )
        )
        model_zero <- suppressWarnings(
          glm.fit(
            Z, probi, weights = weights, offset = offsetz,
            family = binomial(link = link), start = start$zero
          )
        )
        start <- list(count = model_count$coefficients,
                      zero = model_zero$coefficients, theta = model_count$theta)
        mui <- model_count$fitted
        probi <- model_zero$fitted
        probi <- probi /
          (probi + (1 - probi) * dnbinom(0, size = start$theta, mu = mui))
        probi[Y1] <- 0
        ll_new <- loglikfun(c(start$count, start$zero, log(start$theta)))
      }
    }
  }
  fit <- optim(
    par = c(
      start$count, start$zero, if(dist == "negbin") log(start$theta) else NULL
    ),
    fn = loglikfun, gr = gradfun,method = method, control = control,
    hessian = hessian
  )
  coefc <- fit$par[1:kx]
  names(coefc) <- names(start$count) <- colnames(X)
  coefz <- fit$par[(kx + 1):(kx + kz)]
  names(coefz) <- names(start$zero) <- colnames(Z)
  vc <- tryCatch(
    -solve(as.matrix(fit$hessian)), error = function(e) {
      warning(e$message, call = FALSE)
      k <- nrow(as.matrix(fit$hessian))
      return(matrix(NA, k, k))
    }
  )
  if (dist == "negbin") {
    np <- kx + kz + 1
    theta <- as.vector(exp(fit$par[np]))
    SE.logtheta <- as.vector(sqrt(diag(vc)[np]))
    vc <- vc[-np, -np, drop = FALSE]
  } else {
    theta <- NULL
    SE.logtheta <- NULL
  }
  colnames(vc) <- rownames(vc) <- c(
    paste("count", colnames(X), sep = "_"),
    paste("zero", colnames(Z), sep = "_")
  )
  mu <- exp(X %*% coefc + offsetx)[, 1]
  phi <- linkinv(Z %*% coefz + offsetz)[, 1]
  Yhat <- (1 - phi) * mu
  res <- sqrt(weights) * (Y - Yhat)
  nobs <- sum(weights > 0)
  list(
    coefficients = list(count = coefc, zero = coefz),
    residuals = res,
    fitted.values = Yhat,
    optim = fit,
    method = method,
    control = ocontrol,
    start = start,
    weights = if(identical(as.vector(weights),rep.int(1L,n))) NULL else weights,
    offset = list(
      count = if (identical(offsetx, rep.int(0, n))) NULL else offsetx,
      zero = if (identical(offsetz, rep.int(0, n))) NULL else offsetz
    ),
    n = nobs,
    df.null = nobs - 2,
    df.residual = nobs - (kx + kz + (dist == "negbin")),
    theta = theta,
    SE.logtheta = SE.logtheta,
    loglik = fit$value,
    vcov = vc,
    dist = dist,
    link = link,
    linkinv = linkinv,
    converged = fit$convergence < 1,
    contrasts = list(
      count = attr(X, "contrasts"), zero = attr(Z, "contrasts"))
    )
}

predict_zi <- function(object, newdata){
  X <- newdata$x
  Z <- newdata$z
  offsetx <- newdata$offset$count
  offsetz <- newdata$offset$zero
  if (is.null(offsetx)) offsetx <- rep.int(0, NROW(X))
  if (is.null(offsetz)) offsetz <- rep.int(0, NROW(Z))
  mu <- exp(X %*% object$coefficients$count + offsetx)[, 1]
  phi <- object$linkinv(Z %*% object$coefficients$zero + offsetz)[, 1]
  list(mu = mu, phi = phi, y_hat = (1 - phi) * mu)
}

