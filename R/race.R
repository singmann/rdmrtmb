# Distributional quantities of the fitted race: choice probabilities, predicted
# response-time quantiles and random data generation. These are plain numeric
# computations on the trial-level parameters stored in a fit.

## Integration grid for decision times. The substitution t = tau * u / (1 - u)
## maps the positive half line to the unit interval, so a moderate number of
## points covers both the bulk and the exponentially decaying tail.
rdm_time_grid <- function(mp, m = 300L) {
  tau <- mp$t0 + (mp$b - 0.5 * mp$A) / rowMeans(mp$v)
  u <- seq_len(m) / (m + 1)
  list(u = u,
       t = outer(tau, u / (1 - u)),
       w = outer(tau, (1 / (m + 1)) / (1 - u)^2))
}

## Choice probabilities and, optionally, response-time quantiles conditional on
## each accumulator winning.
rdm_race_dist <- function(mp, m = 300L, probs = NULL) {
  n <- length(mp$A)
  K <- ncol(mp$v)
  gr <- rdm_time_grid(mp, m)
  Amat <- matrix(mp$A, n, m)
  bmat <- matrix(mp$b, n, m)

  log_f <- log_S <- vector("list", K)
  for (j in seq_len(K)) {
    vj <- matrix(mp$v[, j], n, m)
    log_f[[j]] <- rdm_log_f(gr$t, Amat, bmat, vj)
    log_S[[j]] <- rdm_log_S(gr$t, Amat, bmat, vj)
  }
  total_S <- Reduce(`+`, log_S)

  prob <- matrix(0, n, K, dimnames = list(NULL, colnames(mp$v)))
  q <- NULL
  if (!is.null(probs)) {
    q <- array(NA_real_, c(n, K, length(probs)),
               dimnames = list(NULL, colnames(mp$v), NULL))
  }
  for (j in seq_len(K)) {
    mass <- exp(log_f[[j]] + total_S - log_S[[j]]) * gr$w
    prob[, j] <- rowSums(mass)
    if (!is.null(q)) {
      for (i in seq_along(probs)) {
        q[, j, i] <- mp$t0 + rdm_inv_cdf(mass, gr$t, probs[i])
      }
    }
  }
  ## the grid integration is accurate to a few digits; normalising makes the
  ## probabilities sum to one exactly
  prob <- prob / rowSums(prob)
  list(prob = prob, quantiles = q)
}

## Quantile of the distribution whose cell masses are given in `mass`, with cell
## times `t`. The CDF is treated as piecewise linear between cell boundaries.
rdm_inv_cdf <- function(mass, t, p) {
  n <- nrow(mass)
  m <- ncol(mass)
  cdf <- t(apply(mass, 1, cumsum))
  total <- cdf[, m]
  target <- pmin(pmax(p * total, 0), total)
  hit <- max.col(cdf >= target, ties.method = "first")

  hi_lo <- cbind(seq_len(n), hit)
  cdf_hi <- cdf[hi_lo]
  t_hi <- t[hi_lo]
  prev <- pmax(hit - 1L, 1L)
  cdf_lo <- ifelse(hit > 1L, cdf[cbind(seq_len(n), prev)], 0)
  t_lo <- ifelse(hit > 1L, t[cbind(seq_len(n), prev)], 0)

  frac <- (target - cdf_lo) / pmax(cdf_hi - cdf_lo, 1e-300)
  t_lo + pmin(pmax(frac, 0), 1) * (t_hi - t_lo)
}

## One draw of the race per trial. Every accumulator gets a start point uniform
## on [0, A] and an inverse-Gaussian first-passage time; the smallest one wins.
rdm_rand_race <- function(mp) {
  n <- length(mp$A)
  K <- ncol(mp$v)
  start <- matrix(stats::runif(n * K), n, K) * mp$A
  distance <- mp$b - start

  fpt <- matrix(0, n, K)
  for (j in seq_len(K)) {
    ## inverse Gaussian with mean distance / v and shape distance^2
    d <- matrix(distance[, j])
    fpt[, j] <- rdm_rig(d / matrix(mp$v[, j]), d^2)
  }
  win <- max.col(-fpt, ties.method = "first")
  list(rt = mp$t0 + fpt[cbind(seq_len(n), win)], response = win)
}

## Inverse-Gaussian draws with mean `mu` and shape `lambda`, using the
## Michael-Schucany-Haas transformation of a squared normal variate.
rdm_rig <- function(mu, lambda) {
  y <- stats::rnorm(length(mu))^2
  x <- mu + mu^2 * y / (2 * lambda) -
    mu / (2 * lambda) * sqrt(4 * mu * lambda * y + mu^2 * y^2)
  ifelse(stats::runif(length(mu)) < mu / (mu + x), x, mu^2 / x)
}
