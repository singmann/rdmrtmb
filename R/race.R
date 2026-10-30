# Distributional quantities of the fitted race: choice probabilities, predicted
# response-time quantiles and random data generation. These are plain numeric
# computations on the trial-level parameters stored in a fit.
#
# All parameters are stored as matrices with one row per trial and one column
# per accumulator, so that every accumulator can have its own values.

## The diffusion constant scales all distances and rates, so the distribution of
## the model with diffusion constant `s` is that of the model with `s = 1` and
## parameters A/s, b/s and v/s.
rdm_scale <- function(mp) {
  s <- if (is.null(mp$s)) 1 else mp$s
  list(A = mp$A / s, b = mp$b / s, v = mp$v / s, t0 = mp$t0)
}

## Integration grid for decision times. The substitution t = tau * u / (1 - u)
## maps the positive half line to the unit interval, so a moderate number of
## points covers both the bulk and the exponentially decaying tail.
rdm_time_grid <- function(tau, m) {
  u <- seq_len(m) / (m + 1)
  list(t = outer(tau, u / (1 - u)),
       w = outer(tau, (1 / (m + 1)) / (1 - u)^2))
}

## Choice probabilities and, optionally, response-time quantiles conditional on
## each accumulator winning.
rdm_race_dist <- function(mp, m = 300L, probs = NULL) {
  sc <- rdm_scale(mp)
  n <- nrow(mp$v)
  K <- ncol(mp$v)
  A <- as_par_matrix(sc$A, n, K)
  b <- as_par_matrix(sc$b, n, K)
  v <- as_par_matrix(sc$v, n, K)
  t0 <- as_par_matrix(sc$t0, n, K)

  ## all accumulators of a trial are integrated on the same grid, which is
  ## scaled by the mean decision time of the fastest one
  mdt <- (b - 0.5 * A) / v
  gr <- rdm_time_grid(apply(mdt, 1, min), m)
  weight <- gr$w
  grid <- gr$t

  log_f <- log_S <- vector("list", K)
  for (j in seq_len(K)) {
    Aj <- matrix(A[, j], n, m)
    bj <- matrix(b[, j], n, m)
    vj <- matrix(v[, j], n, m)
    log_f[[j]] <- rdm_log_f(gr$t, Aj, bj, vj)
    log_S[[j]] <- rdm_log_S(gr$t, Aj, bj, vj)
  }
  total_S <- Reduce(`+`, log_S)

  prob <- matrix(0, n, K, dimnames = list(NULL, colnames(mp$v)))
  q <- NULL
  if (!is.null(probs)) {
    q <- array(NA_real_, c(n, K, length(probs)),
               dimnames = list(NULL, colnames(mp$v), NULL))
  }
  for (j in seq_len(K)) {
    mass <- exp(log_f[[j]] + total_S - log_S[[j]]) * weight
    prob[, j] <- rowSums(mass)
    if (!is.null(q)) {
      for (i in seq_along(probs)) {
        q[, j, i] <- t0[, j] + rdm_inv_cdf(mass, grid, probs[i])
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
  sc <- rdm_scale(mp)
  n <- nrow(mp$v)
  K <- ncol(mp$v)
  A <- as_par_matrix(sc$A, n, K)
  b <- as_par_matrix(sc$b, n, K)
  v <- as_par_matrix(sc$v, n, K)
  t0 <- as_par_matrix(sc$t0, n, K)

  start <- matrix(stats::runif(n * K), n, K) * A
  distance <- b - start

  fpt <- matrix(0, n, K)
  for (j in seq_len(K)) {
    ## inverse Gaussian with mean distance / v and shape distance^2
    d <- matrix(distance[, j])
    fpt[, j] <- rdm_rig(d / matrix(v[, j]), d^2)
  }
  win <- max.col(-fpt, ties.method = "first")
  pick <- cbind(seq_len(n), win)
  list(rt = t0[pick] + fpt[pick], response = win)
}

## Inverse-Gaussian draws with mean `mu` and shape `lambda`, using the
## Michael-Schucany-Haas transformation of a squared normal variate.
rdm_rig <- function(mu, lambda) {
  y <- stats::rnorm(length(mu))^2
  x <- mu + mu^2 * y / (2 * lambda) -
    mu / (2 * lambda) * sqrt(4 * mu * lambda * y + mu^2 * y^2)
  ifelse(stats::runif(length(mu)) < mu / (mu + x), x, mu^2 / x)
}
