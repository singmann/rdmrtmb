# Likelihood of the racing diffusion model (RDM).
#
# The functions in this file are written so that they can be taped by RTMB:
# only vectorised operations supporting automatic differentiation are used.
# They are internal; the user-facing interface is fit_rdm().

#' Log-density of a single accumulator of the RDM
#'
#' Density of the first-passage time of one accumulator with threshold `b`,
#' drift rate `v`, and a start point uniformly distributed on `[0, A]`, with the
#' within-trial diffusion constant fixed at `s = 1`. Obtained by averaging the
#' inverse-Gaussian first-passage density over the start point.
#'
#' @param t Decision time, i.e. the response time minus `t0`. Must be positive.
#' @param A Upper bound of the uniformly distributed start point.
#' @param b Response threshold, must be larger than `A`.
#' @param v Drift rate, must be positive.
#' @return The log density, a vector with the length of `t`.
#' @keywords internal
#' @noRd
rdm_log_f <- function(t, A, b, v) {
  st <- sqrt(t)
  mv <- v * t
  zb <- (b - mv) / st
  zba <- (b - A - mv) / st
  dens <- (mv * (pnorm(zb) - pnorm(zba)) +
             st * (dnorm(zba) - dnorm(zb))) / (A * t)
  log(dens + 1e-300)
}

#' Log survival probability of a single accumulator of the RDM
#'
#' Probability that one accumulator has not reached its threshold by decision
#' time `t`, i.e. `P(T > t)`. Same parameter interpretation as `rdm_log_f()`.
#'
#' @param t Decision time, i.e. the response time minus `t0`. Must be positive.
#' @param A Upper bound of the uniformly distributed start point.
#' @param b Response threshold, must be larger than `A`.
#' @param v Drift rate, must be positive.
#' @return The log survival probability, a vector with the length of `t`.
#' @keywords internal
#' @noRd
rdm_log_S <- function(t, A, b, v) {
  st <- sqrt(t)
  vt <- v * t
  zeta <- function(x) (x + vt) / st     # z of the exponential correction
  # Phi(zeta(k) - 2 v sqrt(t)) is Phi((k - vt)/sqrt(t)), which I1 needs as well
  eta_b <- (b - vt) / st
  eta_ba <- (b - A - vt) / st
  pn_b <- pnorm(eta_b)
  pn_ba <- pnorm(eta_ba)
  # integral of Phi((k - vt)/sqrt(t)) over the start points k in [b - A, b]
  I1 <- st * (eta_b * pn_b + dnorm(eta_b) - eta_ba * pn_ba - dnorm(eta_ba))
  # exp(2vk) * Phi(-zeta(k)) evaluated on the log scale to avoid overflow
  I2 <- (exp(2 * v * b + pnorm(-zeta(b), log.p = TRUE)) -
           exp(2 * v * (b - A) + pnorm(-zeta(b - A), log.p = TRUE))) / (2 * v) +
    (pn_b - pn_ba) / (2 * v)
  log((I1 - I2) / A + 1e-300)
}

#' Log-density of a single accumulator without start-point variability
#'
#' With `A = 0` every accumulator starts at zero, so its first-passage time is
#' inverse Gaussian with mean `b / v` and shape `b^2`. This is the limit of
#' `rdm_log_f()` for `A` going to zero, written out because that function
#' divides by `A`.
#'
#' @param t Decision time, i.e. the response time minus `t0`. Must be positive.
#' @param b Response threshold, must be positive.
#' @param v Drift rate, must be positive.
#' @return The log density, a vector with the length of `t`.
#' @keywords internal
#' @noRd
rdm_log_f_ig <- function(t, b, v) {
  z <- (b - v * t) / sqrt(t)
  log(b) - 0.5 * log(2 * pi) - 1.5 * log(t) - 0.5 * z^2
}

#' Log survival probability of a single accumulator without start-point
#' variability
#'
#' `P(T > t)` for the inverse-Gaussian first-passage time of `rdm_log_f_ig()`.
#' The second term of the survival function is computed on the log scale,
#' because `exp(2 v b)` overflows for larger drift rates and thresholds.
#'
#' @inheritParams rdm_log_f_ig
#' @return The log survival probability, a vector with the length of `t`.
#' @keywords internal
#' @noRd
rdm_log_S_ig <- function(t, b, v) {
  st <- sqrt(t)
  ## survival 1 - F(t) = Phi(z) - exp(2vb) Phi(-(b + vt) / sqrt(t)); the terms
  ## are close for large t, so both are kept on the log scale and subtracted
  ## there, which is also what keeps exp(2vb) from overflowing
  ## pnorm() is called without the namespace, as RTMB can only tape the bare
  ## function
  log_first <- pnorm((b - v * t) / st, log.p = TRUE)
  log_second <- 2 * v * b + pnorm(-(b + v * t) / st, log.p = TRUE)
  RTMB::logspace_sub(log_first, log_second)
}

#' Log-likelihood of one trial of the RDM
#'
#' The accumulator that finishes first determines both the response and the
#' response time, so the likelihood of a trial with response `j` is the density
#' of accumulator `j` at `t` times the survival probabilities of all other
#' accumulators at `t`.
#'
#' @param t Decision time, a matrix with one row per trial and one column per
#'   accumulator, i.e. the response time minus the non-decision time.
#' @param A,b Matrices of start point and threshold, with the same dimensions as
#'   `t`; see `rdm_log_f()`. With `fixed_start = TRUE`, `A` is ignored.
#' @param v Matrix of drift rates with one row per trial and one column per
#'   accumulator.
#' @param response Integer vector of observed responses, coded `1:K`.
#' @param fixed_start If `TRUE`, all accumulators start at zero (`A = 0`), and the
#'   inverse-Gaussian functions are used instead of the start-averaged ones,
#'   which divide by `A`.
#' @return A vector of log likelihoods, one per trial.
#' @keywords internal
#' @noRd
rdm_loglik <- function(t, A, b, v, response, fixed_start = FALSE) {
  K <- ncol(v)
  n <- nrow(v)
  t <- as_par_matrix(t, n, K)
  b <- as_par_matrix(b, n, K)
  win <- cbind(seq_len(n), response)
  ## the density is only needed for the accumulator that won; computing all the
  ## survival probabilities and subtracting the winner's is cheaper than
  ## subsetting the matrix three times
  ## not vapply(): it rejects the representation RTMB uses when taping second
  ## derivatives for the Laplace approximation
  log_S <- if (fixed_start) {
    do.call(cbind, lapply(seq_len(K),
                          \(j) rdm_log_S_ig(t[, j], b[, j], v[, j])))
  } else {
    A <- as_par_matrix(A, n, K)
    do.call(cbind, lapply(seq_len(K),
                          \(j) rdm_log_S(t[, j], A[, j], b[, j], v[, j])))
  }
  log_f <- if (fixed_start) {
    rdm_log_f_ig(t[win], b[win], v[win])
  } else {
    rdm_log_f(t[win], A[win], b[win], v[win])
  }
  log_f + rowSums(log_S) - log_S[win]
}

## parameters that are the same for all accumulators are given as a vector of
## length n and are recycled to one value per accumulator
as_par_matrix <- function(x, n, K) {
  if (!is.null(dim(x))) x else matrix(x, n, K)
}
