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
  eta <- function(x) (x - vt) / st      # z of the first-passage CDF
  zeta <- function(x) (x + vt) / st     # z of the exponential correction
  # integral of Phi((k - vt)/sqrt(t)) over the start points k in [b - A, b]
  I1 <- st * (eta(b) * pnorm(eta(b)) + dnorm(eta(b)) -
                eta(b - A) * pnorm(eta(b - A)) - dnorm(eta(b - A)))
  # exp(2vk) * Phi(-zeta(k)) evaluated on the log scale to avoid overflow
  I2 <- (exp(2 * v * b + pnorm(-zeta(b), log.p = TRUE)) -
           exp(2 * v * (b - A) + pnorm(-zeta(b - A), log.p = TRUE))) / (2 * v) +
    (pnorm(zeta(b) - 2 * v * st) - pnorm(zeta(b - A) - 2 * v * st)) / (2 * v)
  log((I1 - I2) / A + 1e-300)
}

#' Log-likelihood of one trial of the RDM
#'
#' The accumulator that finishes first determines both the response and the
#' response time, so the likelihood of a trial with response `j` is the density
#' of accumulator `j` at `t` times the survival probabilities of all other
#' accumulators at `t`.
#'
#' @param t Decision time, a vector with one element per trial.
#' @param A,b,t0 Trial-wise parameter vectors, see `rdm_log_f()`.
#' @param v Matrix of drift rates with one row per trial and one column per
#'   accumulator.
#' @param response Integer vector of observed responses, coded `1:K`.
#' @return A vector of log likelihoods, one per trial.
#' @keywords internal
#' @noRd
rdm_loglik <- function(t, A, b, v, response) {
  K <- ncol(v)
  n <- length(t)
  ## not vapply(): it rejects the representation RTMB uses when taping second
  ## derivatives for the Laplace approximation
  log_f <- do.call(cbind, lapply(seq_len(K), \(j) rdm_log_f(t, A, b, v[, j])))
  log_S <- do.call(cbind, lapply(seq_len(K), \(j) rdm_log_S(t, A, b, v[, j])))
  out <- numeric(n)
  for (j in seq_len(K)) {
    out <- out + (response == j) * (log_f[, j] + rowSums(log_S[, -j, drop = FALSE]))
  }
  out
}
