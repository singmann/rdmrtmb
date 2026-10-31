#' rdmrtmb: Frequentist Estimation of the Racing Diffusion Model
#'
#' Maximum likelihood estimation of the racing diffusion model of speeded
#' decision making. Every model parameter is described by its own formula,
#' collected with [rdm_formula()]. Random effects are written inside these
#' formulas with lme4-style terms such as `(1 | id)` or `(1 + Condition | id)`
#' and are integrated with the Laplace approximation using automatic
#' differentiation provided by RTMB.
#'
#' The drift rates of a race model necessarily differ between the accumulators:
#' each response alternative accumulates evidence at its own rate. The remaining
#' parameters (the start point range `A`, the threshold `b`, the non-decision
#' time `t0` and the diffusion constant `s`) are properties of the race as a
#' whole and are shared by the accumulators unless a term is placed inside
#' `acc()`. `b = ~ acc(1)` gives one threshold per accumulator, which is a
#' response bias, and `v = ~ Condition + acc(.acc == Word)` adds one coefficient
#' that compares the accumulator with a trial variable; see [acc()].
#'
#' @keywords internal
#' @importFrom stats model.matrix as.formula reformulate setNames
#' @importFrom stats coef vcov confint logLik nobs fitted predict simulate
#' @importFrom stats qnorm rnorm runif
#' @importFrom utils modifyList
#' @importFrom RTMB dnorm pnorm matrix rowSums colSums diag
#' @importMethodsFrom RTMB dnorm pnorm matrix rowSums colSums diag
"_PACKAGE"
