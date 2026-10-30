#' rdmrtmb: Frequentist Estimation of the Racing Diffusion Model
#'
#' Maximum likelihood estimation of the racing diffusion model of speeded
#' decision making. Every model parameter can be described by its own model
#' matrix, and random effects can be specified per parameter with lme4-style
#' terms such as `(1 | id)`. Random effects are integrated with the Laplace
#' approximation using automatic differentiation provided by RTMB.
#'
#' The drift rates of a race model necessarily differ between the accumulators:
#' each response alternative accumulates evidence at its own rate. The remaining
#' parameters (the start point range `A`, the threshold `b`, the non-decision
#' time `t0` and the diffusion constant `s`) are properties of the race as a
#' whole and are shared by the accumulators when they are given a single
#' formula. Any of them can nevertheless be estimated separately for each
#' accumulator by passing a list of formulas, one per accumulator, instead of a
#' single formula; see [fit_rdm()]. A threshold that differs between the
#' response alternatives is a response bias.
#'
#' @keywords internal
#' @importFrom stats model.matrix as.formula setNames
#' @importFrom stats coef vcov confint logLik nobs fitted predict simulate
#' @importFrom stats qnorm rnorm runif
#' @importFrom utils modifyList
#' @importFrom RTMB dnorm pnorm matrix rowSums colSums diag
#' @importMethodsFrom RTMB dnorm pnorm matrix rowSums colSums diag
"_PACKAGE"
