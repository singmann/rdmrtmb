#' rdmrtmb: Frequentist Estimation of the Racing Diffusion Model
#'
#' Maximum likelihood estimation of the racing diffusion model of speeded
#' decision making. Every model parameter can be described by its own model
#' matrix, and random effects can be specified per parameter with lme4-style
#' terms such as `(1 | id)`. Random effects are integrated with the Laplace
#' approximation using automatic differentiation provided by RTMB.
#'
#' @keywords internal
#' @importFrom stats model.matrix as.formula setNames
#' @importFrom stats coef vcov confint logLik nobs fitted predict simulate
#' @importFrom stats qnorm rnorm runif
#' @importFrom utils modifyList
#' @importFrom RTMB dnorm pnorm matrix rowSums colSums diag
#' @importMethodsFrom RTMB dnorm pnorm matrix rowSums colSums diag
"_PACKAGE"
