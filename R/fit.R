# The user-facing estimation function. It assembles the design, the random
# effects and the likelihood, runs the outer optimisation over the fixed
# effects (the random effects are profiled out by RTMB's Laplace
# approximation), and collects everything a fit object needs.

#' Fit the racing diffusion model
#'
#' Maximum likelihood estimation of the racing diffusion model (RDM) with the
#' marginal likelihood approximated by the Laplace method. Random effects are
#' profiled out by `RTMB::MakeADFun()` and the fixed effects are optimised with
#' `stats::nlminb()`.
#'
#' All free parameters enter the model through a log link, so the coefficients
#' returned by [coef()] are on the log scale. The threshold is parameterised as
#' `b = A + exp(eta)`, which guarantees `b > A`; its coefficients therefore
#' describe the log distance between the threshold and the upper bound of the
#' start point. Use `coef(fit, type = "response")` for the back-transformed
#' values.
#'
#' A term inside `acc()` is estimated separately for each accumulator, so
#' `b = ~ acc(1)` gives one threshold distance per accumulator (a response
#' bias). Inside `acc()` the current accumulator is available as `.acc`, and a
#' term that refers to it is a single predictor instead: `acc(.acc == Word)` is
#' one coefficient that is 1 where the accumulator matches the word. See
#' [acc()].
#'
#' @param data A `data.frame` with one row per trial.
#' @param rt,response Response times and responses. Either a vector or the name
#'   of a column of `data`.
#' @param formula An [rdm_formula()] object, or a named list of formulas, one
#'   per free model parameter. Valid names are `"A"`, `"b"`, `"t0"`, `"v"` and
#'   `"s"`. Terms shared by all accumulators are written as usual; terms inside
#'   `acc()` are estimated separately for each accumulator, and `acc()` terms
#'   that refer to the accumulator `.acc` compare it with a trial variable, e.g.
#'   `v = ~ Condition + acc(.acc == Word)`. Random effects are written inside
#'   the formulas as in lme4, e.g. `A = ~ 1 + (1 | id)`. See [acc()] and
#'   [rdm_formula()].
#' @param fixed Named list of parameters held constant on their natural scale,
#'   e.g. `list(s = 1)`. A fixed `b` is the threshold itself, not the distance.
#'   The value can have length 1 (shared) or one value per accumulator.
#' @param contrasts Named list of contrast matrices, passed on to
#'   [stats::model.matrix()]. The contrasts of a term inside `acc()` are given
#'   under its full name, e.g. `list("acc(.acc == Word)" = contr.sum)`.
#' @param prior A list with elements `sd` and `beta` giving the scale of the
#'   half-normal prior on the random-effect standard deviations and of the
#'   normal prior on the fixed effects. Use `NULL` for a flat prior. Weak
#'   priors are needed in practice: with a flat prior the Laplace
#'   approximation is unbounded for large variance components.
#' @param barrier Weight of the smooth penalty for response times shorter than
#'   `t0`. Large values enforce the constraint more strictly.
#' @param control List of control arguments passed to [stats::nlminb()].
#' @param inner.control List of control arguments for the inner optimisation of
#'   the random effects.
#' @param start Optional named vector of starting values for the fixed effects,
#'   on the log scale, e.g. `c("t0:(Intercept)" = log(0.3))`. Names must match
#'   those of the model coefficients.
#' @param sdreport If `TRUE` (the default), standard errors are computed with
#'   `RTMB::sdreport()`.
#' @return An object of class `rdm_fit`, a list with the components
#'   \describe{
#'     \item{`beta`}{the fixed-effect coefficients (log scale)},
#'     \item{`theta`}{the covariance parameters of the random effects (log
#'       standard deviations and Cholesky elements)},
#'     \item{`u`}{the conditional modes of the random effects},
#'     \item{`pars`}{the trial-level model parameters implied by the fit},
#'     \item{`loglik`, `df`, `nobs`}{the Laplace log-likelihood, its degrees of
#'       freedom and the number of trials},
#'     \item{`opt`, `sd`}{the optimiser result and the `sdreport` output},
#'     \item{`design`, `ranef`, `nll`, `obj`}{the underlying model pieces},
#'   }
#'   and the attributes needed by the S3 methods.
#' @seealso [rdm_formula()], [acc()], [coef.rdm_fit()], [confint.rdm_fit()],
#'   [ranef.rdm_fit()], [predict.rdm_fit()], [simulate.rdm_fit()]
#' @export
fit_rdm <- function(data, rt, response, formula, fixed = list(s = 1),
                    prior = list(sd = 0.5, beta = 2),
                    barrier = 1e3, control = list(), inner.control = list(),
                    start = NULL, sdreport = TRUE, contrasts = list()) {
  if (is.character(rt) && length(rt) == 1L) rt <- data[[rt]]
  if (is.character(response) && length(response) == 1L) {
    response <- data[[response]]
  }

  formula <- as_rdm_formula(formula)
  design <- rdm_design(data, rt, response, formula, fixed, contrasts)
  ranef <- rdm_ranef(design)
  nll <- rdm_nll(design, ranef, prior = prior, barrier = barrier)

  if (!is.null(start)) {
    found <- names(start) %in% names(nll$par$beta)
    if (!all(found)) {
      stop("unknown name(s) in `start`: ",
           paste(names(start)[!found], collapse = ", "), ".", call. = FALSE)
    }
    nll$par$beta[names(start)] <- start
  }

  obj <- RTMB::MakeADFun(nll$fn, nll$par, random = nll$random, silent = TRUE,
                         inner.control = inner.control)

  ctrl <- utils::modifyList(list(iter.max = 1000, eval.max = 1000), control)
  opt <- rdm_optimise(obj, ctrl)

  if (opt$convergence != 0) {
    warning("the optimiser reported convergence code ", opt$convergence,
            ": ", opt$message, call. = FALSE)
  }

  full <- obj$env$parList(opt$par)
  mp <- nll$pars(full$beta, full$u)

  ## the log-likelihood itself, without the priors that only serve to make the
  ## optimisation well behaved
  flat <- rdm_nll(design, ranef, prior = list(sd = NULL, beta = NULL),
                  barrier = barrier)
  loglik <- -flat$fn(full)

  grad <- tryCatch(obj$gr(opt$par), error = function(e) rep(NA_real_, length(opt$par)))
  pd_hess <- NA
  sdr <- NULL
  if (isTRUE(sdreport)) {
    sdr <- tryCatch(
      RTMB::sdreport(obj, par.fixed = opt$par, getJointPrecision = TRUE),
      error = function(e) {
        warning("standard errors are not available: ", conditionMessage(e),
                call. = FALSE)
        NULL
      }
    )
    if (!is.null(sdr)) {
      pd_hess <- isTRUE(sdr$pdHess)
      if (!pd_hess) {
        warning("the joint Hessian at the optimum is not positive definite; ",
                "standard errors and Wald intervals are unreliable.",
                call. = FALSE)
      }
    }
  }

  n_theta <- if (length(full$theta)) length(full$theta) else 0L
  structure(list(
    call = match.call(),
    beta = full$beta,
    theta = full$theta,
    u = full$u,
    pars = mp,
    loglik = loglik,
    df = length(full$beta) + n_theta,
    nobs = design$n,
    opt = opt,
    sd = sdr,
    pdHess = pd_hess,
    max_grad = if (all(is.na(grad))) NA_real_ else max(abs(grad)),
    converged = opt$convergence == 0,
    design = design,
    ranef = ranef,
    nll = nll,
    obj = obj,
    prior = prior
  ), class = "rdm_fit")
}

## nlminb, with a BFGS fallback when it does not converge
rdm_optimise <- function(obj, ctrl) {
  ## the inner optimisation of the random effects can produce NaNs while it
  ## searches; those warnings say nothing about the final fit
  quiet <- function(expr) {
    withCallingHandlers(expr, warning = function(w) {
      if (grepl("NA/NaN", conditionMessage(w), fixed = TRUE)) {
        invokeRestart("muffleWarning")
      }
    })
  }
  opt <- quiet(stats::nlminb(obj$par, obj$fn, obj$gr, control = ctrl))
  opt$method <- "nlminb"
  if (opt$convergence != 0) {
    alt <- try(quiet(stats::optim(opt$par, obj$fn, obj$gr, method = "BFGS",
                                  control = list(maxit = ctrl$iter.max))),
               silent = TRUE)
    if (!inherits(alt, "try-error") && is.finite(alt$value) &&
          (!is.finite(opt$objective) || alt$value < opt$objective)) {
      opt <- list(par = alt$par, objective = alt$value,
                  convergence = alt$convergence, message = alt$message,
                  iterations = unname(alt$counts[1]), method = "BFGS")
    }
  }
  opt
}
