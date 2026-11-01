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
#' start point. Use `coef(fit, type = "ratio")` for the exponentiated
#' coefficients, which are multiplicative effects rather than the parameters,
#' and [fitted.rdm_fit()] for the trial-level parameters themselves.
#'
#' A term inside `acc()` is estimated separately for each accumulator, so
#' `b = ~ acc(1)` gives one threshold distance per accumulator (a response
#' bias). Inside `acc()` the current accumulator is available as `.acc`, and a
#' term that refers to it is a single predictor instead: `acc(.acc == Word)` is
#' one coefficient that is 1 where the accumulator matches the word. See
#' [acc()].
#'
#' The left-hand side of a random-effect term is built like the fixed part, so
#' `acc()` can be used there as well and the random effects can mirror the fixed
#' effects: with `v = ~ acc(Condition) + (acc(Condition) | id)` every
#' accumulator-and-condition effect of the fixed part has a random counterpart.
#' The accumulator variables `response` and `.acc` are available there too.
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
#'   the formulas as in lme4, e.g. `A = ~ 1 + (1 | id)`. A term can be wrapped
#'   in `diag()`, `homdiag()` or `cs()` to choose its covariance structure, as
#'   in glmmTMB, e.g. `v = ~ acc(1) + diag(0 + Condition | id)`. See [acc()]
#'   and [rdm_formula()].
#' @param fixed Named list of parameters held constant on their natural scale,
#'   e.g. `list(s = 1)`. A fixed `b` is the threshold itself, not the distance.
#'   The value can have length 1 (shared) or one value per accumulator. `A` may
#'   also be fixed at zero for all accumulators, which drops the start-point
#'   variability; that is the same model as leaving `A` out of the formula.
#' @param contrasts Named list of contrast matrices, passed on to
#'   [stats::model.matrix()]. The contrasts of a term inside `acc()` are given
#'   under its full name, e.g. `list("acc(.acc == Word)" = contr.sum)`.
#' @param anchor Named list choosing which accumulator the diffusion constant
#'   `s` is anchored at, e.g. `list(s = "Word")`. A free `s` scales `A`, `b` and
#'   `v`, so the scale of the model has to be pinned down somewhere. With
#'   `s = ~ acc(1)` the diffusion constant of the first response level is fixed
#'   at 1 and the other ones are estimated relative to it; `anchor` picks a
#'   different accumulator. Anchoring one accumulator is only enough when at
#'   least one of `A`, `b` and `v` is shared by all accumulators.
#' @param prior A list with elements `sd` and `beta` giving the scale of the
#'   half-normal prior on the random-effect standard deviations and of the
#'   normal prior on the fixed effects (and on the correlations of `cs()`
#'   blocks). Use `NULL` for a flat prior. Weak
#'   priors are needed in practice: with a flat prior the Laplace
#'   approximation is unbounded for large variance components.
#' @param barrier Weight of the smooth penalty for response times shorter than
#'   `t0`. Large values enforce the constraint more strictly.
#' @param control List of control arguments passed to [stats::nlminb()].
#' @param inner.control List of control arguments for the inner optimisation of
#'   the random effects.
#' @param start Starting values. Either a named vector for the fixed effects, on
#'   the log scale, e.g. `c("t0:(Intercept)" = log(0.3))`, or a previous
#'   `rdm_fit`, or a list with elements `beta`, `theta` and `u`. Restarting from
#'   a previous fit is much faster than fitting from scratch, which matters for
#'   larger models where every Laplace evaluation costs seconds. Coefficients in
#'   a vector that the model does not have are an error; coefficients the vector
#'   leaves out keep their default start.
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
#'     \item{`inner_eig`}{the smallest and largest eigenvalue of the Hessian of
#'       the joint objective with respect to the random effects, at the
#'       conditional modes; `NULL` for a model without random effects, or with
#'       more than 2000 of them. A non-positive smallest eigenvalue, or a
#'       smallest / largest ratio below 1e-8, triggers a warning: the
#'       conditional modes are then not a proper mode of the integrand and the
#'       Laplace value and its standard errors cannot be trusted},
#'     \item{`outer_eig`}{the smallest and largest eigenvalue of the Hessian of
#'       the Laplace-approximated objective with respect to the fixed effects,
#'       from the covariance that `sdreport()` returns; `NULL` without
#'       `sdreport`, or when that covariance is not positive definite. A
#'       smallest / largest ratio below 1e-8 triggers a warning, because some
#'       combination of coefficients is then barely identified},
#'     \item{`design`, `ranef`, `nll`, `obj`}{the underlying model pieces},
#'   }
#'   and the attributes needed by the S3 methods.
#' @seealso [rdm_formula()], [acc()], [coef.rdm_fit()], [confint.rdm_fit()],
#'   [ranef.rdm_fit()], [predict.rdm_fit()], [simulate.rdm_fit()]
#' @export
rdm_fit <- function(data, rt, response, formula, fixed = list(s = 1),
                    prior = list(sd = 0.5, beta = 2),
                    barrier = 1e3, control = list(), inner.control = list(),
                    start = NULL, sdreport = TRUE, contrasts = list(),
                    anchor = list()) {
  formula <- as_rdm_formula(formula)
  design <- rdm_design(data, rt, response, formula, fixed, contrasts, anchor)
  ranef <- rdm_ranef(design)
  nll <- rdm_nll(design, ranef, prior = prior, barrier = barrier)

  nll$par <- start_values(start, nll)

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

  ## the inner Hessian: if it is indefinite or nearly singular the conditional
  ## modes are not a proper mode and the Laplace value is unreliable
  inner_eig <- rdm_inner_eigen(obj)
  if (!is.null(inner_eig)) {
    ratio <- inner_eig[["min"]] / inner_eig[["max"]]
    if (inner_eig[["min"]] <= 0) {
      warning("the inner Hessian for the random effects is not positive definite ",
              "(smallest eigenvalue ", format(inner_eig[["min"]], digits = 3),
              "); the Laplace approximation is unreliable.", call. = FALSE)
    } else if (ratio < 1e-8) {
      warning("the inner Hessian for the random effects is close to singular ",
              "(smallest / largest eigenvalue ", format(ratio, digits = 3),
              "); the Laplace approximation is unreliable.", call. = FALSE)
    }
  }

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

  ## the outer Hessian: a large ratio between its largest and smallest
  ## eigenvalue means that some combination of coefficients is barely identified
  outer_eig <- rdm_outer_eigen(sdr)
  if (!is.null(outer_eig) &&
      outer_eig[["min"]] / outer_eig[["max"]] < 1e-8) {
    warning("the fixed effects are close to collinear ",
            "(smallest / largest eigenvalue of the outer Hessian ",
            format(outer_eig[["min"]] / outer_eig[["max"]], digits = 3),
            "); their standard errors are unreliable.", call. = FALSE)
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
    inner_eig = inner_eig,
    outer_eig = outer_eig,
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

## The extreme eigenvalues of the Hessian of the joint objective with respect to
## the random effects, at the conditional modes. A nearly singular or indefinite
## matrix means the modes are not a proper mode of the integrand, which makes the
## Laplace value path dependent and its standard errors meaningless.
## `obj` must have been evaluated at the optimum for the modes to be current.
rdm_inner_eigen <- function(obj, n_max = 2000) {
  n_u <- length(obj$env$random)
  if (n_u == 0L || n_u > n_max) return(NULL)
  ev <- tryCatch({
    H <- obj$env$spHess(obj$env$last.par, random = TRUE)
    eigen(as.matrix(H), symmetric = TRUE, only.values = TRUE)$values
  }, error = function(e) NULL)
  if (is.null(ev) || !length(ev)) return(NULL)
  c(min = min(ev), max = max(ev))
}

## The extreme eigenvalues of the Hessian of the Laplace-approximated objective
## with respect to the fixed effects, taken from the fixed-effect covariance
## that `sdreport()` returns. A large ratio between the largest and the smallest
## eigenvalue means that some combination of coefficients is barely identified.
rdm_outer_eigen <- function(sdr, n_max = 500) {
  V <- if (is.null(sdr)) NULL else sdr$cov.fixed
  if (is.null(V) || !length(V) || nrow(V) == 0L || nrow(V) > n_max) return(NULL)
  ev <- tryCatch(eigen(as.matrix(V), symmetric = TRUE, only.values = TRUE)$values,
                 error = function(e) NULL)
  if (is.null(ev) || any(!is.finite(ev)) || any(ev <= 0)) return(NULL)
  ev <- 1 / ev
  c(min = min(ev), max = max(ev))
}

## starting values from a vector of fixed effects, a previous fit, or a list with
## elements `beta`, `theta` and `u`
start_values <- function(start, nll) {
  par <- nll$par
  if (is.null(start)) return(par)
  if (inherits(start, "rdm_fit")) {
    start <- list(beta = start$beta, theta = start$theta, u = start$u)
  }
  if (!is.list(start)) start <- list(beta = start)
  unknown <- setdiff(names(start), names(par))
  if (length(unknown)) {
    stop("unknown element(s) in `start`: ", paste(unknown, collapse = ", "),
         ".", call. = FALSE)
  }
  for (nm in names(start)) {
    x <- start[[nm]]
    if (is.null(x) || is.null(par[[nm]])) next
    named <- !is.null(names(x)) && all(names(x) != "")
    ## the fixed effects may be given by name, one or a few at a time; the
    ## covariance parameters and the random effects have to be complete
    if (named && nm == "beta") {
      found <- names(x) %in% names(par[[nm]])
      if (!all(found)) {
        stop("unknown name(s) in `start$", nm, "`: ",
             paste(names(x)[!found], collapse = ", "), ".", call. = FALSE)
      }
      par[[nm]][names(x)] <- x
      next
    }
    if (length(x) != length(par[[nm]])) {
      stop("`start$", nm, "` has ", length(x), " values but the model has ",
           length(par[[nm]]), ".", call. = FALSE)
    }
    par[[nm]][] <- x
  }
  par
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
