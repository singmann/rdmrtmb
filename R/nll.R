# The negative log-likelihood of the RDM as a function that can be taped by
# RTMB, together with starting values and the bookkeeping that maps the
# parameter vector onto the model parameters.

#' Build the negative log-likelihood of an RDM model
#'
#' Combines the design matrices of `rdm_design()` and the random-effect terms of
#' `rdm_ranef()` into an objective function suitable for `RTMB::MakeADFun()`.
#'
#' The model parameters are linked to the linear predictors with a log link, and
#' the threshold is modelled as `b = A + distance`, so that `b > A` holds by
#' construction. Effects on `b` therefore describe the log threshold distance.
#' A smooth barrier keeps the decision time `rt - t0` positive; it is zero for
#' admissible values and grows linearly outside them.
#'
#' A free diffusion constant `s` scales `A`, `b` and `v`, so it is only
#' identified when one of them is fixed; `rdm_design()` checks this.
#'
#' @param design An `rdm_design` object.
#' @param ranef An `rdm_ranef` object, or `NULL` for a model without random
#'   effects.
#' @param prior A list with elements `sd` and `beta` giving the scale of the
#'   half-normal prior on the random-effect standard deviations and of the
#'   normal prior on the fixed effects. Use `NULL` elements for flat priors.
#'   Weak priors are needed in practice: without them the Laplace
#'   approximation is unbounded for large variance components.
#' @param barrier Weight of the smooth barrier for inadmissible decision times.
#' @return A list with the objective `fn`, the starting values `par`, the name
#'   of the random-effects parameter `random` (or `NULL`), and `pars()`, which
#'   returns the model parameters for a given parameter vector.
#' @keywords internal
#' @noRd
rdm_nll <- function(design, ranef = NULL,
                    prior = list(sd = 0.5, beta = 2), barrier = 1e3) {
  n <- design$n
  K <- design$K

  ## coefficient layout
  present <- rdm_pars[!vapply(rdm_pars, \(nm) is.null(design[[nm]]), TRUE)]
  beta_names <- character()
  idx <- list()
  for (nm in present) {
    spec <- design[[nm]]
    if (!spec$free) {
      idx[[nm]] <- integer()
      next
    }
    idx[[nm]] <- length(beta_names) + seq_along(spec$coef_names)
    beta_names <- c(beta_names, spec$coef_names)
  }

  ## random-effect layout: one contiguous block of parameters per term
  blocks <- if (is.null(ranef)) list() else ranef$blocks
  u_idx <- vector("list", length(blocks))
  theta_idx <- vector("list", length(blocks))
  theta_names <- character()
  n_u <- 0L
  n_theta <- 0L
  for (b in seq_along(blocks)) {
    if (!design[[blocks[[b]]$par]]$free) {
      stop("random effects are not allowed for the fixed parameter `",
           blocks[[b]]$par, "`.", call. = FALSE)
    }
    u_idx[[b]] <- n_u + seq_len(blocks[[b]]$n_lev * blocks[[b]]$q)
    n_u <- n_u + blocks[[b]]$n_lev * blocks[[b]]$q
    theta_idx[[b]] <- n_theta + seq_len(blocks[[b]]$n_theta)
    n_theta <- n_theta + blocks[[b]]$n_theta
    theta_names <- c(theta_names,
                     paste0(blocks[[b]]$label, ".", seq_len(blocks[[b]]$n_theta)))
  }

  ## starting values
  beta <- setNames(rep(0, length(beta_names)), beta_names)
  for (nm in present) {
    if (!design[[nm]]$free) next
    beta[idx[[nm]]] <- coef_starts(design[[nm]]$coef_names, nm, design)
  }
  theta <- setNames(rep(0, length(theta_names)), theta_names)
  for (b in seq_along(blocks)) {
    theta[theta_idx[[b]][seq_len(blocks[[b]]$q)]] <- log(0.1)
  }
  par <- list(beta = beta)
  if (length(theta)) par$theta <- theta
  if (length(blocks)) par$u <- setNames(rep(0, n_u), paste0("u", seq_len(n_u)))

  diagonal <- unlist(lapply(seq_along(blocks),
                            \(b) theta_idx[[b]][seq_len(blocks[[b]]$q)]))

  ## linear predictors for a given parameter vector
  lin_fun <- function(beta, u) {
    lin <- vector("list", length(present))
    names(lin) <- present
    for (nm in present) {
      spec <- design[[nm]]
      val <- if (spec$free) drop(spec$X %*% beta[idx[[nm]]]) else spec$value
      for (b in seq_along(blocks)) {
        if (identical(blocks[[b]]$par, nm)) {
          val <- val + rdm_ranef_contrib(blocks[[b]], u[u_idx[[b]]])
        }
      }
      lin[[nm]] <- val
    }
    lin
  }

  ## model parameters implied by a parameter vector
  pars_fun <- function(beta, u = NULL) {
    lin <- lin_fun(beta, if (is.null(u)) par$u else u)
    ## free parameters are exponentiated, fixed ones are used as given
    A <- if (design$A$free) exp(lin$A) else lin$A
    b <- if (design$b$free) A + exp(lin$b) else lin$b   # distance, so b > A
    t0 <- if (design$t0$free) exp(lin$t0) else lin$t0
    s <- if (is.null(lin$s)) rep(1, n) else if (design$s$free) exp(lin$s) else lin$s
    v <- if (design$v$free) exp(lin$v) else lin$v
    list(A = A, b = b, t0 = t0, s = s,
         v = matrix(v / rep(s, times = K), nrow = n, ncol = K))
  }

  fn <- function(par) {
    mp <- pars_fun(par$beta, par$u)
    t <- design$rt - mp$t0

    ## barrier keeping the decision time positive
    nll <- barrier * sum(smooth_violation(1e3 * t))
    t_use <- sqrt(t^2 + 1e-12)
    nll <- nll - sum(rdm_loglik(t_use, mp$A / mp$s, mp$b / mp$s, mp$v,
                                design$response))

    if (!is.null(prior$beta)) {
      nll <- nll - sum(dnorm(par$beta, 0, prior$beta, log = TRUE))
    }
    for (b in seq_along(blocks)) {
      nll <- nll + rdm_ranef_penalty(blocks[[b]], par$u[u_idx[[b]]],
                                     par$theta[theta_idx[[b]]])
    }
    if (!is.null(prior$sd) && length(diagonal)) {
      sd <- exp(par$theta[diagonal])
      nll <- nll - sum(dnorm(sd, 0, prior$sd, log = TRUE) + log(sd))
    }
    nll
  }

  list(fn = fn, par = par,
       random = if (length(blocks)) "u" else NULL,
       pars = pars_fun,
       beta_names = beta_names, theta_names = theta_names,
       index = list(beta = idx, theta = theta_idx, u = u_idx),
       design = design, ranef = ranef, present = present)
}

## a smooth version of max(0, x), exactly zero for positive x
smooth_violation <- function(x) {
  0.5 * (sqrt(x^2 + 1e-6) - x)
}

## starting values for the coefficients of one parameter
coef_starts <- function(names, par, design) {
  out <- rep(0, length(names))
  intercept <- grepl("\\(Intercept\\)", names)
  default <- switch(par,
    A = 0.5,
    b = 1,
    t0 = 0.5 * min(design$rt),        # must stay below all response times
    v = 2,
    s = 1
  )
  out[intercept] <- log(default)
  ## drift rates: accumulators with more responses start with a larger rate
  if (par == "v" && any(intercept)) {
    prop <- pmax(tabulate(design$response, design$K) / design$n, 0.5 / design$n)
    level <- sub("^.*\\[(.*)\\].*$", "\\1", names[intercept])
    hit <- match(level, design$response_levels)
    out[intercept] <- ifelse(is.na(hit), log(2), log(2 * prop[hit] / mean(prop)))
  }
  out
}
