# S3 methods for fitted RDM models.

#' Coefficients of a fitted RDM
#'
#' @param object An `rdm_fit` object.
#' @param type `"link"` (default) for the coefficients as estimated, `"ratio"`
#'   for `exp()` of each coefficient. All parameters enter the model through a
#'   log link, so the latter is a multiplicative effect rather than the
#'   parameter: it is the factor by which the parameter changes for a one-unit
#'   change in its covariate, and it equals an actual parameter only where all
#'   of that parameter's other covariates are zero. Use [fitted.rdm_fit()] for
#'   the trial-level parameters. For the threshold `b` the ratio refers to the
#'   distance between the threshold and `A`.
#' @param full If `TRUE`, the covariance parameters of the random effects are
#'   appended.
#' @param ... Unused.
#' @export
coef.rdm_fit <- function(object, type = c("link", "ratio"), full = FALSE,
                         ...) {
  type <- match.arg(type)
  out <- object$beta
  if (full && length(object$theta)) out <- c(out, object$theta)
  if (type == "ratio") out <- exp(out)
  out
}

#' Fixed effects of a fitted model
#'
#' A generic with a method for `rdm_fit` objects, see [coef.rdm_fit()].
#'
#' @param object A fitted model object.
#' @param ... Arguments passed to methods.
#' @export
fixef <- function(object, ...) UseMethod("fixef")

#' @rdname coef.rdm_fit
#' @export
fixef.rdm_fit <- function(object, type = c("link", "ratio"), ...) {
  coef.rdm_fit(object, type = match.arg(type), full = FALSE)
}

#' Variance-covariance matrix of a fitted RDM
#'
#' The covariance of the fixed effects, computed from the Hessian of the
#' Laplace-approximated marginal likelihood.
#'
#' @param object An `rdm_fit` object.
#' @param full If `TRUE`, the covariance parameters of the random effects are
#'   included, so that the matrix describes all estimated parameters.
#' @param ... Unused.
#' @export
vcov.rdm_fit <- function(object, full = FALSE, ...) {
  V <- object$sd$cov.fixed
  if (is.null(V)) {
    stop("standard errors are not available: the fit was made with ",
         "`sdreport = FALSE` or the Hessian could not be computed.",
         call. = FALSE)
  }
  nms <- c(names(object$beta), names(object$theta))
  dimnames(V) <- list(nms, nms)
  if (!full) {
    nb <- length(object$beta)
    V <- V[seq_len(nb), seq_len(nb), drop = FALSE]
  }
  V
}

#' Wald confidence intervals for a fitted RDM
#'
#' @param object An `rdm_fit` object.
#' @param parm Names or indices of the coefficients; defaults to all of them.
#' @param level Confidence level.
#' @param type `"link"` (default) for intervals on the scale of the
#'   coefficients, `"ratio"` for intervals on the multiplicative scale, obtained
#'   by exponentiating the limits (see [coef.rdm_fit()]). For `b` the interval
#'   refers to the distance between the threshold and `A`.
#' @param ... Unused.
#' @export
confint.rdm_fit <- function(object, parm, level = 0.95,
                            type = c("link", "ratio"), ...) {
  type <- match.arg(type)
  est <- object$beta
  se <- sqrt(diag(vcov(object)))
  if (missing(parm)) {
    parm <- seq_along(est)
  } else if (is.character(parm)) {
    hit <- pmatch(parm, names(est))
    if (anyNA(hit)) {
      stop("unknown coefficient(s): ", paste(parm[is.na(hit)], collapse = ", "),
           ".", call. = FALSE)
    }
    parm <- hit
  }
  z <- stats::qnorm(1 - (1 - level) / 2)
  out <- cbind(est[parm] - z * se[parm], est[parm] + z * se[parm])
  if (type == "ratio") out <- exp(out)
  a <- (1 - level) / 2
  dimnames(out) <- list(names(est)[parm],
                        paste(format(100 * c(a, 1 - a), trim = TRUE), "%"))
  out
}

#' Random effects of a fitted model
#'
#' A generic with a method for `rdm_fit` objects, see [ranef.rdm_fit()].
#'
#' @param object A fitted model object.
#' @param ... Arguments passed to methods.
#' @export
ranef <- function(object, ...) UseMethod("ranef")

#' Random effects of a fitted RDM
#'
#' Conditional modes of the random effects, that is the estimates of the
#' deviations of each level from the population value. They are on the scale of
#' the linear predictor, so they are added to the fixed effects to obtain the
#' level-specific coefficients; `type = "ratio"` exponentiates them, which
#' turns them into the factors relating each level's parameter to the
#' population value, not into parameter values themselves. Use [fitted.rdm_fit()]
#' for the trial-level parameters.
#'
#' @param object An `rdm_fit` object.
#' @param type `"link"` (default) or `"ratio"`.
#' @param ... Unused.
#' @return A data frame with the columns `term`, `group`, `effect` and
#'   `estimate`, or `NULL` for a model without random effects.
#' @export
ranef.rdm_fit <- function(object, type = c("link", "ratio"), ...) {
  type <- match.arg(type)
  blocks <- object$ranef$blocks
  if (!length(blocks)) return(NULL)
  out <- vector("list", length(blocks))
  for (b in seq_along(blocks)) {
    blk <- blocks[[b]]
    u <- object$u[object$nll$index$u[[b]]]
    U <- matrix(u, blk$n_lev, blk$q)
    if (type == "ratio") U <- exp(U)
    out[[b]] <- data.frame(
      term = blk$label,
      group = rep(blk$levels, times = blk$q),
      effect = rep(blk$col_names, each = blk$n_lev),
      estimate = as.vector(U),
      stringsAsFactors = FALSE
    )
  }
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out
}

#' Log-likelihood of a fitted RDM
#'
#' The value is the Laplace approximation of the marginal log-likelihood,
#' evaluated at the estimates and without the optimisation priors. It contains
#' the Laplace correction term, which grows as the random-effect standard
#' deviations approach zero; the approximation is therefore least trustworthy
#' when the estimated variance components are small, and likelihood-based model
#' comparisons should be treated with corresponding caution.
#'
#' @param object An `rdm_fit` object.
#' @param ... Unused.
#' @return An object of class `logLik`.
#' @export
logLik.rdm_fit <- function(object, ...) {
  structure(object$loglik, df = object$df, nobs = object$nobs,
            class = "logLik")
}

#' @export
nobs.rdm_fit <- function(object, ...) object$nobs

#' Trial-level model parameters of a fitted RDM
#'
#' @param object An `rdm_fit` object.
#' @param ... Unused.
#' @return A data frame with the columns `A`, `b`, `t0`, `s` and the drift
#'   rates. When all parameters are shared by the accumulators there is one row
#'   per trial and one column of drift rates per accumulator. When a parameter
#'   varies over accumulators the data frame has one row per trial and
#'   accumulator, as does the `long` component of the design, and the drift rate
#'   is in the column `v`.
#' @export
fitted.rdm_fit <- function(object, ...) {
  mp <- object$pars
  lv <- object$design$response_levels
  ## the drift rates are the only parameter with one column per accumulator in
  ## the wide format; the others need one row per trial and accumulator
  expanded <- any(vapply(c("A", "b", "t0", "s"),
                         \(p) isTRUE(object$design[[p]]$per_accumulator), TRUE))
  out <- if (expanded) {
    out <- object$design$long
    out$A <- as.vector(mp$A)
    out$b <- as.vector(mp$b)
    out$t0 <- as.vector(mp$t0)
    out$s <- as.vector(mp$s)
    out$v <- as.vector(mp$v)
    out
  } else {
    out <- data.frame(A = mp$A[, 1], b = mp$b[, 1], t0 = mp$t0[, 1],
                      s = mp$s[, 1], check.names = FALSE)
    v <- mp$v
    colnames(v) <- paste0("v", lv)
    cbind(out, as.data.frame(v))
  }
  rownames(out) <- NULL
  out
}

#' Predictions from a fitted RDM
#'
#' Predictions for the trials the model was fitted to, since the design matrices
#' are tied to the data.
#'
#' @param object An `rdm_fit` object.
#' @param type `"parameters"` (default) for the trial-level model parameters,
#'   `"probs"` for the probability of each response, or `"quantiles"` for
#'   response-time quantiles conditional on each response.
#' @param probs Quantiles to report for `type = "quantiles"`.
#' @param grid Number of points on the integration grid used for `"probs"` and
#'   `"quantiles"`.
#' @param ... Unused.
#' @return For `"parameters"` a data frame as returned by [fitted.rdm_fit()],
#'   for `"probs"` a matrix with one row per trial and one column per response,
#'   and for `"quantiles"` an array with dimensions trial, response and
#'   quantile.
#' @export
predict.rdm_fit <- function(object, type = c("parameters", "probs", "quantiles"),
                            probs = c(0.1, 0.3, 0.5, 0.7, 0.9), grid = 300L,
                            ...) {
  type <- match.arg(type)
  if (type == "parameters") return(fitted.rdm_fit(object))
  mp <- object$pars
  colnames(mp$v) <- object$design$response_levels
  q <- if (type == "quantiles") probs else NULL
  res <- rdm_race_dist(mp, m = grid, probs = q)
  if (type == "probs") return(res$prob)
  dimnames(res$quantiles)[[3]] <- paste0(100 * probs, "%")
  res$quantiles
}

#' Simulate responses from a fitted RDM
#'
#' Draws new responses and response times from the fitted model for the
#' covariates of the original data.
#'
#' @param object An `rdm_fit` object.
#' @param nsim Number of simulated data sets.
#' @param seed Optional seed for the random number generator.
#' @param ... Unused.
#' @return A data frame with the columns of the original data plus the
#'   simulated `rt` and `response`. With `nsim > 1` the data sets are stacked
#'   and identified by the column `.sim`.
#' @export
simulate.rdm_fit <- function(object, nsim = 1, seed = NULL, ...) {
  if (!is.null(seed)) set.seed(seed)
  levels_response <- object$design$response_levels
  sims <- vector("list", nsim)
  for (s in seq_len(nsim)) {
    draw <- rdm_rand_race(object$pars)
    d <- object$design$data
    d$rt <- draw$rt
    d$response <- factor(levels_response[draw$response], levels = levels_response)
    if (nsim > 1) d$.sim <- s
    sims[[s]] <- d
  }
  out <- do.call(rbind, sims)
  rownames(out) <- NULL
  out
}

## one line describing the spread of the inner Hessian's eigenvalues
rdm_format_eig <- function(eig) {
  if (is.null(eig)) return("not available")
  paste0(format(eig[["min"]], digits = 3), " to ", format(eig[["max"]], digits = 3),
         " (ratio ", format(eig[["min"]] / eig[["max"]], digits = 3), ")")
}

#' @export
print.rdm_fit <- function(x, ...) {
  cat("<rdm_fit> racing diffusion model\n")
  cat("  data: ", x$nobs, " trials, ", x$design$K, " accumulators (",
      paste(x$design$response_levels, collapse = ", "), ")\n", sep = "")
  cat("  fixed effects: ", length(x$beta), ", random-effect parameters: ",
      length(x$theta), ", random effects: ", length(x$u), "\n", sep = "")
  cat("  Laplace log-likelihood: ", format(x$loglik, digits = 7),
      " (df = ", x$df, ")\n", sep = "")
  if (!is.null(x$inner_eig)) {
    cat("  inner Hessian: ", rdm_format_eig(x$inner_eig), "\n", sep = "")
  }
  if (!is.null(x$outer_eig)) {
    cat("  outer Hessian: ", rdm_format_eig(x$outer_eig), "\n", sep = "")
  }
  cat("  optimiser: ", x$opt$method, ", convergence ", x$opt$convergence,
      ", max |gradient| ", format(x$max_grad, digits = 3), "\n", sep = "")
  if (!is.null(x$design$s$anchor)) {
    cat("  diffusion constant s anchored at 1 for ", x$design$s$anchor$level,
        "\n", sep = "")
  }
  cat("  fixed effects (log scale):\n")
  print(round(x$beta, 4))
  invisible(x)
}

#' @export
summary.rdm_fit <- function(object, ...) {
  beta <- object$beta
  se <- tryCatch(sqrt(diag(vcov(object))),
                 error = function(e) rep(NA_real_, length(beta)))
  z <- beta / se
  coefs <- cbind(Estimate = beta, `Std. Error` = unname(se),
                 `z value` = unname(z),
                 `Pr(>|z|)` = unname(2 * stats::pnorm(-abs(z))))
  varcov <- tryCatch(rdm_varcov(object), error = function(e) NULL)
  structure(list(call = object$call, nobs = object$nobs, K = object$design$K,
                 coefficients = coefs, variance = varcov,
                 loglik = object$loglik, df = object$df,
                 inner_eig = object$inner_eig,
                 outer_eig = object$outer_eig,
                 AIC = -2 * object$loglik + 2 * object$df,
                 converged = object$converged, pdHess = object$pdHess,
                 max_grad = object$max_grad, method = object$opt$method),
            class = "summary.rdm_fit")
}

#' @export
print.summary.rdm_fit <- function(x, ...) {
  cat("Call:\n")
  print(x$call)
  cat("\nRacing diffusion model: ", x$nobs, " trials, ", x$K,
      " accumulators\n", sep = "")
  if (!is.null(x$variance)) {
    cat("\nRandom effects:\n")
    tab <- x$variance
    tab$estimate <- round(tab$estimate, 4)
    tab$se <- round(tab$se, 4)
    ## one header per block, so that a long term label is printed once
    for (lab in unique(tab$term)) {
      cat(" ", lab, "\n", sep = "")
      print(tab[tab$term == lab, c("parameter", "estimate", "se")],
            row.names = FALSE)
    }
  }
  cat("\nFixed effects (log scale):\n")
  stats::printCoefmat(x$coefficients, na.print = "NA")
  cat("\nLaplace log-likelihood:", format(x$loglik, digits = 7),
      " on", x$df, "df,  AIC:", format(x$AIC, digits = 7), "\n")
  if (!is.null(x$inner_eig)) {
    cat("Inner Hessian:", rdm_format_eig(x$inner_eig), "\n")
  }
  if (!is.null(x$outer_eig)) {
    cat("Outer Hessian:", rdm_format_eig(x$outer_eig), "\n")
  }
  cat("Optimiser:", x$method, "| converged:", x$converged,
      "| max |gradient|:", format(x$max_grad, digits = 3), "\n")
  if (isFALSE(x$pdHess)) {
    cat("Warning: the joint Hessian is not positive definite; standard errors ",
        "are unreliable.\n", sep = "")
  }
  invisible(x)
}

## Standard deviations and correlations of the random-effect covariance
## matrices, with delta-method standard errors.
rdm_varcov <- function(object) {
  blocks <- object$ranef$blocks
  if (!length(blocks) || is.null(object$sd)) return(NULL)
  nb <- length(object$beta)
  V <- object$sd$cov.fixed
  out <- vector("list", length(blocks))
  for (b in seq_along(blocks)) {
    blk <- blocks[[b]]
    ti <- object$nll$index$theta[[b]]
    th <- object$theta[ti]
    q <- blk$q
    Vb <- V[nb + ti, nb + ti, drop = FALSE]
    C <- rdm_chol(th, q, blk$struct)
    C <- C %*% t(C)
    sds <- sqrt(diag(C))
    rows <- data.frame(
      term = blk$label,
      parameter = paste0("sd(", blk$col_names, ")"),
      estimate = sds,
      se = sds * sqrt(diag(Vb)[rdm_sd_index(q, blk$struct)]),
      stringsAsFactors = FALSE
    )
    ## only the unstructured and compound-symmetric blocks have correlations
    if (q > 1L && blk$struct %in% c("us", "cs")) {
      extra <- list()
      for (i in 2:q) {
        for (j in seq_len(i - 1L)) {
          cor_fun <- function(x) {
            L <- rdm_chol(x, q, blk$struct)
            Cx <- L %*% t(L)
            Cx[i, j] / sqrt(Cx[i, i] * Cx[j, j])
          }
          extra[[length(extra) + 1L]] <- data.frame(
            term = blk$label,
            parameter = paste0("cor(", blk$col_names[i], ", ",
                               blk$col_names[j], ")"),
            estimate = cor_fun(th),
            se = delta_se(cor_fun, th, Vb),
            stringsAsFactors = FALSE
          )
        }
      }
      rows <- rbind(rows, do.call(rbind, extra))
    }
    out[[b]] <- rows
  }
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out
}

## Standard error of a smooth function of `par` from its covariance `V`
delta_se <- function(fun, par, V, eps = 1e-6) {
  val <- fun(par)
  J <- matrix(0, length(val), length(par))
  for (k in seq_along(par)) {
    step <- eps * max(1, abs(par[k]))
    up <- par
    up[k] <- up[k] + step
    dn <- par
    dn[k] <- dn[k] - step
    J[, k] <- (fun(up) - fun(dn)) / (2 * step)
  }
  sqrt(pmax(diag(J %*% V %*% t(J)), 0))
}
