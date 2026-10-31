# A quick first look at a hierarchical model: fit every subject on its own and
# summarise the estimates. The per-subject fits take a fraction of a second
# each, so they give both a picture of the group distribution and starting
# values for the slower hierarchical fit.

#' Starting values from independent per-subject fits
#'
#' Fits the same fixed-effect model to each subject separately, with the
#' random-effect terms removed, and summarises the resulting group distribution.
#' The result can be passed as the `start` argument of [fit_rdm()], which saves
#' a lot of time on larger models: an independent fit takes a fraction of a
#' second, whereas every evaluation of the hierarchical likelihood can take
#' seconds.
#'
#' The fixed effects of the hierarchical model start at the mean of the
#' per-subject estimates and the random-effect standard deviations start at
#' their between-subject standard deviation; correlations start at zero. Both
#' are on the log scale, which is where the random effects are additive. The
#' random effects start at each subject's deviation from the group mean, so the
#' first evaluation of the likelihood is already close to a per-subject
#' solution. Coefficients that a per-subject fit cannot estimate, such as a
#' covariate that only varies between subjects, or that not every subject could
#' estimate, keep their default start.
#'
#' @inheritParams fit_rdm
#' @param id The grouping variable, i.e. the variable a random effect is
#'   grouped by: a vector, or the name of a column of `data`.
#' @param prior The prior for the per-subject fits. The weak prior on the fixed
#'   effects is kept, and the one on the standard deviations is not needed
#'   because these fits have no random effects. Without a prior the per-subject
#'   likelihood can have degenerate maxima when a subject has few trials, where
#'   a large start-point range, threshold and drift rate trade off against each
#'   other.
#' @param ... Further arguments passed on to [fit_rdm()] for the per-subject
#'   fits.
#' @return An object of class `rdm_start`, a list with the elements `beta`,
#'   `theta` and `u`, or as many of them as the model has. The per-subject
#'   coefficient matrix is stored in the `coef` attribute, and the per-subject
#'   fits in the `fits` attribute.
#' @seealso [fit_rdm()]
#' @export
rdm_subject_start <- function(data, rt, response, formula, id,
                              fixed = list(s = 1), contrasts = list(),
                              anchor = list(), prior = list(sd = NULL, beta = 2),
                              control = list(), ...) {
  if (!is.data.frame(data)) stop("`data` must be a data.frame.", call. = FALSE)
  rt <- rdm_column(rt, data, "rt")
  response <- rdm_column(response, data, "response")
  id <- rdm_column(id, data, "id")
  if (length(id) != nrow(data)) {
    stop("`id` must have one value per row of `data`.", call. = FALSE)
  }
  formula <- as_rdm_formula(formula)
  ## the same fixed effects, without the random-effect terms
  sub_formula <- lapply(unclass(formula), reformulas::nobars)

  rows <- split(seq_len(nrow(data)), as.character(id))
  fit_one <- function(i) {
    tryCatch(
      suppressWarnings(
        fit_rdm(data[i, , drop = FALSE], rt = rt[i], response = response[i],
                formula = sub_formula, fixed = fixed, contrasts = contrasts,
                anchor = anchor, prior = prior, control = control,
                sdreport = FALSE, ...)
      ),
      error = function(e) NULL
    )
  }
  fits <- lapply(rows, fit_one)
  ok <- !vapply(fits, is.null, TRUE)
  if (sum(ok) < 2L) {
    stop("fewer than two subjects could be fitted; try fitting one subject ",
         "with fit_rdm() to see the error.", call. = FALSE)
  }

  ## the hierarchical model decides which coefficients and which blocks exist
  design <- rdm_design(data, rt, response, formula, fixed, contrasts, anchor)
  ranef <- rdm_ranef(design)
  nll <- rdm_nll(design, ranef)
  beta <- nll$par$beta

  coefs <- matrix(NA_real_, length(fits), length(beta),
                  dimnames = list(names(rows), names(beta)))
  for (s in which(ok)) {
    b <- fits[[s]]$beta
    coefs[s, names(b)] <- b
  }
  mean_coef <- colMeans(coefs, na.rm = TRUE)
  known <- is.finite(mean_coef)
  beta[known] <- mean_coef[known]

  theta <- nll$par$theta
  for (j in seq_along(ranef$blocks)) {
    block <- ranef$blocks[[j]]
    pos <- nll$index$theta[[j]]
    sds <- vapply(block$col_names, \(nm) between_sd(coefs[, nm]), 0)
    shared <- block$struct %in% c("homdiag", "cs")
    values <- if (shared) sqrt(mean(sds^2, na.rm = TRUE)) else sds
    sd_pos <- rdm_sd_index(block$q, block$struct)
    for (k in seq_along(sd_pos)) {
      ## a variance of zero, or no estimate at all, keeps the default
      if (!is.finite(values[k]) || values[k] <= 0) next
      theta[pos[sd_pos[k]]] <- log(values[k])
    }
  }

  ## the random effects start at each subject's own deviation from the group
  ## mean, which is the scale `u` is on. With zero deviations every subject gets
  ## the group parameter values, and a shared `t0` is then above the shortest
  ## response time of some subject, which the barrier punishes so hard that the
  ## inner optimisation fails.
  u <- nll$par$u
  if (!is.null(u)) {
    for (j in seq_along(ranef$blocks)) {
      block <- ranef$blocks[[j]]
      dev <- matrix(0, block$n_lev, block$q)
      row <- match(block$levels, rownames(coefs))
      for (k in seq_len(block$q)) {
        nm <- block$col_names[k]
        if (!nm %in% colnames(coefs) || is.na(beta[nm])) next
        x <- coefs[row, nm] - beta[nm]
        dev[, k] <- ifelse(is.na(x), 0, x)
      }
      u[nll$index$u[[j]]] <- as.numeric(dev)
    }
  }

  start <- list(beta = beta, theta = theta, u = u)
  structure(start[!vapply(start, is.null, TRUE)],
            class = "rdm_start",
            coef = coefs[ok, , drop = FALSE],
            fits = fits[ok],
            n_failed = sum(!ok))
}

## the standard deviation of a coefficient across subjects, or NA for fewer than
## two estimates
between_sd <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) < 2L) NA_real_ else stats::sd(x)
}

#' @export
print.rdm_start <- function(x, ...) {
  coefs <- attr(x, "coef")
  cat("<rdm_start> from ", nrow(coefs), " subject fits", sep = "")
  if (attr(x, "n_failed") > 0) {
    cat(", ", attr(x, "n_failed"), " failed", sep = "")
  }
  cat("\n\nGroup distribution (back-transformed):\n")
  tab <- data.frame(
    estimate = exp(colMeans(coefs, na.rm = TRUE)),
    sd = exp(apply(coefs, 2, between_sd)),
    row.names = NULL, check.names = FALSE
  )
  tab <- cbind(coefficient = colnames(coefs), tab)
  print(tab, row.names = FALSE)
  cat("\nThe estimate is exp() of the mean of the log-scale coefficients and",
      "the standard deviation is exp() of their standard deviation, i.e. the",
      "multiplicative between-subject variation.\n")
  invisible(x)
}
