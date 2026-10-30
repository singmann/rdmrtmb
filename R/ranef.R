# Random effects for the RDM.
#
# Random effects are specified per parameter with lme4-style terms such as
# `(1 | id)` or `(1 + cond | id)`. Every term becomes its own block with a
# Cholesky-parameterised covariance matrix, so `(1 | id) + (0 + cond | id)`
# gives independent effects while `(1 + cond | id)` gives correlated ones.
# Several grouping factors may be combined in one formula, which allows crossed
# random effects.

#' Parse the random-effects terms of an RDM model
#'
#' @param random A one-sided formula (applied to every parameter) or a named
#'   list of one-sided formulas, one per parameter. Each formula may contain
#'   several terms, e.g. `~(1 | id) + (1 | item)`.
#' @param data The trial-level `data.frame` of the model.
#' @param design An `rdm_design` object.
#' @return An object of class `rdm_ranef`.
#' @keywords internal
#' @noRd
rdm_ranef <- function(random, data, design) {
  if (is.null(random)) return(structure(list(blocks = list(), npar = 0L),
                                        class = c("rdm_ranef", "list")))
  if (inherits(random, "formula")) {
    ## a single formula applies to every free parameter of the model
    used <- rdm_pars[vapply(rdm_pars, \(nm) !is.null(design[[nm]]) &&
                              design[[nm]]$free, TRUE)]
    random <- stats::setNames(rep(list(random), length(used)), used)
  }
  if (!is.list(random) || is.null(names(random))) {
    stop("`random` must be a formula or a named list of formulas.", call. = FALSE)
  }
  unknown <- setdiff(names(random), rdm_pars)
  if (length(unknown)) {
    stop("unknown parameter(s) in `random`: ", paste(unknown, collapse = ", "), ".",
         call. = FALSE)
  }
  for (par in names(random)) {
    if (is.null(design[[par]])) {
      stop("`random` specifies effects for `", par, "`, which is not part of the ",
           "model.", call. = FALSE)
    }
    if (!design[[par]]$free) {
      stop("`", par, "` is held constant, so it cannot have random effects.",
           call. = FALSE)
    }
  }

  blocks <- list()
  for (par in names(random)) {
    form <- random[[par]]
    if (!inherits(form, "formula") || length(form) != 2L) {
      stop("the random-effects formula for `", par, "` must be one-sided, e.g. ",
           "~(1 | id).", call. = FALSE)
    }
    ## random effects on a parameter that varies over trials and accumulators
    ## are evaluated on the expanded data
    re_data <- if (isTRUE(design[[par]]$expanded)) design$long else design$data
    terms <- reformulas::findbars(form)
    if (!length(terms)) {
      stop("the random-effects formula for `", par, "` has no `(term | group)` ",
           "specification.", call. = FALSE)
    }
    for (term in terms) {
      blocks[[length(blocks) + 1L]] <- ranef_block(term, re_data, par)
    }
  }

  npar <- sum(vapply(blocks, \(b) b$n_theta, 1L))
  structure(list(blocks = blocks, npar = npar, random = random),
            class = c("rdm_ranef", "list"))
}

## one `(expr | group)` term
ranef_block <- function(term, data, par) {
  if (!identical(as.character(term[[1]]), "|")) {
    stop("only single-bar terms such as (1 | id) are supported.", call. = FALSE)
  }
  group_expr <- term[[3]]
  re_expr <- term[[2]]

  group <- tryCatch(eval(group_expr, data, parent.frame()),
                    error = function(e) {
                      stop("the grouping variable `", deparse1(group_expr),
                           "` was not found in `data`.", call. = FALSE)
                    })
  if (length(group) != nrow(data)) {
    stop("the grouping variable `", deparse1(group_expr),
         "` must have one value per trial.", call. = FALSE)
  }
  if (anyNA(group)) stop("the grouping variable `", deparse1(group_expr),
                         "` contains missing values.", call. = FALSE)
  ## keep the original type for sorting so that numeric ids are ordered
  ## numerically rather than as strings
  if (is.factor(group)) group <- as.character(group)
  levels <- sort(unique(group))
  group <- as.character(group)
  levels <- as.character(levels)
  group_idx <- match(group, levels)

  Z <- stats::model.matrix(stats::as.formula(paste("~", deparse1(re_expr))), data)
  if (ncol(Z) == 0) {
    stop("the random-effects term `", deparse1(term), "` has no columns.",
         call. = FALSE)
  }
  ## a term with only an intercept column is stored as a single column of ones
  q <- ncol(Z)
  list(par = par,
       label = paste0(par, ":(", deparse1(term), ")"),
       Z = matrix(as.numeric(Z), nrow = nrow(Z), ncol = ncol(Z)),
       col_names = colnames(Z),
       group = group_idx,
       group_expr = deparse1(group_expr),
       levels = levels,
       n_lev = length(levels),
       n_row = nrow(Z),
       q = q,
       n_theta = as.integer(q * (q + 1L) / 2L),
       correlated = q > 1L)
}

## Position of every element of a q x q matrix in the vector of lower-triangular
## parameters: the first q entries are the log standard deviations, the rest are
## the off-diagonal elements. Elements above the diagonal get a dummy index.
lower_map <- function(q) {
  map <- rep(1L, q * q)
  lt <- which(lower.tri(diag(q)), arr.ind = TRUE)
  ord <- order(lt[, "col"], lt[, "row"])
  map[lt[ord, "row"] + (lt[ord, "col"] - 1L) * q] <- q + seq_len(nrow(lt))
  map
}

## Cholesky factor of the random-effect covariance matrix from the
## unconstrained parameters theta. Assignment into a numeric matrix would drop
## the automatic-differentiation class, so the factor is assembled from a
## constant index map instead.
rdm_chol <- function(theta, q) {
  if (q == 1L) return(matrix(exp(theta[1]), 1, 1))
  mask <- lower.tri(diag(q))
  unit <- matrix(theta[lower_map(q)], q, q) * mask + diag(q)
  ## scaling the rows by the standard deviations gives L with diag(L) = sd
  unit * exp(theta[seq_len(q)])
}

## Solve L Z = U for a lower triangular L, vectorised over the rows of U.
## RTMB does not support solve() for AD matrices, so forward substitution is
## spelled out explicitly.
tri_solve <- function(L, U) {
  q <- ncol(L)
  if (q == 1L) return(U / L[1, 1])
  z <- vector("list", q)
  z[[1]] <- U[, 1] / L[1, 1]
  for (k in seq(2, q)) {
    acc <- U[, k]
    for (i in seq_len(k - 1L)) acc <- acc - L[k, i] * z[[i]]
    z[[k]] <- acc / L[k, k]
  }
  do.call(cbind, z)
}

## contribution of one block to the linear predictor, one value per row of the
## block's design matrix
rdm_ranef_contrib <- function(block, u) {
  U <- matrix(u, block$n_lev, block$q)
  rowSums(block$Z * U[block$group, , drop = FALSE])
}

## -log density of the random effects of one block, including the
## normalising constant
rdm_ranef_penalty <- function(block, u, theta) {
  q <- block$q
  L <- rdm_chol(theta, q)
  Z <- tri_solve(L, matrix(u, block$n_lev, q))
  0.5 * sum(Z^2) + block$n_lev * sum(theta[seq_len(q)]) +
    block$n_lev * q / 2 * log(2 * pi)
}

#' @export
print.rdm_ranef <- function(x, ...) {
  if (!length(x$blocks)) {
    cat("<rdm_ranef> no random effects\n")
    return(invisible(x))
  }
  cat("<rdm_ranef> ", length(x$blocks), " term(s), ", x$npar,
      " covariance parameter(s)\n", sep = "")
  for (b in x$blocks) {
    cat("  ", b$label, ": ", b$n_lev, " levels, ", b$q, " effect(s)",
        if (b$correlated) ", correlated" else "", "\n", sep = "")
  }
  invisible(x)
}
