# Random effects for the RDM.
#
# Random effects are written inside the parameter formulas with lme4-style
# terms such as `(1 | id)` or `(1 + cond | id)`. Every term becomes its own
# block with a Cholesky-parameterised covariance matrix, so
# `(1 | id) + (0 + cond | id)` gives independent effects while `(1 + cond | id)`
# gives correlated ones. Several grouping factors may be combined in one
# formula, which allows crossed random effects.
#
# The covariance structure is chosen by wrapping the term in a function, as in
# glmmTMB:
#
#   * `us` (the default, written as a plain term) estimates all variances and
#     correlations;
#   * `diag` gives one variance per column of the term and no correlations,
#     which is the only way to get independent effects for a factor, because
#     `(0 + f | id)` still estimates a covariance between the levels of `f`;
#   * `homdiag` gives a single variance that is shared by all columns;
#   * `cs` (compound symmetry) gives a single variance and a single correlation
#     between any two columns.
#
# The left-hand side of the bar is built like the fixed part of the formula, so
# `acc()` and the accumulator variables `response` and `.acc` can be used there
# too. `(acc(Condition) | id)` therefore gives the same set of effects as
# `acc(Condition)` in the formula, one per accumulator and condition.

## the covariance structures a term can have, by the name of the wrapper
rdm_structs <- c(us = "us", diag = "diag", homdiag = "homdiag", cs = "cs")

## number of covariance parameters of a block with `q` columns
rdm_n_theta <- function(q, struct = "us") {
  if (identical(struct, "homdiag")) return(1L)
  if (identical(struct, "cs")) return(if (q > 1L) 2L else 1L)
  if (identical(struct, "diag")) return(as.integer(q))
  as.integer(q * (q + 1L) / 2L)
}

## positions in the parameter vector of a block that are log standard
## deviations; the remaining ones are correlations
rdm_sd_index <- function(q, struct = "us") {
  if (identical(struct, "homdiag") || identical(struct, "cs")) return(1L)
  seq_len(q)
}

## The random-effect terms of one parameter formula, each with its covariance
## structure.
rdm_ranef_terms <- function(form) {
  env <- environment(form)
  if (is.null(env)) env <- parent.frame()
  labels <- attr(stats::terms(form), "term.labels")
  out <- list()
  for (label in labels) {
    if (!grepl("|", label, fixed = TRUE)) next
    expr <- str2lang(label)
    struct <- "us"
    if (!as.character(expr[[1]]) %in% c("|", "||")) {
      name <- as.character(expr[[1]])
      if (!name %in% names(rdm_structs)[-1] || !is.call(expr[[2]])) {
        stop("unsupported covariance structure `", name, "()` in the term `",
             label, "`. Supported are `(term | group)` and the wrappers ",
             paste0("`", names(rdm_structs)[-1], "(term | group)`",
                    collapse = ", "), ".", call. = FALSE)
      }
      struct <- rdm_structs[[name]]
      expr <- expr[[2]]
    }
    ## findbars() expands double bars into independent single-bar terms
    for (bar in reformulas::findbars(eval(call("~", expr), env))) {
      out[[length(out) + 1L]] <- list(term = bar, struct = struct)
    }
  }
  out
}

#' Parse the random-effects terms of an RDM model
#'
#' The random effects are read from the `(term | group)` parts of the parameter
#' formulas, which `rdm_design()` has already split off with
#' `reformulas::findbars()`. A term of a parameter that is held constant is
#' rejected. Printing the result shows every block with its grouping factor, the
#' number of levels, the number of effects in the block and its covariance
#' structure, which is a quick way to check how many covariance parameters a
#' model has.
#'
#' @param design An `rdm_design` object.
#' @return An object of class `rdm_ranef`, a list with one block per term, each
#'   with its design matrix `Z`, the grouping index, the number of effects `q`,
#'   the covariance structure `struct` and the number of covariance parameters
#'   `n_theta`.
#' @seealso [fit_rdm()], [rdm_design()]
#' @export
rdm_ranef <- function(design) {
  terms <- design$ranef_terms
  if (!length(terms)) {
    return(structure(list(blocks = list(), npar = 0L),
                     class = c("rdm_ranef", "list")))
  }

  blocks <- list()
  for (par in names(terms)) {
    if (is.null(design[[par]]) || !design[[par]]$free) {
      stop("`", par, "` is held constant, so it cannot have random effects.",
           call. = FALSE)
    }
    ## the design matrices live on the expanded data, so the random effects are
    ## evaluated there as well
    for (term in terms[[par]]) {
      blocks[[length(blocks) + 1L]] <- ranef_block(term$term, design, par,
                                                   term$struct)
    }
  }

  npar <- sum(vapply(blocks, \(b) b$n_theta, 1L))
  structure(list(blocks = blocks, npar = npar), class = c("rdm_ranef", "list"))
}

## one `(expr | group)` term. The effects are built with the same machinery as
## the fixed effects, so the accumulator terms of the fixed part can be repeated
## in the random part: `acc(Condition)` in a bar gives one effect per
## accumulator and condition, exactly like `acc(Condition)` in the formula.
ranef_block <- function(term, design, par, struct = "us") {
  if (!identical(as.character(term[[1]]), "|")) {
    stop("only single-bar terms such as (1 | id) are supported.", call. = FALSE)
  }
  group_expr <- term[[3]]
  re_expr <- term[[2]]
  ## the design matrices live on the expanded data, so the random effects are
  ## evaluated there as well, where `.acc` and `response` are available
  data <- design$long
  env <- environment(design$formula[[par]])
  if (is.null(env)) env <- parent.frame()

  group <- tryCatch(eval(group_expr, data, env),
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

  ## the columns of the block, one row per row of the expanded data; only the
  ## values are kept, the names are stored separately
  spec <- par_spec(stats::as.formula(paste("~", deparse1(re_expr)), env = env),
                   design$data, data, design$n, design$K, par,
                   design$response_levels, design$contrasts,
                   what = paste0("the random-effects term (", deparse1(term), ")"))
  Z <- matrix(as.numeric(spec$X), nrow = nrow(spec$X), ncol = ncol(spec$X))
  q <- ncol(Z)
  wrapped <- !identical(struct, "us")
  list(par = par,
       label = if (wrapped) paste0(par, ":(", struct, "(", deparse1(term), "))")
               else paste0(par, ":(", deparse1(term), ")"),
       struct = struct,
       Z = Z,
       col_names = spec$coef_names,
       group = group_idx,
       group_expr = deparse1(group_expr),
       levels = levels,
       n_lev = length(levels),
       n_row = nrow(Z),
       q = q,
       n_theta = rdm_n_theta(q, struct),
       correlated = q > 1L && identical(struct, "us"))
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
rdm_chol <- function(theta, q, struct = "us") {
  ## one variance shared by all columns
  if (identical(struct, "homdiag")) return(diag(q) * exp(theta[1]))
  ## compound symmetry: the Cholesky factor of (1 - rho) I + rho J is lower
  ## triangular with equal elements below the diagonal of every column, so it
  ## can be assembled from the constant masks below
  if (identical(struct, "cs") && q > 1L) {
    sigma <- exp(theta[1])
    a <- 1
    b <- cs_rho(theta[2], q)
    L <- matrix(0, q, q)
    for (k in seq_len(q)) {
      dk <- sqrt(a)
      L <- L + basis_diag(q, k) * (sigma * dk)
      if (k < q) L <- L + basis_lower(q, k) * (sigma * b / dk)
      a0 <- a
      a <- a0 - b^2 / a0
      b <- b - b^2 / a0
    }
    return(L)
  }
  if (q == 1L) return(matrix(exp(theta[1]), 1, 1))
  ## a diagonal covariance matrix: one standard deviation per column and no
  ## off-diagonal elements, so there is nothing to fill in
  if (identical(struct, "diag")) return(diag(q) * exp(theta[seq_len(q)]))
  mask <- lower.tri(diag(q))
  unit <- matrix(theta[lower_map(q)], q, q) * mask + diag(q)
  ## scaling the rows by the standard deviations gives L with diag(L) = sd
  unit * exp(theta[seq_len(q)])
}

## Map an unconstrained parameter to the correlation of a compound-symmetric
## covariance matrix with `q` columns. The matrix is positive definite exactly
## for rho in (-1/(q-1), 1), so the map covers that interval and theta = 0
## gives a correlation of zero. The boundary is at theta = -Inf, which lets the
## prior on theta keep the correlation away from the singular matrix.
cs_rho <- function(theta, q) {
  lower <- -1 / (q - 1)
  span <- 1 - lower
  ## the logistic link is written out because RTMB differentiates exp() but not
  ## stats::plogis()
  lower + span / (1 + exp(-(theta + stats::qlogis(-lower / span))))
}

## the inverse of cs_rho(), used for starting values
cs_theta <- function(rho, q) {
  lower <- -1 / (q - 1)
  span <- 1 - lower
  stats::qlogis((rho - lower) / span) - stats::qlogis(-lower / span)
}

## constant masks used to assemble the Cholesky factor of a compound-symmetric
## covariance matrix: the k-th diagonal position and the k-th column below it
basis_diag <- function(q, k) {
  m <- matrix(0, q, q)
  m[k, k] <- 1
  m
}

basis_lower <- function(q, k) {
  m <- matrix(0, q, q)
  if (k < q) m[(k + 1L):q, k] <- 1
  m
}

## the log determinant of the covariance matrix of a block, which is twice the
## sum of the logs of the diagonal of its Cholesky factor
rdm_chol_logdet <- function(theta, q, struct = "us") {
  if (identical(struct, "homdiag")) return(q * theta[1])
  if (identical(struct, "cs")) {
    if (q == 1L) return(theta[1])
    rho <- cs_rho(theta[2], q)
    return(q * theta[1] +
             0.5 * ((q - 1) * log1p(-rho) + log1p((q - 1) * rho)))
  }
  sum(theta[seq_len(q)])
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
  L <- rdm_chol(theta, q, block$struct)
  Z <- tri_solve(L, matrix(u, block$n_lev, q))
  0.5 * sum(Z^2) + block$n_lev * rdm_chol_logdet(theta, q, block$struct) +
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
        if (!identical(b$struct, "us")) paste0(", ", b$struct)
        else if (b$correlated) ", correlated" else "", "\n", sep = "")
  }
  invisible(x)
}
