# Design matrices for the RDM.
#
# Every model parameter (A, b, t0, v, s) gets its own model matrix, built from
# a user-supplied formula and the data. A formula can be given once, in which
# case the parameter is shared by all accumulators, or as a list with one
# formula per accumulator, in which case every accumulator gets its own
# coefficients.

rdm_pars <- c("A", "b", "t0", "v", "s")

#' Build the design matrices of an RDM model
#'
#' Turns one formula per model parameter into a model matrix. The design is
#' stored in accumulator-major order for parameters with one formula per
#' accumulator, that is rows `(j - 1) * n + 1:n` belong to accumulator `j`.
#'
#' @param data A `data.frame` with one row per trial.
#' @param rt Numeric vector of response times.
#' @param response Vector of observed responses; the sorted unique values define
#'   the accumulators and are coded `1:K`.
#' @param formula Named list of formulas, one per free parameter. Valid names
#'   are `"A"`, `"b"`, `"t0"`, `"v"` and `"s"`. Every formula may either be a
#'   single formula, in which case the parameter is the same for all
#'   accumulators, or a list of `K` formulas, one per accumulator. A single
#'   formula for the drift rates is evaluated on the expanded data, so the
#'   column `response` is available and refers to the accumulator.
#' @param fixed Named list of parameters that are held constant, on their
#'   natural scale, e.g. `list(s = 1)`. A fixed `b` is the threshold itself, not
#'   the threshold distance. The value can have length 1 (shared) or `K` (one
#'   value per accumulator).
#' @param matchfun Optional function of the expanded data that returns a logical
#'   vector, `TRUE` where the accumulator of a row is the one the stimulus calls
#'   for, e.g. `function(d) d$S == d$lR`. The expanded data always contains the
#'   accumulator as the factor `lR`; with `matchfun` the logical is added as the
#'   factor `lM`, so that formulas can use the latent match.
#' @return An object of class `rdm_design`.
#' @keywords internal
#' @noRd
rdm_design <- function(data, rt, response, formula, fixed = list(s = 1),
                       matchfun = NULL) {
  if (!is.data.frame(data)) stop("`data` must be a data.frame.", call. = FALSE)
  n <- nrow(data)
  if (length(rt) != n) stop("`rt` must have one value per row of `data`.", call. = FALSE)
  if (length(response) != n) {
    stop("`response` must have one value per row of `data`.", call. = FALSE)
  }
  if (anyNA(rt) || anyNA(response)) stop("`rt` and `response` cannot be missing.", call. = FALSE)
  if (any(rt <= 0)) stop("all response times must be positive.", call. = FALSE)

  ## sort numerically for numeric responses, so that 10 does not come before 2
  levels_response <- sort(unique(if (is.factor(response)) as.character(response) else response))
  levels_response <- as.character(levels_response)
  K <- length(levels_response)
  if (K < 2) stop("the RDM needs at least two response alternatives.", call. = FALSE)
  resp_idx <- match(as.character(response), levels_response)

  if (length(formula) && (is.null(names(formula)) || any(names(formula) == ""))) {
    stop("`formula` must be a named list, e.g. list(A = ~1, v = ~1).", call. = FALSE)
  }
  unknown <- setdiff(names(formula), rdm_pars)
  if (length(unknown)) {
    stop("unknown parameter(s) in `formula`: ", paste(unknown, collapse = ", "),
         ". Valid names are ", paste(rdm_pars, collapse = ", "), ".", call. = FALSE)
  }
  unknown <- setdiff(names(fixed), rdm_pars)
  if (length(unknown)) {
    stop("unknown parameter(s) in `fixed`: ", paste(unknown, collapse = ", "), ".",
         call. = FALSE)
  }
  both <- intersect(names(formula), names(fixed))
  if (length(both)) {
    stop("parameter(s) ", paste(both, collapse = ", "),
         " appear in both `formula` and `fixed`.", call. = FALSE)
  }
  if (is.null(formula$v) && !"v" %in% names(fixed)) {
    stop("the drift rates `v` need a formula.", call. = FALSE)
  }

  ## expand the data to one row per trial and accumulator
  long <- data[rep(seq_len(n), K), , drop = FALSE]
  long[["response"]] <- factor(rep(levels_response, each = n), levels = levels_response)
  ## the accumulator as a factor, as it is called in the EMC2 modelling framework
  long[["lR"]] <- long[["response"]]
  if (!is.null(matchfun)) {
    if (!is.function(matchfun)) stop("`matchfun` must be a function.", call. = FALSE)
    m <- matchfun(long)
    if (!is.logical(m) || length(m) != nrow(long)) {
      stop("`matchfun` must return a logical vector with one value per trial and ",
           "accumulator.", call. = FALSE)
    }
    long[["lM"]] <- factor(m, levels = c(FALSE, TRUE))
  }

  design <- list(n = n, K = K, rt = rt, response = resp_idx,
                 response_levels = levels_response, data = data, long = long,
                 formula = formula, fixed = fixed)

  for (par in rdm_pars) {
    design[[par]] <- if (par %in% names(fixed)) {
      fixed_spec(fixed[[par]], n, par, K)
    } else if (par %in% names(formula)) {
      par_spec(formula[[par]], data, long, n, K, par, levels_response)
    } else {
      NULL
    }
  }

  check_identified(design)
  structure(design, class = c("rdm_design", "list"))
}

## the design of one free parameter
##
## `X` has one row per trial for a shared parameter and one row per trial and
## accumulator for a parameter that is expanded to the accumulators
par_spec <- function(form, data, long, n, K, par, levels) {
  if (inherits(form, "formula")) {
    ## a single drift-rate formula is evaluated on the expanded data so that the
    ## accumulator is available as the column `response`
    obs <- if (par == "v") long else data
    X <- stats::model.matrix(form, obs)
    if (ncol(X) == 0) {
      stop("the formula for `", par, "` has no coefficients (did you use `~ 0`?).",
           call. = FALSE)
    }
    if (nrow(X) != nrow(obs)) {
      stop("internal error: wrong number of rows for `", par, "`.", call. = FALSE)
    }
    rownames(X) <- NULL
    return(list(free = TRUE, X = X, coef_names = paste0(par, ":", colnames(X)),
                contrasts = attr(X, "contrasts"), formula = form,
                expanded = nrow(X) != n, per_accumulator = FALSE))
  }

  if (!is.list(form) || length(form) != K ||
        !all(vapply(form, inherits, TRUE, "formula"))) {
    stop("the formula for `", par, "` must be a single formula or a list of ",
         K, " formulas, one per accumulator.", call. = FALSE)
  }
  blocks <- lapply(seq_len(K), \(j) stats::model.matrix(form[[j]], data))
  p_j <- vapply(blocks, ncol, 1L)
  if (any(p_j == 0)) {
    stop("every formula for `", par, "` needs at least one coefficient.",
         call. = FALSE)
  }
  X <- matrix(0, nrow = n * K, ncol = sum(p_j))
  colnames(X) <- unlist(lapply(seq_len(K), \(j) {
    paste0(par, "[", levels[j], "]:", colnames(blocks[[j]]))
  }))
  start <- cumsum(c(0, p_j))
  for (j in seq_len(K)) {
    rows <- (j - 1) * n + seq_len(n)
    X[rows, start[j] + seq_len(p_j[j])] <- blocks[[j]]
  }
  list(free = TRUE, X = X, coef_names = colnames(X), contrasts = NULL,
       formula = form, expanded = TRUE, per_accumulator = TRUE)
}

## a parameter held constant: either one value for all accumulators or one value
## per accumulator, recycled to the required length
fixed_spec <- function(value, n, par, K) {
  if (!is.numeric(value) || anyNA(value) || any(value <= 0)) {
    stop("the fixed value for `", par, "` must be a positive number.", call. = FALSE)
  }
  if (length(value) == 1) {
    value <- rep(value, n)
    expanded <- FALSE
  } else if (length(value) == K) {
    value <- rep(value, each = n)
    expanded <- TRUE
  } else {
    stop("the fixed value for `", par, "` must have length 1 or ", K, ".",
         call. = FALSE)
  }
  list(free = FALSE, value = value, coef_names = character(), formula = NULL,
       contrasts = NULL, expanded = expanded, per_accumulator = expanded)
}

## the diffusion constant s is only identified if the scale of the model is
## pinned down by at least one of A, b and v
check_identified <- function(design) {
  if (!isTRUE(design$s$free)) return(invisible(TRUE))
  pinned <- intersect(names(design$fixed), c("A", "b", "v"))
  if (!length(pinned)) {
    stop("a free `s` (the diffusion constant) is not identified: `s` scales ",
         "`A`, `b` and `v`, so at least one of them must be fixed, e.g. ",
         "fixed = list(A = 0.5).", call. = FALSE)
  }
  invisible(TRUE)
}

#' @export
print.rdm_design <- function(x, ...) {
  cat("<rdm_design> ", x$n, " trials, ", x$K, " accumulators (",
      paste(x$response_levels, collapse = ", "), ")\n", sep = "")
  for (par in rdm_pars) {
    spec <- x[[par]]
    if (is.null(spec)) {
      cat("  ", par, ": not in the model\n", sep = "")
    } else if (!spec$free) {
      if (spec$expanded) {
        cat("  ", par, ": fixed per accumulator at ",
            paste(unique(spec$value), collapse = ", "), "\n", sep = "")
      } else {
        cat("  ", par, ": fixed at ", unique(spec$value)[1], "\n", sep = "")
      }
    } else {
      cat("  ", par, ": ", ncol(spec$X), " coefficient(s)",
          if (spec$per_accumulator) ", one set per accumulator" else "",
          ": ", paste(spec$coef_names, collapse = ", "), "\n", sep = "")
    }
  }
  invisible(x)
}
