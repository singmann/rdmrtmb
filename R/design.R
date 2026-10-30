# Design matrices for the RDM.
#
# Every model parameter (A, b, t0, v, s) gets its own model matrix, built from
# a user-supplied formula and the data. The drift rates additionally vary across
# accumulators, so their design has one row per trial and accumulator.

rdm_pars <- c("A", "b", "t0", "v", "s")

#' Build the design matrices of an RDM model
#'
#' Turns one formula per model parameter into a model matrix. The design is
#' stored in accumulator-major order for the drift rates, that is rows
#' `(j - 1) * n + 1:n` belong to accumulator `j`.
#'
#' @param data A `data.frame` with one row per trial.
#' @param rt Numeric vector of response times.
#' @param response Vector of observed responses; the sorted unique values define
#'   the accumulators and are coded `1:K`.
#' @param formula Named list of formulas, one per free parameter. Valid names
#'   are `"A"`, `"b"`, `"t0"`, `"v"` and `"s"`. For `v` either a single formula
#'   is given, in which case the column `response` refers to the accumulator,
#'   or a list of `K` formulas with one formula per accumulator.
#' @param fixed Named list of parameters that are held constant, on their
#'   natural scale, e.g. `list(s = 1)`. A fixed `b` is the threshold itself, not
#'   the threshold distance. A constant for `v` can be a vector of length `K`.
#' @return An object of class `rdm_design`.
#' @keywords internal
#' @noRd
rdm_design <- function(data, rt, response, formula, fixed = list(s = 1)) {
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

  ## drift rates: expand the data to one row per trial and accumulator
  long <- data[rep(seq_len(n), K), , drop = FALSE]
  long[["response"]] <- factor(rep(levels_response, each = n), levels = levels_response)

  design <- list(n = n, K = K, rt = rt, response = resp_idx,
                 response_levels = levels_response, data = data, long = long,
                 formula = formula, fixed = fixed)

  for (par in c("A", "b", "t0", "s")) {
    design[[par]] <- if (par %in% names(fixed)) {
      fixed_spec(fixed[[par]], n, par)
    } else if (par %in% names(formula)) {
      free_spec(formula[[par]], data, n, par)
    } else {
      NULL
    }
  }

  design[["v"]] <- if ("v" %in% names(fixed)) {
    fixed_spec(fixed[["v"]], n, "v", K = K)
  } else if (inherits(formula[["v"]], "formula")) {
    free_spec(formula[["v"]], long, n * K, "v")
  } else {
    v_spec_list(formula[["v"]], data, n, K, levels_response)
  }

  check_identified(design)
  structure(design, class = c("rdm_design", "list"))
}

## a free parameter: a model matrix with one row per observation
free_spec <- function(form, data, n, par) {
  if (!inherits(form, "formula")) {
    stop("the formula for `", par, "` must be a formula, e.g. ~1.", call. = FALSE)
  }
  X <- stats::model.matrix(form, data)
  if (nrow(X) != n) stop("internal error: wrong number of rows for `", par, "`.",
                         call. = FALSE)
  if (ncol(X) == 0) {
    stop("the formula for `", par, "` has no coefficients (did you use `~ 0`?).",
         call. = FALSE)
  }
  rownames(X) <- NULL
  list(free = TRUE, X = X, coef_names = paste0(par, ":", colnames(X)),
       contrasts = attr(X, "contrasts"), formula = form)
}

## a parameter held constant, recycled to the required length
fixed_spec <- function(value, n, par, K = NULL) {
  if (!is.numeric(value) || anyNA(value) || any(value <= 0)) {
    stop("the fixed value for `", par, "` must be a positive number.", call. = FALSE)
  }
  ## without K the constant is recycled over trials; for the drift rates it can
  ## also be given once per accumulator and is then recycled over trials
  if (is.null(K)) {
    if (length(value) == 1) value <- rep(value, n)
    if (length(value) != n) {
      stop("the fixed value for `", par, "` must have length 1.", call. = FALSE)
    }
  } else if (length(value) == 1) {
    value <- rep(value, n * K)
  } else if (length(value) == K) {
    value <- rep(value, each = n)
  } else {
    stop("the fixed value for `", par, "` must have length 1 or ", K, ".",
         call. = FALSE)
  }
  list(free = FALSE, value = value, coef_names = character(),
       formula = NULL, contrasts = NULL)
}

## drift rates given as one formula per accumulator, stacked accumulator-major
v_spec_list <- function(forms, data, n, K, levels_response) {
  if (!is.list(forms) || length(forms) != K || !all(vapply(forms, inherits, TRUE, "formula"))) {
    stop("`v` must be a single formula or a list of ", K,
         " formulas, one per accumulator.", call. = FALSE)
  }
  blocks <- lapply(seq_len(K), \(j) stats::model.matrix(forms[[j]], data))
  p_j <- vapply(blocks, ncol, 1L)
  if (any(p_j == 0)) {
    stop("every drift-rate formula needs at least one coefficient.", call. = FALSE)
  }
  X <- matrix(0, nrow = n * K, ncol = sum(p_j))
  colnames(X) <- unlist(lapply(seq_len(K), \(j) {
    paste0("v[", levels_response[j], "]:", colnames(blocks[[j]]))
  }))
  start <- cumsum(c(0, p_j))
  for (j in seq_len(K)) {
    rows <- (j - 1) * n + seq_len(n)
    X[rows, start[j] + seq_len(p_j[j])] <- blocks[[j]]
  }
  list(free = TRUE, X = X, coef_names = colnames(X), contrasts = NULL,
       formula = forms, per_accumulator = TRUE)
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
      cat("  ", par, ": fixed at ", unique(spec$value)[1], "\n", sep = "")
    } else {
      cat("  ", par, ": ", ncol(spec$X), " coefficient(s): ",
          paste(spec$coef_names, collapse = ", "), "\n", sep = "")
    }
  }
  invisible(x)
}
