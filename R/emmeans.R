# Support for emmeans.
#
# Every model parameter has its own linear predictor, so emmeans summarises one
# parameter at a time, named with `dpar`:
#
#   emmeans(fit, ~ session, dpar = "v")
#
# `recover_data()` supplies the data of that parameter's reference grid and
# `emm_basis()` rebuilds its model matrix on the grid, using the contrasts that
# were in force when the model was fitted. Terms inside `acc()` are evaluated on
# the grid, so the accumulator (`.acc`) takes part in the grid like any other
# factor.

## the parameter to summarise: `dpar` if given, else the only free parameter
rdm_dpar <- function(object, dpar = NULL) {
  free <- names(object$nll$index$beta)
  free <- free[vapply(free, function(p) length(object$nll$index$beta[[p]]) > 0L,
                      TRUE)]
  if (is.null(dpar)) {
    if (length(free) == 1L) return(free)
    stop("`emmeans()` summarises one parameter at a time: the model estimates ",
         paste(free, collapse = ", "), ". Name the one you want with `dpar`, ",
         "e.g. emmeans(fit, ~ session, dpar = \"v\").", call. = FALSE)
  }
  if (!is.character(dpar) || length(dpar) != 1L || !dpar %in% free) {
    stop("`dpar` must be one of ", paste(free, collapse = ", "), ".",
         call. = FALSE)
  }
  dpar
}

## the variables of the reference grid: the shared terms of the formula, the
## variables inside acc(), and the accumulator itself when acc() is used
rdm_grid_vars <- function(form) {
  labels <- attr(stats::terms(form), "term.labels")
  shared <- character()
  acc_vars <- ".acc"
  has_acc <- FALSE
  for (lab in labels) {
    exprs <- acc_calls(lab)
    if (!length(exprs)) {
      shared <- c(shared, lab)
    } else {
      has_acc <- TRUE
      for (e in exprs) acc_vars <- c(acc_vars, all.vars(str2lang(e)))
    }
  }
  vars <- unique(c(shared, if (has_acc) acc_vars))
  vars
}

#' @exportS3Method emmeans::recover_data
recover_data.rdm_fit <- function(object, dpar = NULL, ...) {
  ## `ref_grid()` replaces errors from `recover_data()` by a generic message, so
  ## an unusable `dpar` is left for `emm_basis()`, where the message survives
  dpar <- tryCatch(rdm_dpar(object, dpar), error = function(e) NULL)
  forms <- lapply(object$design$formula, reformulas::nobars)
  forms <- if (is.null(dpar)) forms else forms[dpar]
  vars <- unique(unlist(lapply(forms, rdm_grid_vars)))
  d <- object$design$long[, vars, drop = FALSE]
  ## emmeans needs at least one variable to build the grid from
  if (!length(vars)) {
    vars <- ".rdm_dummy"
    d[[vars]] <- 1
  }
  ## emmeans builds the grid from the factor levels, so character variables
  ## become factors and the accumulator keeps the order of the responses
  for (v in vars) {
    d[[v]] <- if (v == ".acc") {
      factor(d[[v]], levels = object$design$response_levels)
    } else if (is.character(d[[v]])) {
      factor(d[[v]])
    } else {
      d[[v]]
    }
  }
  attr(d, "terms") <- stats::terms(stats::reformulate(vars), data = d)
  attr(d, "call") <- object$call
  attr(d, "predictors") <- vars
  attr(d, "responses") <- character()
  d
}

#' @exportS3Method emmeans::emm_basis
emm_basis.rdm_fit <- function(object, trms, xlev, grid, dpar = NULL, ...) {
  dpar <- rdm_dpar(object, dpar)
  design <- object$design
  spec <- design[[dpar]]
  form <- reformulas::nobars(spec$formula)
  levels <- design$response_levels
  K <- design$K

  ## the grid has one row per cell; the design wants the data with one row per
  ## trial and accumulator, in the order rdm_design() uses
  if (".acc" %in% names(grid)) {
    wide_vars <- setdiff(names(grid), c(".acc", ".wgt."))
    ## only the accumulator varies, so there is a single trial-level cell
    wide <- if (length(wide_vars)) unique(grid[, wide_vars, drop = FALSE]) else
      grid[1L, wide_vars, drop = FALSE]
    ## the fit used the original types, the grid has factors
    for (nm in intersect(wide_vars, names(design$data))) {
      orig <- design$data[[nm]]
      if (is.character(orig)) {
        wide[[nm]] <- as.character(wide[[nm]])
      } else if (is.numeric(orig) && is.factor(wide[[nm]])) {
        wide[[nm]] <- as.numeric(as.character(wide[[nm]]))
      }
    }
    n <- nrow(wide)
    long <- wide[rep(seq_len(n), K), , drop = FALSE]
    long$.acc <- rep(unique(design$long$.acc), each = n)
    key <- function(x) {
      if (!length(wide_vars)) return(rep("", nrow(x)))
      do.call(paste, c(as.list(x[, wide_vars, drop = FALSE]), list(sep = "\r")))
    }
    rows <- (match(as.character(grid$.acc), levels) - 1L) * n +
      match(key(grid), key(wide))
  } else {
    wide <- grid
    long <- grid
    n <- nrow(grid)
    rows <- seq_len(n)
  }

  ## fits made before the contrasts were stored still carry them on the matrix
  used <- spec$contrasts_used %||% attr(spec$X, "contrasts")
  rebuilt <- par_spec(form, wide, long, n, K, dpar, levels, design$contrasts,
                      check_rank = FALSE, contrasts_used = used)
  missing <- setdiff(spec$coef_names, colnames(rebuilt$X))
  if (length(missing)) {
    stop("the coefficients ", paste(missing, collapse = ", "),
         " could not be rebuilt from the reference grid.", call. = FALSE)
  }
  X <- rebuilt$X[rows, spec$coef_names, drop = FALSE]

  idx <- object$nll$index$beta[[dpar]]
  V <- tryCatch(vcov(object)[idx, idx, drop = FALSE], error = function(e) {
    stop("emmeans needs the standard errors, which are not available: ",
         conditionMessage(e), call. = FALSE)
  })
  list(X = X, bhat = unname(object$beta[idx]),
       nbasis = estimability::all.estble, V = unname(V),
       dffun = function(k, dfargs) Inf, dfargs = list(df = Inf),
       misc = list(tran = "log", inv.lbl = "response", sigma = NA))
}
