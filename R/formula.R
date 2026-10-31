# The user-facing formula constructor. One formula per model parameter, with
# the random effects written inside the formulas, as in brms::bf().

#' Model formulas for the RDM
#'
#' Collects one formula per model parameter and is the `formula` argument of
#' [fit_rdm()]. Every formula has two parts:
#'
#' * the fixed effects, which may contain [acc()] terms, e.g.
#'   `v = ~ Condition + acc(.acc == Word)`;
#' * the random effects, written in the lme4 style (as in brms), e.g.
#'   `(1 | id)`, `(1 + Condition | id)` or `(1 | id) + (1 | item)`;
#' * a term can be wrapped in a covariance structure, as in glmmTMB:
#'   `diag()` gives one variance per column and no correlations, `homdiag()`
#'   gives a single variance shared by all columns, and `cs()` (compound
#'   symmetry) gives a single variance and a single correlation between any two
#'   columns. For example, `diag(0 + Word | id)` estimates one variance per
#'   level of a factor, where `(0 + Word | id)` would estimate a covariance
#'   between the levels. The correlation of a `cs()` term is constrained to the
#'   values that keep the covariance matrix positive definite and is held away
#'   from the singular boundary by the prior on the covariance parameters.
#'
#' A parameter that is left out takes its default: no start-point variability
#' (`A = 0`), a threshold one unit above the start point (`b = A + 1`) and no
#' drift-rate bias. `v` is required, because the drift rates are what drives
#' the race. `s` is fixed to 1 by [fit_rdm()] unless a formula is given for it.
#' A free `s` needs the scale of the model pinned down: with `s = ~ acc(1)` the
#' diffusion constant of the first response level is held at 1 and the others
#' are estimated relative to it (see the `anchor` argument of [fit_rdm()]),
#' and a formula that leaves an accumulator out, e.g.
#' `s = ~ 0 + acc(.acc != "Word")`, fixes that one instead.
#'
#' @param A,b,t0,v,s One-sided formulas named after the model parameters, or
#'   `NULL` for a parameter that is left out of the model. A two-sided formula
#'   such as `y ~ 1` is accepted and only its right-hand side is used.
#' @return An object of class `rdm_formula`, a named list of formulas.
#' @seealso [fit_rdm()], [acc()]
#' @export
rdm_formula <- function(A = NULL, b = NULL, t0 = NULL, v = NULL, s = NULL) {
  args <- list(A = A, b = b, t0 = t0, v = v, s = s)
  as_rdm_formula(args[!vapply(args, is.null, TRUE)])
}

## normalise the formula argument of fit_rdm(): an `rdm_formula`, a named list
## of formulas, or nothing
as_rdm_formula <- function(x) {
  if (inherits(x, "rdm_formula")) return(x)
  if (!is.list(x) || is.null(names(x)) || any(names(x) == "")) {
    stop("`formula` must be an `rdm_formula()` or a named list of one formula ",
         "per parameter, e.g. list(A = ~1, v = ~1).", call. = FALSE)
  }
  unknown <- setdiff(names(x), rdm_pars)
  if (length(unknown)) {
    stop("unknown parameter(s) in `formula`: ", paste(unknown, collapse = ", "),
         ". Valid names are ", paste(rdm_pars, collapse = ", "), ".", call. = FALSE)
  }
  for (par in names(x)) {
    form <- x[[par]]
    if (!inherits(form, "formula")) {
      stop("the formula for `", par, "` must be a formula, e.g. ~1.", call. = FALSE)
    }
    ## a response on the left-hand side is dropped, as in brms
    if (length(form) == 3L) {
      env <- environment(form)
      if (is.null(env)) env <- parent.frame()
      form <- eval(call("~", form[[3]]), env)
    }
    x[[par]] <- form
  }
  structure(x[intersect(rdm_pars, names(x))], class = c("rdm_formula", "list"))
}

#' @export
print.rdm_formula <- function(x, ...) {
  cat("<rdm_formula>\n")
  for (par in names(x)) {
    fixed <- reformulas::nobars(x[[par]])
    terms <- rdm_ranef_terms(x[[par]])
    cat("  ", par, " ~ ", deparse1(fixed[[length(fixed)]]), sep = "")
    if (length(terms)) {
      shown <- vapply(terms, function(t) {
        bar <- paste0("(", deparse1(t$term), ")")
        if (identical(t$struct, "us")) bar else paste0(t$struct, bar)
      }, "")
      cat(" + ", paste(shown, collapse = " + "), sep = "")
    }
    cat("\n")
  }
  invisible(x)
}
