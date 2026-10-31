# Design matrices for the RDM.
#
# Every model parameter (A, b, t0, v, s) gets its own model matrix, built from
# a user-supplied formula and the expanded data, which has one row per trial and
# accumulator. Terms in the formula are either shared by all accumulators or
# estimated separately for each accumulator, which is expressed with the
# operator acc().

rdm_pars <- c("A", "b", "t0", "v", "s")

#' Mark a term as a property of the accumulator
#'
#' Used inside the formulas of [fit_rdm()]. `acc(1)` gives one intercept per
#' accumulator, and `acc(x)` estimates the effect of `x` separately for each
#' accumulator. The current accumulator is available inside `acc()` as `.acc`,
#' which allows comparisons with trial variables: `acc(.acc == Word)` is one
#' coefficient that is 1 where the accumulator matches the word.
#'
#' @param x A term of a model formula.
#' @return `x`, unchanged. The operator is interpreted when the design matrices
#'   are built.
#' @export
acc <- function(x) x

#' Build the design matrices of an RDM model
#'
#' Turns one formula per model parameter into a model matrix over the expanded
#' data, that is one row per trial and accumulator. The formula can contain
#' three kinds of terms:
#'
#' * shared terms, such as `~ Condition`, are the same for all accumulators;
#' * `acc(...)` terms are estimated separately for each accumulator;
#' * `acc(...)` terms that refer to the accumulator `.acc`, such as
#'   `acc(.acc == Word)`, are a single predictor that compares the accumulator
#'   with a trial variable.
#'
#' @param data A `data.frame` with one row per trial.
#' @param rt Numeric vector of response times.
#' @param response Vector of observed responses; the sorted unique values define
#'   the accumulators and are coded `1:K`.
#' @param formula Named list of formulas, one per free parameter. Valid names
#'   are `"A"`, `"b"`, `"t0"`, `"v"` and `"s"`, see the details. Random-effect
#'   terms such as `(1 | id)`, `diag(1 + Condition | id)` or `cs(1 + Condition
#'   | id)` may appear in any of the formulas; they are split off here and read
#'   by `rdm_ranef()`. Their left-hand side is built like the fixed part, so
#'   `acc()` and the accumulator variables `response` and `.acc` can be used
#'   there too.
#' @param fixed Named list of parameters that are held constant, on their
#'   natural scale, e.g. `list(s = 1)`. A fixed `b` is the threshold itself, not
#'   the threshold distance. The value can have length 1 (shared) or `K` (one
#'   value per accumulator).
#' @param contrasts Named list of contrast matrices, passed on to
#'   [stats::model.matrix()]. Contrasts of a term inside `acc()` are given under
#'   its full name, e.g. `list("acc(.acc == Word)" = contr.sum)`.
#' @return An object of class `rdm_design`.
#' @keywords internal
#' @noRd
rdm_design <- function(data, rt, response, formula, fixed = list(s = 1),
                       contrasts = list(), anchor = list()) {
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
  not_formula <- !vapply(formula, inherits, TRUE, "formula")
  if (any(not_formula)) {
    stop("the formula for `", names(formula)[which(not_formula)[1]],
         "` must be a formula, e.g. ~1.", call. = FALSE)
  }

  ## expand the data to one row per trial and accumulator
  long <- data[rep(seq_len(n), K), , drop = FALSE]
  long[["response"]] <- factor(rep(levels_response, each = n), levels = levels_response)
  ## the accumulator itself, in the type of the response, for use inside acc()
  long[[".acc"]] <- rep(if (is.numeric(response)) as.numeric(levels_response) else
                          levels_response, each = n)

  ## the random-effect terms are the `(term | group)` parts of the formulas;
  ## they are kept on the side, the remaining fixed part builds the model matrix
  ranef_terms <- lapply(formula, rdm_ranef_terms)
  ranef_terms <- ranef_terms[vapply(ranef_terms, length, 1L) > 0L]

  design <- list(n = n, K = K, rt = rt, response = resp_idx,
                 response_levels = levels_response, data = data, long = long,
                 formula = formula, fixed = fixed, contrasts = contrasts,
                 ranef_terms = ranef_terms)

  for (par in rdm_pars) {
    design[[par]] <- if (par %in% names(fixed)) {
      fixed_spec(fixed[[par]], n, par, K)
    } else if (par %in% names(formula)) {
      par_spec(reformulas::nobars(formula[[par]]), data, long, n, K, par,
               levels_response, contrasts)
    } else {
      NULL
    }
  }

  design <- anchor_s(design, anchor)
  check_identified(design)
  structure(design, class = c("rdm_design", "list"))
}

## the design of one parameter: shared terms, terms estimated per accumulator
## and terms that compare the accumulator with a trial variable. `what` names
## the term in the error messages; `rdm_ranef()` uses it for the random-effect
## terms, which go through the same machinery.
par_spec <- function(form, data, long, n, K, par, levels, contrasts, what = par) {
  if (!inherits(form, "formula")) {
    stop("the formula for `", what, "` must be a formula, e.g. ~1.", call. = FALSE)
  }
  env <- environment(form)
  if (is.null(env)) env <- parent.frame()
  tt <- stats::terms(form)
  labels <- attr(tt, "term.labels")
  new_labels <- labels
  is_stacked <- rep(FALSE, length(labels))
  stacked <- list()
  mapping <- list()

  for (i in seq_along(labels)) {
    exprs <- acc_calls(labels[i])
    if (!length(exprs)) next
    if (length(exprs) > 1L) {
      stop("only one acc() term per formula term is supported, but `", labels[i],
           "` has ", length(exprs), ".", call. = FALSE)
    }
    expr <- exprs[1]
    relates <- grepl("\\.acc", expr) || grepl("\\bresponse\\b", expr)
    if (relates) {
      ## the same comparison may appear in several terms, e.g. alone and in an
      ## interaction, and then it is the same variable
      if (is.null(mapping[[expr]])) {
        mapping[[expr]] <- paste0(".", par, "_acc", length(mapping) + 1L)
      }
      new_labels[i] <- acc_map(labels[i], \(e) mapping[[e]])
    } else if (identical(labels[i], paste0("acc(", expr, ")"))) {
      is_stacked[i] <- TRUE
      stacked[[length(stacked) + 1L]] <- list(expr = expr, label = labels[i])
    } else {
      stop("the term `", labels[i], "` combines acc() with other terms, which ",
           "is not supported. Write the acc() term on its own; to let an effect ",
           "vary over the accumulators write it inside acc(), e.g. ",
           "`acc(Condition)` instead of `Condition:acc(1)`.", call. = FALSE)
    }
  }

  ## the accumulator-relative terms are evaluated on the expanded data
  for (expr in names(mapping)) {
    val <- eval(parse(text = expr), long, env)
    ## a logical becomes a 0/1 indicator; passed on as a logical it would be
    ## turned into a factor by model.matrix, which codes only the first such
    ## term as full dummies and the rest against it
    if (is.logical(val)) val <- as.numeric(val)
    if (length(val) == 1L) val <- rep(val, n * K)
    if (length(val) != n * K) {
      stop("the term `acc(", expr, ")` must have one value per trial and ",
           "accumulator.", call. = FALSE)
    }
    long[[mapping[[expr]]]] <- val
  }

  ## the intercept of the parameter is dropped when an acc() term brings its own
  inner <- lapply(stacked, \(s) stats::terms(as.formula(paste("~", s$expr), env = env)))
  has_inner_intercept <- any(vapply(inner, \(x) attr(x, "intercept") == 1L, TRUE))
  intercept <- attr(tt, "intercept") == 1L && !has_inner_intercept

  shared <- reformulate(new_labels[!is_stacked], intercept = intercept, env = env)
  X <- stats::model.matrix(shared, long, contrasts.arg = map_contrasts(contrasts, mapping))
  if (ncol(X)) {
    colnames(X) <- paste0(par, ":", pretty_names(colnames(X), mapping))
  }

  ## terms estimated separately for each accumulator: the design of the term is
  ## evaluated once per trial and replicated, with its own coefficients per
  ## accumulator
  if (length(stacked)) {
    blocks <- lapply(stacked, function(s) {
      Xs <- stats::model.matrix(as.formula(paste("~", s$expr), env = env), data,
                                contrasts.arg = map_contrasts(contrasts, list(),
                                                              s$label, s$expr))
      p <- ncol(Xs)
      out <- matrix(0, n * K, K * p)
      cn <- character(K * p)
      for (j in seq_len(K)) {
        out[(j - 1) * n + seq_len(n), (j - 1) * p + seq_len(p)] <- Xs
        cn[(j - 1) * p + seq_len(p)] <- paste0(par, "[", levels[j], "]:", colnames(Xs))
      }
      colnames(out) <- cn
      out
    })
    X <- cbind(X, do.call(cbind, blocks))
  }

  rownames(X) <- NULL
  if (ncol(X) == 0) {
    stop("the formula for `", what, "` has no coefficients (did you use `~ 0`?).",
         call. = FALSE)
  }
  rank <- qr(X)$rank
  if (rank < ncol(X)) {
    stop("the model matrix of `", what, "` is rank deficient (", rank, " of ",
         ncol(X), " coefficients are identified): ", paste(colnames(X), collapse = ", "),
         ". Remove the redundant terms.", call. = FALSE)
  }

  relates <- grepl("\\.acc", labels) | grepl("\\bresponse\\b", labels)
  list(free = TRUE, X = X, coef_names = colnames(X), formula = form,
       per_accumulator = length(stacked) > 0L || any(relates))
}

## a parameter held constant: either one value for all accumulators or one value
## per accumulator, recycled to the required length
fixed_spec <- function(value, n, par, K) {
  if (!is.numeric(value) || anyNA(value) || any(value <= 0)) {
    stop("the fixed value for `", par, "` must be a positive number.", call. = FALSE)
  }
  if (length(value) == 1) {
    value <- rep(value, n)
    per_accumulator <- FALSE
  } else if (length(value) == K) {
    value <- rep(value, each = n)
    per_accumulator <- TRUE
  } else {
    stop("the fixed value for `", par, "` must have length 1 or ", K, ".",
         call. = FALSE)
  }
  list(free = FALSE, value = value, coef_names = character(), formula = NULL,
       per_accumulator = per_accumulator)
}

## the expressions inside the acc() calls of one term label
acc_calls <- function(label) {
  out <- character()
  start <- 1L
  repeat {
    hit <- acc_find(label, start)
    if (is.null(hit)) break
    out <- c(out, substr(label, hit$open + 1L, hit$close - 1L))
    start <- hit$close + 1L
    if (start > nchar(label)) break
  }
  out
}

## replace every acc() call in a term label by the value returned from `fun`
acc_map <- function(label, fun) {
  out <- label
  start <- 1L
  repeat {
    hit <- acc_find(out, start)
    if (is.null(hit)) break
    expr <- substr(out, hit$open + 1L, hit$close - 1L)
    repl <- fun(expr)
    before <- substr(out, 1L, hit$open - 4L)
    out <- paste0(before, repl, substr(out, hit$close + 1L, nchar(out)))
    start <- nchar(before) + nchar(repl) + 1L
    if (start > nchar(out)) break
  }
  out
}

## position of the parentheses of the next acc() call at or after `start`
acc_find <- function(label, start = 1L) {
  if (start > nchar(label)) return(NULL)
  m <- regexpr("(^|[^[:alnum:]_.])acc\\(", substring(label, start))
  if (m < 0) return(NULL)
  open <- start + m + attr(m, "match.length") - 2L
  depth <- 0L
  i <- open
  repeat {
    ch <- substr(label, i, i)
    if (ch == "(") depth <- depth + 1L
    if (ch == ")") {
      depth <- depth - 1L
      if (depth == 0L) break
    }
    i <- i + 1L
    if (i > nchar(label)) {
      stop("unbalanced parentheses in the term `", label, "`.", call. = FALSE)
    }
  }
  list(open = open, close = i)
}

## contrasts for terms inside acc() are given under the name of the whole term,
## e.g. "acc(.acc == Word)"; they are stored under the variable name
map_contrasts <- function(contrasts, mapping, label = NULL, expr = NULL) {
  if (!length(contrasts)) return(NULL)
  out <- contrasts
  nm <- names(out)
  for (e in names(mapping)) {
    hit <- nm == paste0("acc(", e, ")")
    nm[hit] <- mapping[[e]]
  }
  if (!is.null(expr)) {
    hit <- nm == paste0("acc(", expr, ")")
    nm[hit] <- expr
  }
  names(out) <- nm
  out
}

## replace the internal variable names of the accumulator terms by the term the
## user wrote
pretty_names <- function(names, mapping) {
  for (expr in names(mapping)) {
    names <- gsub(mapping[[expr]], paste0("acc(", expr, ")"), names, fixed = TRUE)
  }
  names
}

## fix the diffusion constant of one accumulator at 1, which pins down the scale
## of the model. With `s = ~ acc(1)` every accumulator is free, and one of them
## is anchored; `anchor` chooses which one. A formula that leaves an accumulator
## out already anchors that one and is left alone.
anchor_s <- function(design, anchor = list()) {
  unknown <- setdiff(names(anchor), "s")
  if (length(unknown)) {
    stop("`anchor` only applies to `s`, but it names ",
         paste(unknown, collapse = ", "), ".", call. = FALSE)
  }
  wanted <- as.character(anchor$s)
  if (length(wanted) > 1L) {
    stop("`anchor$s` must name a single accumulator.", call. = FALSE)
  }
  spec <- design$s
  if (is.null(spec) || !isTRUE(spec$free)) {
    if (length(wanted)) {
      stop("`anchor$s` was given, but `s` has no formula: it is already fixed. ",
           "Give `s` a formula, e.g. s = ~ acc(1).", call. = FALSE)
    }
    return(design)
  }
  ## a single `s` shared by all accumulators is the scale of the whole model,
  ## which can only be pinned by a fixed A, b or v; check_identified() does that
  if (!isTRUE(spec$per_accumulator)) {
    if (length(wanted)) {
      stop("`anchor$s` was given, but `s` is one value shared by all ",
           "accumulators, so there is no accumulator to anchor. Write ",
           "s = ~ acc(1) to estimate one per accumulator.", call. = FALSE)
    }
    return(design)
  }
  acc <- as.character(design$long$.acc)
  rows <- split(seq_along(acc), acc)
  zero <- vapply(rows, \(i) all(spec$X[i, , drop = FALSE] == 0), TRUE)
  if (any(zero)) {
    if (length(wanted)) {
      stop("`anchor$s` was given, but the formula of `s` already fixes the ",
           "diffusion constant of ", paste(names(zero)[zero], collapse = ", "),
           " at 1.", call. = FALSE)
    }
    return(design)
  }
  ## a single anchor is only enough when some other parameter ties the
  ## accumulators together; otherwise every s is absorbed into its own A, b, v
  per_acc <- vapply(c("A", "b", "v"),
                    \(p) isTRUE(design[[p]]$per_accumulator), TRUE)
  if (all(per_acc)) {
    stop("`s` cannot be estimated per accumulator: A, b and v all vary over the ",
         "accumulators, so each diffusion constant is absorbed into the scale of ",
         "its own accumulator. Hold at least one of A, b and v constant over the ",
         "accumulators.", call. = FALSE)
  }
  levs <- design$response_levels
  ## only the per-accumulator intercepts of acc(1) are single columns that can
  ## be dropped; anything else would need a linear constraint on the coefficients
  if (ncol(spec$X) != design$K || !all(abs(rowSums(spec$X) - 1) < 1e-12)) {
    stop("cannot tell which accumulator to anchor in the formula of `s`. ",
         "Write it as `s = ~ acc(1)` to anchor an accumulator automatically, or ",
         "leave one out, e.g. s = ~ 0 + acc(.acc != \"", levs[1], "\").",
         call. = FALSE)
  }
  lvl <- if (length(wanted)) wanted else levs[1]
  if (!lvl %in% levs) {
    stop("`anchor$s` must be one of the accumulators: ",
         paste(levs, collapse = ", "), ".", call. = FALSE)
  }
  drop <- match(lvl, levs)
  design$s$anchor <- list(level = lvl, name = spec$coef_names[drop], value = 1)
  design$s$X <- spec$X[, -drop, drop = FALSE]
  design$s$coef_names <- spec$coef_names[-drop]
  design
}

## the diffusion constant s is only identified if the scale of the model is
## pinned down, either by fixing one of A, b and v or by leaving one accumulator
## out of the formula of s, which fixes its diffusion constant at 1
check_identified <- function(design) {
  if (!isTRUE(design$s$free)) return(invisible(TRUE))
  if (length(intersect(names(design$fixed), c("A", "b", "v")))) {
    return(invisible(TRUE))
  }
  acc <- as.character(design$long$.acc)
  pinned <- vapply(split(seq_along(acc), acc), \(rows) {
    all(design$s$X[rows, , drop = FALSE] == 0)
  }, TRUE)
  if (any(pinned)) return(invisible(TRUE))
  stop("a free `s` (the diffusion constant) is not identified: `s` scales ",
       "`A`, `b` and `v`, so the scale of the model must be pinned down. Either ",
       "fix one of A, b and v, e.g. fixed = list(A = 0.5), or leave one ",
       "accumulator out of the formula of `s`, e.g. ",
       "s = ~ 0 + acc(.acc != \"red\"), which fixes its value at 1.",
       call. = FALSE)
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
      if (spec$per_accumulator) {
        cat("  ", par, ": fixed per accumulator at ",
            paste(unique(spec$value), collapse = ", "), "\n", sep = "")
      } else {
        cat("  ", par, ": fixed at ", unique(spec$value)[1], "\n", sep = "")
      }
    } else {
      cat("  ", par, ": ", ncol(spec$X), " coefficient(s)",
          if (spec$per_accumulator) ", varying over the accumulators" else "",
          ": ", paste(spec$coef_names, collapse = ", "), "\n", sep = "")
      if (!is.null(spec$anchor)) {
        cat("  ", par, ": anchored at 1 for ", spec$anchor$level,
            " (its coefficient is not estimated)\n", sep = "")
      }
    }
  }
  invisible(x)
}
