test_that("rdm_formula collects one formula per parameter", {
  f <- rdm_formula(A = ~1, b = ~cond, t0 = ~1, v = ~acc(1))
  expect_s3_class(f, "rdm_formula")
  expect_identical(names(f), c("A", "b", "t0", "v"))
  expect_true(all(vapply(f, inherits, TRUE, "formula")))

  ## parameters that are left out are dropped
  f <- rdm_formula(v = ~1, s = ~1)
  expect_identical(names(f), c("v", "s"))
  expect_null(f$A)

  ## the order of the arguments does not matter
  expect_identical(names(rdm_formula(v = ~1, A = ~1)), c("A", "v"))
})

test_that("random-effect terms survive the constructor", {
  f <- rdm_formula(A = ~1 + (1 | id), v = ~cond + (1 + cond | id))
  expect_length(reformulas::findbars(f$A), 1L)
  expect_length(reformulas::findbars(f$v), 1L)
  expect_identical(reformulas::nobars(f$A)[[2]], 1)
})

test_that("a two-sided formula keeps only its right-hand side", {
  f <- rdm_formula(v = response ~ cond)
  expect_identical(f$v[[2]], as.name("cond"))
  expect_identical(attr(stats::terms(f$v), "term.labels"), "cond")
})

test_that("formulas are validated", {
  expect_error(rdm_formula(v = "cond"), "must be a formula")
  expect_error(rdm_formula(w = ~1), "unused argument")
  expect_error(as_rdm_formula(list(~1)), "named list")
  expect_error(as_rdm_formula(list(w = ~1)), "unknown parameter")
})

test_that("the print method shows the fixed and random parts", {
  out <- capture.output(print(rdm_formula(A = ~1 + (1 | id), v = ~cond)))
  expect_identical(out[1], "<rdm_formula>")
  expect_match(out[2], "^  A ~ 1 \\+ \\(1 \\| id\\)$")
  expect_match(out[3], "^  v ~ cond$")

  out <- capture.output(print(rdm_formula(v = ~ acc(1) + diag(0 + cond | id))))
  expect_match(out[2], "diag\\(0 \\+ cond \\| id\\)")

  out <- capture.output(print(rdm_formula(v = ~ homdiag(1 | id) + cs(1 | id))))
  expect_match(out[2], "homdiag\\(1 \\| id\\)")
  expect_match(out[2], "cs\\(1 \\| id\\)")
})
