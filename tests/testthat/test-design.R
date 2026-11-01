set.seed(1)
example_data <- data.frame(
  rt = runif(6, 0.5, 1.5),
  resp = rep(c("a", "b"), 3),
  cond = factor(rep(c("x", "y"), each = 3)),
  id = rep(1:3, each = 2)
)

test_that("a design has one model matrix per parameter", {
  d <- rdm_design(example_data, example_data$rt, example_data$resp,
                  formula = list(A = ~1, b = ~cond, t0 = ~1, v = ~1))
  expect_s3_class(d, "rdm_design")
  expect_identical(d$n, 6L)
  expect_identical(d$K, 2L)
  expect_identical(d$response_levels, c("a", "b"))
  expect_identical(d$response, c(1L, 2L, 1L, 2L, 1L, 2L))
  ## the matrices have one row per trial and accumulator
  expect_identical(dim(d$A$X), c(12L, 1L))
  expect_identical(dim(d$b$X), c(12L, 2L))
  expect_identical(d$b$coef_names, c("b:(Intercept)", "b:condy"))
  ## a shared drift formula gives one coefficient for all accumulators
  expect_identical(d$v$coef_names, "v:(Intercept)")
  expect_false(d$v$per_accumulator)
  ## s is fixed at 1 by default
  expect_false(d$s$free)
  expect_identical(unique(d$s$value), 1)
})

test_that("rt and response may be given by column name", {
  f <- list(A = ~1, b = ~cond, t0 = ~1, v = ~1)
  by_name <- rdm_design(example_data, rt = "rt", response = "resp", formula = f)
  by_value <- rdm_design(example_data, example_data$rt, example_data$resp, formula = f)
  expect_identical(by_name$rt, by_value$rt)
  expect_identical(by_name$response, by_value$response)
  ## the same for fit_rdm(), which passes the names on
  expect_error(fit_rdm(example_data, rt = "rrt", response = "resp", formula = f),
               "not a column of `data`")
})

test_that("acc(1) gives one intercept per accumulator", {
  d <- rdm_design(example_data, example_data$rt, example_data$resp,
                  formula = list(A = ~1, t0 = ~1, b = ~acc(1), v = ~1))
  expect_true(d$b$per_accumulator)
  expect_identical(d$b$coef_names, c("b[a]:(Intercept)", "b[b]:(Intercept)"))
  expect_identical(dim(d$b$X), c(12L, 2L))
  ## rows 1:6 belong to accumulator a, rows 7:12 to accumulator b
  expect_equal(d$b$X[1:6, 1], rep(1, 6))
  expect_equal(d$b$X[7:12, 2], rep(1, 6))
  expect_equal(d$b$X[1:6, 2], rep(0, 6))
  expect_equal(d$b$X[7:12, 1], rep(0, 6))
  ## the shared intercept is dropped, the acc() term brings its own
  expect_false(d$A$per_accumulator)
  expect_false(d$s$free)
})

test_that("acc() can contain a term with several coefficients", {
  d <- rdm_design(example_data, example_data$rt, example_data$resp,
                  formula = list(A = ~1, t0 = ~1, v = ~acc(cond)))
  expect_identical(d$v$coef_names,
                   c("v[a]:(Intercept)", "v[a]:condy",
                     "v[b]:(Intercept)", "v[b]:condy"))
  expect_identical(dim(d$v$X), c(12L, 4L))
  expect_equal(d$v$X[7:12, 4], c(0, 0, 0, 1, 1, 1))
  expect_equal(d$v$X[1:6, 4], rep(0, 6))
})

test_that("a drift formula can describe the accumulator via `response`", {
  d <- rdm_design(example_data, example_data$rt, example_data$resp,
                  formula = list(A = ~1, t0 = ~1, v = ~ 0 + response))
  expect_identical(d$v$coef_names, c("v:responsea", "v:responseb"))
  ## exactly one accumulator per row is active
  expect_true(all(rowSums(d$v$X) == 1))
  expect_equal(d$v$X[1:6, 1], rep(1, 6))
  expect_equal(d$v$X[7:12, 2], rep(1, 6))
})

test_that("acc(1) and `response` describe the same linear predictor", {
  d1 <- rdm_design(example_data, example_data$rt, example_data$resp,
                   formula = list(A = ~1, t0 = ~1, v = ~acc(1)))
  d2 <- rdm_design(example_data, example_data$rt, example_data$resp,
                   formula = list(A = ~1, t0 = ~1, v = ~ 0 + response))
  expect_equal(as.numeric(d1$v$X %*% c(2, 3)),
               as.numeric(d2$v$X %*% c(2, 3)))
})

test_that("a term that refers to .acc is a single predictor", {
  d <- rdm_design(example_data, example_data$rt, example_data$resp,
                  formula = list(A = ~1, t0 = ~1, v = ~ acc(.acc == resp)))
  expect_identical(d$v$coef_names, c("v:(Intercept)", "v:acc(.acc == resp)"))
  ## row 1 is a response "a" from accumulator a, row 2 a "b" from accumulator a
  expect_equal(unname(d$v$X[1, 2]), 1)
  expect_equal(unname(d$v$X[2, 2]), 0)
  expect_equal(sum(d$v$X[, 2]), 6)
})

test_that("a comparison inside acc() is a 0/1 indicator", {
  ## a logical is coded 0/1 and not, as model.matrix would, as a factor: with
  ## two such terms the second would otherwise be coded against the first
  d <- rdm_design(example_data, example_data$rt, example_data$resp,
                  formula = list(A = ~1, t0 = ~1,
                                 v = ~ 0 + acc(.acc == resp) + acc(.acc == "a")))
  expect_identical(d$v$coef_names,
                   c("v:acc(.acc == resp)", "v:acc(.acc == \"a\")"))
  expect_setequal(unique(d$v$X[, 1]), c(0, 1))
  expect_setequal(unique(d$v$X[, 2]), c(0, 1))
  ## the coding can be written out, with the opposite sign for the match
  d2 <- rdm_design(example_data, example_data$rt, example_data$resp,
                   formula = list(A = ~1, t0 = ~1,
                                  v = ~ acc(2 * (.acc == resp) - 1)))
  expect_equal(unname(d2$v$X[, 2]), unname(2 * d$v$X[, 1] - 1))
})

test_that("a free s is anchored at one accumulator", {
  levs <- sort(unique(as.character(example_data$resp)))
  fml <- list(A = ~1, t0 = ~1, v = ~ acc(1), s = ~ acc(1))
  d <- rdm_design(example_data, example_data$rt, example_data$resp, formula = fml,
                  fixed = list())
  ## the first accumulator is held at 1, the others are estimated relative to it
  expect_identical(d$s$anchor$level, levs[1])
  expect_identical(d$s$coef_names, paste0("s[", levs[2], "]:(Intercept)"))
  s <- exp(drop(d$s$X %*% 0.3))
  expect_equal(unique(s[d$long$.acc == levs[1]]), 1)
  expect_equal(unique(s[d$long$.acc == levs[2]]), exp(0.3))

  ## another accumulator can be anchored instead
  d2 <- rdm_design(example_data, example_data$rt, example_data$resp, formula = fml,
                   fixed = list(), anchor = list(s = levs[2]))
  expect_identical(d2$s$anchor$level, levs[2])
  expect_identical(d2$s$coef_names, paste0("s[", levs[1], "]:(Intercept)"))
  s2 <- exp(drop(d2$s$X %*% -0.4))
  expect_equal(unique(s2[d2$long$.acc == levs[2]]), 1)
  expect_equal(unique(s2[d2$long$.acc == levs[1]]), exp(-0.4))

  ## leaving an accumulator out of the formula fixes that one
  d3 <- rdm_design(example_data, example_data$rt, example_data$resp,
                   formula = list(A = ~1, t0 = ~1, v = ~ acc(1),
                                  s = as.formula(paste0(
                                    "~ 0 + acc(.acc == \"", levs[1], "\")"))),
                   fixed = list())
  expect_null(d3$s$anchor)
  expect_identical(d3$s$coef_names, paste0("s:acc(.acc == \"", levs[1], "\")"))
  s3 <- exp(drop(d3$s$X %*% 0.7))
  expect_equal(unique(s3[d3$long$.acc == levs[1]]), exp(0.7))
  expect_equal(unique(s3[d3$long$.acc != levs[1]]), 1)
})

test_that("anchoring s fails when no parameter pins the scale", {
  ## A, b and v all vary over the accumulators, so every s is absorbed
  fml <- list(A = ~ acc(1), t0 = ~1, b = ~ acc(1), v = ~ acc(1), s = ~ acc(1))
  expect_error(rdm_design(example_data, example_data$rt, example_data$resp,
                          formula = fml, fixed = list()),
               "cannot be estimated per accumulator")
  ## anchor for a parameter that is not s, or for an s that is already fixed
  fml <- list(A = ~1, t0 = ~1, v = ~1, s = ~ acc(1))
  expect_error(rdm_design(example_data, example_data$rt, example_data$resp,
                          formula = fml, fixed = list(), anchor = list(v = "a")),
               "only applies to `s`")
  expect_error(rdm_design(example_data, example_data$rt, example_data$resp,
                          formula = list(A = ~1, t0 = ~1, v = ~1),
                          fixed = list(), anchor = list(s = "a")),
               "already fixed")
  ## an s formula that the anchor cannot be dropped from
  expect_error(rdm_design(example_data, example_data$rt, example_data$resp,
                          formula = list(A = ~1, t0 = ~1, v = ~1,
                                         s = ~ 1 + acc(.acc == "a")),
                          fixed = list()),
               "cannot tell which accumulator")
  ## a single s shared by all accumulators has no accumulator to anchor
  expect_error(rdm_design(example_data, example_data$rt, example_data$resp,
                          formula = list(A = ~1, t0 = ~1, v = ~1, s = ~1),
                          fixed = list(), anchor = list(s = "a")),
               "shared by all accumulators")
})

test_that("an accumulator-relative term can be interacted with other terms", {
  d <- rdm_design(example_data, example_data$rt, example_data$resp,
                  formula = list(A = ~1, t0 = ~1,
                                 v = ~ cond * acc(.acc == resp)))
  expect_identical(d$v$coef_names,
                   c("v:(Intercept)", "v:condy", "v:acc(.acc == resp)",
                     "v:condy:acc(.acc == resp)"))
})

test_that("redundant and unsupported acc() terms are rejected", {
  expect_error(rdm_design(example_data, example_data$rt, example_data$resp,
                          formula = list(A = ~1, t0 = ~1,
                                         v = ~ cond + acc(cond))),
               "rank deficient")
  expect_error(rdm_design(example_data, example_data$rt, example_data$resp,
                          formula = list(A = ~1, t0 = ~1, v = ~ cond:acc(1))),
               "combines acc")
})

test_that("fixed parameters are stored as constants", {
  d <- rdm_design(example_data, example_data$rt, example_data$resp,
                  formula = list(A = ~1, b = ~1, t0 = ~1, v = ~1),
                  fixed = list(s = 2))
  expect_false(d$s$free)
  expect_identical(unique(d$s$value), 2)
  ## a constant drift rate is recycled over accumulators
  d <- rdm_design(example_data, example_data$rt, example_data$resp,
                  formula = list(A = ~1, b = ~1, t0 = ~1),
                  fixed = list(s = 1, v = c(2, 1)))
  expect_false(d$v$free)
  expect_identical(length(d$v$value), 12L)
})

test_that("the start-point variability can be fixed at zero", {
  d <- rdm_design(example_data, example_data$rt, example_data$resp,
                  formula = list(b = ~1, t0 = ~1, v = ~1),
                  fixed = list(A = 0))
  expect_false(d$A$free)
  expect_identical(unique(d$A$value), 0)
  ## only for all accumulators, and no other parameter may be zero
  expect_error(
    rdm_design(example_data, example_data$rt, example_data$resp,
               formula = list(t0 = ~1, v = ~1),
               fixed = list(A = c(0, 0.2))),
    "zero for all accumulators"
  )
  expect_error(
    rdm_design(example_data, example_data$rt, example_data$resp,
               formula = list(b = ~1, v = ~1),
               fixed = list(t0 = 0)),
    "positive number"
  )
})

test_that("a fixed parameter can take one value per accumulator", {
  d <- rdm_design(example_data, example_data$rt, example_data$resp,
                  formula = list(A = ~1, t0 = ~1, v = ~1),
                  fixed = list(s = 1, b = c(1.2, 0.8)))
  expect_true(d$b$per_accumulator)
  expect_equal(d$b$value, rep(c(1.2, 0.8), each = 6))
})

test_that("a free diffusion constant needs another parameter to fix the scale", {
  expect_error(
    rdm_design(example_data, example_data$rt, example_data$resp,
               formula = list(A = ~1, b = ~1, t0 = ~1, v = ~1, s = ~1),
               fixed = list()),
    "not identified"
  )
  d <- rdm_design(example_data, example_data$rt, example_data$resp,
                  formula = list(b = ~1, t0 = ~1, v = ~1, s = ~1),
                  fixed = list(A = 0.5))
  expect_true(d$s$free)
  expect_true(d$A$free == FALSE)
})

test_that("the design validates its input", {
  expect_error(rdm_design(example_data, example_data$rt, example_data$resp,
                          formula = list(A = ~1, t0 = ~1)),
               "drift rates")
  expect_error(rdm_design(example_data, example_data$rt, example_data$resp,
                          formula = list(foo = ~1, v = ~1)),
               "unknown parameter")
  expect_error(rdm_design(example_data, example_data$rt, example_data$resp,
                          formula = list(A = ~1, v = ~1), fixed = list(A = 1)),
               "both `formula` and `fixed`")
  expect_error(rdm_design(example_data, example_data$rt, example_data$resp,
                          formula = list(A = ~ 0, t0 = ~1, v = ~1)),
               "no coefficients")
  expect_error(rdm_design(example_data, example_data$rt, example_data$resp[1:5],
                          formula = list(A = ~1, t0 = ~1, v = ~1)),
               "one value per row")
  expect_error(rdm_design(example_data, example_data$rt, example_data$resp,
                          formula = list(A = ~1, t0 = ~1, v = list(~1))),
               "must be a formula")
})

test_that("the print method summarises the design", {
  d <- rdm_design(example_data, example_data$rt, example_data$resp,
                  formula = list(A = ~1, b = ~cond, t0 = ~1, v = ~acc(1)))
  out <- capture.output(print(d))
  expect_true(any(grepl("6 trials, 2 accumulators", out)))
  expect_true(any(grepl("b: 2 coefficient", out)))
  expect_true(any(grepl("v: 2 coefficient.*varying over the accumulators", out)))
})

test_that("the print method shows the unique rows of each model matrix", {
  d <- rdm_design(example_data, example_data$rt, example_data$resp,
                  formula = list(A = ~1, b = ~cond, t0 = ~1, v = ~acc(1)))
  out <- capture.output(print(d))
  ## the variables of the formula label the rows, next to the coefficients
  expect_true(any(grepl("cond.*condy", out)))
  expect_true(any(grepl("acc.*v\\[a\\]:\\(Intercept\\)", out)))
  ## only the unique rows are printed: two rows for b, not six trials
  b_rows <- out[seq(grep("b: 2 coefficient", out) + 1L, length.out = 3L)]
  expect_identical(sum(grepl("^\\s+[xy]\\s", b_rows)), 2L)
  ## and at most `max_unique` of them
  out1 <- capture.output(print(d, max_unique = 1))
  expect_true(any(grepl("more of 2 unique rows", out1)))
  expect_error(print(d, max_unique = 0), "positive number")
})

test_that("the print method labels acc() terms with their variables", {
  d <- rdm_design(example_data, example_data$rt, example_data$resp,
                  formula = list(A = ~1, b = ~ 1 + acc(.acc == "a"), t0 = ~1,
                                 v = ~ acc(cond)))
  out <- capture.output(print(d))
  ## the comparison is a variable of the design, without the acc() wrapper
  expect_true(any(grepl("\\.acc == \"a\" \\(Intercept\\)", out)))
  ## and the accumulator comes last, so that it varies fastest
  expect_true(any(grepl("^\\s+cond\\s+acc\\s+v\\[a\\]", out)))
  expect_true(any(grepl("^\\s+x\\s+a\\s+1", out)))
  expect_true(any(grepl("^\\s+x\\s+b\\s+0", out)))
})
