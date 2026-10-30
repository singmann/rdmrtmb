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
  expect_identical(dim(d$A$X), c(6L, 1L))
  expect_identical(dim(d$b$X), c(6L, 2L))
  expect_identical(d$b$coef_names, c("b:(Intercept)", "b:condy"))
  ## a single drift formula shares one coefficient across accumulators
  expect_identical(d$v$coef_names, "v:(Intercept)")
  ## s is fixed at 1 by default
  expect_false(d$s$free)
  expect_identical(unique(d$s$value), 1)
})

test_that("drift rates given as a list are stacked accumulator-major", {
  d <- rdm_design(example_data, example_data$rt, example_data$resp,
                  formula = list(A = ~1, t0 = ~1,
                                 v = list(~1, ~cond)))
  X <- d$v$X
  expect_identical(dim(X), c(12L, 3L))
  ## rows 1:6 belong to accumulator a, rows 7:12 to accumulator b
  expect_equal(X[1:6, 1], rep(1, 6))
  expect_equal(X[7:12, 2], rep(1, 6))
  expect_equal(unname(X[7:12, 1:2]), cbind(rep(0, 6), rep(1, 6)))
  expect_equal(X[7:12, 3], c(0, 0, 0, 1, 1, 1))
})

test_that("a single drift formula describes the accumulator via `response`", {
  d <- rdm_design(example_data, example_data$rt, example_data$resp,
                  formula = list(A = ~1, t0 = ~1, v = ~ 0 + response))
  expect_identical(colnames(d$v$X), c("responsea", "responseb"))
  expect_identical(d$v$coef_names, c("v:responsea", "v:responseb"))
  ## exactly one accumulator per row is active
  expect_true(all(rowSums(d$v$X) == 1))
  expect_equal(d$v$X[1:6, 1], rep(1, 6))
  expect_equal(d$v$X[7:12, 2], rep(1, 6))
})

test_that("the two drift interfaces describe the same linear predictor", {
  d1 <- rdm_design(example_data, example_data$rt, example_data$resp,
                   formula = list(A = ~1, t0 = ~1, v = list(~1, ~1)))
  d2 <- rdm_design(example_data, example_data$rt, example_data$resp,
                   formula = list(A = ~1, t0 = ~1, v = ~ 0 + response))
  expect_equal(as.numeric(d1$v$X %*% c(2, 3)),
               as.numeric(d2$v$X %*% c(2, 3)))
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
               "list of 2")
})

test_that("the print method summarises the design", {
  d <- rdm_design(example_data, example_data$rt, example_data$resp,
                  formula = list(A = ~1, b = ~cond, t0 = ~1, v = ~1))
  out <- capture.output(print(d))
  expect_true(any(grepl("6 trials, 2 accumulators", out)))
  expect_true(any(grepl("b: 2 coefficient", out)))
})

test_that("every parameter can be given one formula per accumulator", {
  d <- rdm_design(example_data, example_data$rt, example_data$resp,
                  formula = list(A = ~1, t0 = ~1, b = list(~1, ~cond),
                                 v = list(~1, ~1)))
  expect_true(d$b$per_accumulator)
  expect_identical(d$b$coef_names,
                   c("b[a]:(Intercept)", "b[b]:(Intercept)", "b[b]:condy"))
  expect_identical(dim(d$b$X), c(12L, 3L))
  ## rows 1:6 belong to accumulator a, rows 7:12 to accumulator b
  expect_equal(d$b$X[1:6, 1], rep(1, 6))
  expect_equal(d$b$X[7:12, 2], rep(1, 6))
  expect_equal(d$b$X[7:12, 3], c(0, 0, 0, 1, 1, 1))
  expect_equal(d$b$X[1:6, 3], rep(0, 6))
  expect_false(d$A$per_accumulator)
  expect_false(d$s$free)   # the default fixed s
})

test_that("a fixed parameter can take one value per accumulator", {
  d <- rdm_design(example_data, example_data$rt, example_data$resp,
                  formula = list(A = ~1, t0 = ~1, v = ~1),
                  fixed = list(s = 1, b = c(1.2, 0.8)))
  expect_true(d$b$expanded)
  expect_equal(d$b$value, rep(c(1.2, 0.8), each = 6))
})

test_that("matchfun adds the latent match factor to the expanded data", {
  dat <- data.frame(rt = runif(6, 0.5, 1.5), resp = rep(c("a", "b"), 3),
                    S = rep(c("a", "b"), 3))
  d <- rdm_design(dat, dat$rt, dat$resp, formula = list(A = ~1, t0 = ~1, v = ~ S),
                  matchfun = function(x) x$S == x$lR)
  expect_identical(levels(d$long$lM), c("FALSE", "TRUE"))
  expect_identical(as.character(d$long$lM[1:6]), as.character(dat$S == "a"))
  expect_identical(as.character(d$long$lM[7:12]), as.character(dat$S == "b"))
  expect_identical(levels(d$long$lR), c("a", "b"))
  expect_error(rdm_design(dat, dat$rt, dat$resp, formula = list(A = ~1, t0 = ~1, v = ~1),
                          matchfun = function(x) TRUE),
               "one value per trial and")
})
