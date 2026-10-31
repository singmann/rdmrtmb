example_data <- data.frame(
  rt = rep(1, 8),
  resp = rep(c("a", "b"), 4),
  cond = factor(rep(c("x", "y"), 4)),
  id = rep(1:4, each = 2),
  item = rep(1:2, 4)
)
design <- rdm_design(example_data, example_data$rt, example_data$resp,
                     formula = list(A = ~1, b = ~1, t0 = ~1, v = ~1))

test_that("a formula is applied to every parameter", {
  r <- rdm_ranef(~(1 | id), example_data, design)
  expect_s3_class(r, "rdm_ranef")
  ## all free parameters: A, b, t0 and v (s is fixed at 1 by default)
  expect_identical(length(r$blocks), 4L)
  expect_identical(r$npar, 4L)
  expect_identical(vapply(r$blocks, \(b) b$par, ""),
                   c("A", "b", "t0", "v"))
})

test_that("terms are parsed into blocks with the right size", {
  r <- rdm_ranef(list(A = ~(1 | id)), example_data, design)
  b <- r$blocks[[1]]
  expect_identical(b$q, 1L)
  expect_identical(b$n_lev, 4L)
  expect_identical(b$n_theta, 1L)
  expect_false(b$correlated)
  ## the design matrices live on the expanded data, one row per trial and
  ## accumulator, so the random effects are evaluated there too
  expect_identical(b$group, rep(rep(1:4, each = 2), times = 2))
  expect_equal(b$Z, matrix(1, 16, 1))

  ## a correlated intercept and slope: one covariance parameter per element of
  ## the lower triangle
  r <- rdm_ranef(list(A = ~(1 + cond | id)), example_data, design)
  b <- r$blocks[[1]]
  expect_identical(b$q, 2L)
  expect_identical(b$n_theta, 3L)
  expect_true(b$correlated)
  expect_equal(b$Z[, 2], rep(as.numeric(example_data$cond == "y"), times = 2))

  ## slopes without an intercept
  r <- rdm_ranef(list(A = ~(0 + cond | id)), example_data, design)
  expect_identical(r$blocks[[1]]$q, 2L)
})

test_that("crossed random effects become separate blocks", {
  r <- rdm_ranef(list(A = ~(1 | id) + (1 | item)), example_data, design)
  expect_identical(length(r$blocks), 2L)
  expect_identical(r$npar, 2L)
  expect_identical(vapply(r$blocks, \(b) b$n_lev, 1L), c(4L, 2L))
  expect_match(r$blocks[[2]]$label, "item")
})

test_that("random-effect contributions are the sum of the block terms", {
  r <- rdm_ranef(list(A = ~(1 | id) + (0 + cond | item)), example_data, design)
  u1 <- matrix(c(-1, 0, 1, 2), ncol = 1)
  u2 <- matrix(c(0.5, -0.5, 1, -1), ncol = 2)
  c1 <- rdm_ranef_contrib(r$blocks[[1]], u1)
  c2 <- rdm_ranef_contrib(r$blocks[[2]], u2)
  expect_equal(c1, rep(rep(c(-1, 0, 1, 2), each = 2), times = 2))
  ## (0 + cond | item) gives one effect per condition and item
  expect_equal(c2, rep(u2[cbind(example_data$item, as.integer(example_data$cond))],
                       times = 2))
})

test_that("the Cholesky factorisation is correct", {
  ## theta holds the log standard deviations first, then the off-diagonals
  expect_equal(rdm_chol(log(2), 1), matrix(2))
  L <- rdm_chol(c(log(2), log(0.5), 0.3), 2)
  expect_equal(L, matrix(c(2, 0.15, 0, 0.5), 2))
  expect_equal(diag(rdm_chol(c(log(2), log(0.5), 0.3, 0.1, 0.2, -0.3), 3)),
               c(2, 0.5, exp(0.3)))
  expect_true(all(eigen(L %*% t(L), symmetric = TRUE, only.values = TRUE)$values > 0))
})

test_that("the triangular solve is correct", {
  L <- rdm_chol(c(log(2), log(0.5), 0.3), 2)
  U <- matrix(c(1, 2, 3, 4, -1, 0.5, 2, -2), ncol = 2)
  expect_equal(tri_solve(L, U), t(solve(L) %*% t(U)))
  expect_equal(tri_solve(rdm_chol(log(2), 1), U), U / 2)
})

test_that("the random-effect penalty is a multivariate normal density", {
  theta <- c(log(0.8), log(0.6), 0.4)
  u <- matrix(c(-1, 2, 0.5, -0.5, 1, -1.5, 0, 1), ncol = 2)
  block <- list(q = 2L, n_lev = 4L)
  L <- rdm_chol(theta, 2)
  Sigma <- L %*% t(L)
  manual <- 0
  for (i in 1:4) {
    manual <- manual + 0.5 * drop(t(u[i, ]) %*% solve(Sigma) %*% u[i, ]) +
      0.5 * log(det(Sigma)) + log(2 * pi)
  }
  expect_equal(rdm_ranef_penalty(block, u, theta), manual)
})

test_that("random effects are validated", {
  no_a <- rdm_design(example_data, example_data$rt, example_data$resp,
                     formula = list(b = ~1, t0 = ~1, v = ~1), fixed = list(A = 0.5))
  expect_error(rdm_ranef(list(A = ~(1 | id)), example_data, no_a),
               "held constant")
  expect_error(rdm_ranef(list(foo = ~(1 | id)), example_data, design),
               "unknown parameter")
  expect_error(rdm_ranef(~(1 | nope), example_data, design), "not found")
  expect_error(rdm_ranef(~1, example_data, design), "random-effects formula")
})

test_that("the print method summarises the random effects", {
  out <- capture.output(print(rdm_ranef(list(A = ~(1 + cond | id)), example_data,
                                        design)))
  expect_true(any(grepl("1 term", out)))
  expect_true(any(grepl("correlated", out)))
  expect_true(any(grepl("no random effects", capture.output(print(
    rdm_ranef(NULL, example_data, design))))))
})

test_that("drift-rate random effects are evaluated over accumulators", {
  ## a common shift across accumulators
  r <- rdm_ranef(list(v = ~(1 | id)), example_data, design)
  b <- r$blocks[[1]]
  expect_identical(b$q, 1L)
  expect_identical(b$n_row, 16L)                    # 8 trials x 2 accumulators
  expect_identical(b$group, rep(rep(1:4, each = 2), 2))

  ## accumulator-specific effects: one column per accumulator
  r <- rdm_ranef(list(v = ~(0 + response | id)), example_data, design)
  b <- r$blocks[[1]]
  expect_identical(b$q, 2L)
  expect_identical(b$n_theta, 3L)
  expect_equal(rowSums(b$Z), rep(1, 16))
  expect_equal(b$Z[1:8, 1], rep(1, 8))              # first block: accumulator a
  expect_equal(b$Z[9:16, 2], rep(1, 8))
})
