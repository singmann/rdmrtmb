example_data <- data.frame(
  rt = rep(1, 8),
  resp = rep(c("a", "b"), 4),
  cond = factor(rep(c("x", "y"), 4)),
  id = rep(1:4, each = 2),
  item = rep(1:2, 4)
)
design <- rdm_design(example_data, example_data$rt, example_data$resp,
                     formula = list(A = ~1, b = ~1, t0 = ~1, v = ~1))

## a design whose parameter formulas carry random-effect terms
design_re <- function(...) {
  rdm_design(example_data, example_data$rt, example_data$resp,
             formula = rdm_formula(...))
}

test_that("random-effect terms are split off the parameter formulas", {
  d <- design_re(A = ~1 + (1 | id), v = ~1)
  expect_identical(names(d$ranef_terms), "A")
  expect_length(d$ranef_terms$A, 1L)
  ## the bars do not change the fixed part of the model matrix
  expect_identical(colnames(d$A$X), colnames(design$A$X))
  expect_identical(d$A$X, design$A$X)

  r <- rdm_ranef(d)
  expect_s3_class(r, "rdm_ranef")
  expect_identical(length(r$blocks), 1L)
  expect_identical(r$npar, 1L)
  expect_identical(vapply(r$blocks, \(b) b$par, ""), "A")

  ## without bars there are no random effects at all
  expect_identical(length(rdm_ranef(design)$blocks), 0L)
})

test_that("terms are parsed into blocks with the right size", {
  r <- rdm_ranef(design_re(A = ~1 + (1 | id), v = ~1))
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
  r <- rdm_ranef(design_re(A = ~1 + (1 + cond | id), v = ~1))
  b <- r$blocks[[1]]
  expect_identical(b$q, 2L)
  expect_identical(b$n_theta, 3L)
  expect_true(b$correlated)
  expect_equal(b$Z[, 2], rep(as.numeric(example_data$cond == "y"), times = 2))

  ## slopes without an intercept
  r <- rdm_ranef(design_re(A = ~1 + (0 + cond | id), v = ~1))
  expect_identical(r$blocks[[1]]$q, 2L)
})

test_that("crossed random effects become separate blocks", {
  r <- rdm_ranef(design_re(A = ~1 + (1 | id) + (1 | item), v = ~1))
  expect_identical(length(r$blocks), 2L)
  expect_identical(r$npar, 2L)
  expect_identical(vapply(r$blocks, \(b) b$n_lev, 1L), c(4L, 2L))
  expect_match(r$blocks[[2]]$label, "item")
})

test_that("double bars give independent blocks", {
  ## (1 + cond || id) is (1 | id) + (0 + cond | id): the intercept and the
  ## condition effects are uncorrelated, so there is one covariance parameter
  ## more than for (1 + cond | id)
  r <- rdm_ranef(design_re(A = ~1 + (1 + cond || id), v = ~1))
  expect_identical(length(r$blocks), 2L)
  expect_identical(vapply(r$blocks, \(b) b$q, 1L), c(1L, 2L))
  expect_identical(r$npar, 4L)

  correlated <- rdm_ranef(design_re(A = ~1 + (1 + cond | id), v = ~1))
  expect_identical(length(correlated$blocks), 1L)
  expect_identical(correlated$npar, 3L)
})

test_that("diag() gives a diagonal covariance matrix", {
  ## a factor gets one variance per level and no correlations
  r <- rdm_ranef(design_re(A = ~1 + diag(0 + cond | id), v = ~1))
  b <- r$blocks[[1]]
  expect_identical(b$q, 2L)
  expect_identical(b$n_theta, 2L)
  expect_false(b$correlated)
  expect_identical(b$struct, "diag")
  expect_identical(b$label, "A:(diag(0 + cond | id))")

  ## the correlated version needs one parameter more
  expect_identical(rdm_ranef(design_re(A = ~1 + (0 + cond | id), v = ~1))$npar, 3L)

  ## diag() can be combined with other terms
  r <- rdm_ranef(design_re(A = ~1 + (1 | id) + diag(1 + cond | item), v = ~1))
  expect_identical(length(r$blocks), 2L)
  expect_identical(vapply(r$blocks, \(b) b$n_theta, 1L), c(1L, 2L))
  expect_identical(r$npar, 3L)

  ## double bars work inside diag() as well
  r <- rdm_ranef(design_re(A = ~1 + diag(1 + cond || id), v = ~1))
  expect_identical(length(r$blocks), 2L)
  expect_true(all(vapply(r$blocks, \(b) b$struct, "") == "diag"))
})

test_that("the diag structure is part of the fixed part of the formula", {
  ## nobars() drops the whole term, wrapper included
  d <- design_re(A = ~1 + diag(0 + cond | id), v = ~1)
  expect_identical(colnames(d$A$X), "A:(Intercept)")
})

test_that("unknown covariance structures are rejected", {
  expect_error(rdm_ranef(design_re(A = ~1 + us(1 + cond | id), v = ~1)),
               "unsupported covariance structure")
})

test_that("random-effect contributions are the sum of the block terms", {
  r <- rdm_ranef(design_re(A = ~1 + (1 | id) + (0 + cond | item), v = ~1))
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
  ## the diagonal structure has no off-diagonal elements
  expect_equal(rdm_chol(c(log(2), log(0.5)), 2, "diag"),
               diag(c(2, 0.5)))
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

test_that("the diag penalty is a row of independent normal densities", {
  block <- list(q = 2L, n_lev = 4L, struct = "diag")
  theta <- c(log(0.8), log(0.5))
  u <- matrix(c(-1, 2, 0.5, -0.5, 1, -1.5, 0, 1), ncol = 2)
  manual <- -sum(stats::dnorm(u, 0, rep(exp(theta), each = 4), log = TRUE))
  expect_equal(rdm_ranef_penalty(block, u, theta), manual)
})

test_that("random effects are validated", {
  ## the grouping variable has to be in the data
  expect_error(rdm_ranef(design_re(A = ~1 + (1 | nope), v = ~1)), "not found")

  ## a parameter that is held constant cannot have random effects; the terms
  ## are put in by hand because rdm_design() rejects the formula itself
  no_a <- rdm_design(example_data, example_data$rt, example_data$resp,
                     formula = list(b = ~1, t0 = ~1, v = ~1), fixed = list(A = 0.5))
  no_a$ranef_terms <- list(A = reformulas::findbars(~(1 | id)))
  expect_error(rdm_ranef(no_a), "held constant")
})

test_that("the print method summarises the random effects", {
  out <- capture.output(print(rdm_ranef(design_re(A = ~1 + (1 + cond | id),
                                                v = ~1))))
  expect_true(any(grepl("1 term", out)))
  expect_true(any(grepl("correlated", out)))
  expect_true(any(grepl("diagonal", capture.output(
    print(rdm_ranef(design_re(A = ~1 + diag(0 + cond | id), v = ~1)))))))
  expect_true(any(grepl("no random effects",
                        capture.output(print(rdm_ranef(design))))))
})

test_that("drift-rate random effects are evaluated over accumulators", {
  ## a common shift across accumulators
  r <- rdm_ranef(design_re(v = ~1 + (1 | id)))
  b <- r$blocks[[1]]
  expect_identical(b$q, 1L)
  expect_identical(b$n_row, 16L)                    # 8 trials x 2 accumulators
  expect_identical(b$group, rep(rep(1:4, each = 2), 2))

  ## accumulator-specific effects: one column per accumulator
  r <- rdm_ranef(design_re(v = ~1 + (0 + response | id)))
  b <- r$blocks[[1]]
  expect_identical(b$q, 2L)
  expect_identical(b$n_theta, 3L)
  expect_equal(rowSums(b$Z), rep(1, 16))
  expect_equal(b$Z[1:8, 1], rep(1, 8))              # first block: accumulator a
  expect_equal(b$Z[9:16, 2], rep(1, 8))
})
