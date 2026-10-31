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

test_that("homdiag shares one variance and cs adds one correlation", {
  r <- rdm_ranef(design_re(A = ~1 + homdiag(0 + cond | id), v = ~1))
  b <- r$blocks[[1]]
  expect_identical(b$q, 2L)
  expect_identical(b$n_theta, 1L)
  expect_false(b$correlated)
  expect_identical(b$label, "A:(homdiag(0 + cond | id))")

  r <- rdm_ranef(design_re(A = ~1 + cs(0 + cond | id), v = ~1))
  expect_identical(r$blocks[[1]]$n_theta, 2L)
  expect_identical(r$blocks[[1]]$struct, "cs")

  ## a single-column term has no correlation to estimate
  expect_identical(rdm_ranef(design_re(A = ~1 + cs(1 | id), v = ~1))$npar, 1L)

  ## the wrappers work for several terms and group factors
  r <- rdm_ranef(design_re(A = ~1 + homdiag(1 | id) + cs(1 + cond | item),
                           v = ~1))
  expect_identical(vapply(r$blocks, \(b) b$n_theta, 1L), c(1L, 2L))
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

test_that("acc() can be used on the left of the bar", {
  ## acc(Condition) in the bar gives the same columns as acc(Condition) in the
  ## fixed part: one intercept and one contrast per accumulator
  d <- design_re(A = ~1, v = ~ acc(cond) + (acc(cond) | id))
  b <- rdm_ranef(d)$blocks[[1]]
  expect_identical(b$q, 4L)
  expect_identical(b$n_theta, 10L)
  expect_identical(b$col_names, colnames(d$v$X))
  ## row j of Z is the rotated column of the fixed design of that row
  expect_equal(b$Z[, 1], as.numeric(d$v$X[, 1]))

  ## the covariance wrappers apply to the larger block
  expect_identical(rdm_ranef(design_re(A = ~1, v = ~ acc(cond) +
                                         diag(acc(cond) | id)))$blocks[[1]]$n_theta,
                   4L)

  ## terms that compare the accumulator with a trial variable are available
  ## here as well, with as many columns as in the fixed part
  d <- design_re(A = ~1, v = ~ acc(.acc == "b") + (acc(.acc == "b") | id))
  b <- rdm_ranef(d)$blocks[[1]]
  expect_identical(b$q, ncol(d$v$X))
  expect_identical(b$col_names, colnames(d$v$X))

  ## the accumulator is available through `response` as well
  expect_identical(rdm_ranef(design_re(A = ~1, v = ~1 + (0 + response | id)))$blocks[[1]]$q,
                   2L)

  ## a bar that repeats a shared term is rank deficient
  expect_error(rdm_ranef(design_re(A = ~1, v = ~ acc(cond) + (cond + acc(cond) | id))),
               "rank deficient")
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

test_that("the compound-symmetry correlation covers the admissible range", {
  ## theta = 0 is no correlation, whatever the number of columns
  expect_equal(vapply(c(2, 3, 5), \(q) cs_rho(0, q), 0), rep(0, 3))
  ## the bounds are the singular matrices
  expect_equal(cs_rho(-Inf, 3), -1/2)
  expect_equal(cs_rho(Inf, 3), 1)
  expect_equal(cs_rho(c(-1, 0.5, 2), 3), cs_rho(cs_theta(cs_rho(c(-1, 0.5, 2), 3), 3), 3))
})

test_that("the compound-symmetry factor reproduces the covariance matrix", {
  sigma <- 1.5
  rho <- 0.4
  q <- 4
  th <- c(log(sigma), cs_theta(rho, q))
  L <- rdm_chol(th, q, "cs")
  Sigma <- sigma^2 * ((1 - rho) * diag(q) + rho * matrix(1, q, q))
  expect_equal(L %*% t(L), Sigma)
  ## the log determinant shortcut is the sum of the logs of the diagonal
  expect_equal(rdm_chol_logdet(th, q, "cs"), sum(log(diag(L))))
  expect_equal(2 * rdm_chol_logdet(th, q, "cs"), log(det(Sigma)))

  ## a negative correlation is allowed as long as the matrix stays positive
  ## definite
  th <- c(log(sigma), cs_theta(-0.2, q))
  L <- rdm_chol(th, q, "cs")
  Sigma <- sigma^2 * ((1 + 0.2) * diag(q) - 0.2 * matrix(1, q, q))
  expect_equal(L %*% t(L), Sigma)

  ## the homogeneous diagonal structure
  expect_equal(rdm_chol(log(2), 3, "homdiag"), diag(2, 3))
  expect_equal(rdm_chol_logdet(log(2), 3, "homdiag"), 3 * log(2))
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

test_that("the homdiag and cs penalties are multivariate normal densities", {
  q <- 3
  u <- matrix(c(-1, 2, 0.5, -0.5, 1, -1.5, 0, 1, 2), ncol = q)
  manual_penalty <- function(theta, block, Sigma) {
    sum(0.5 * rowSums((u %*% solve(Sigma)) * u) + 0.5 * log(det(Sigma)) +
          q / 2 * log(2 * pi))
  }

  ## homogeneous diagonal: one variance for all three columns
  theta <- log(0.7)
  block <- list(q = q, n_lev = 3L, struct = "homdiag")
  expect_equal(rdm_ranef_penalty(block, u, theta),
               manual_penalty(theta, block, 0.7^2 * diag(q)))

  ## compound symmetry: common variance and common correlation
  sigma <- 1.2
  rho <- 0.5
  theta <- c(log(sigma), cs_theta(rho, q))
  block <- list(q = q, n_lev = 3L, struct = "cs")
  Sigma <- sigma^2 * ((1 - rho) * diag(q) + rho * matrix(1, q, q))
  expect_equal(rdm_ranef_penalty(block, u, theta), manual_penalty(theta, block, Sigma))
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
  out <- capture.output(print(rdm_ranef(design_re(A = ~1 + diag(0 + cond | id),
                                                 v = ~1))))
  expect_true(any(grepl(", diag", out)))
  out <- capture.output(print(rdm_ranef(design_re(A = ~1 + cs(0 + cond | id),
                                                 v = ~1))))
  expect_true(any(grepl(", cs", out)))
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
