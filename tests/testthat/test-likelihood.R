# The likelihood in R/rdm-likelihood.R is compared against the implementation
# in rtdists. Densities should agree to about 1e-7; the CDF implied by
# numerically integrating our density is compared against pRDM, which is itself
# computed by numerical integration and therefore only agrees to about 1e-5.

skip_if_not_installed("rtdists")

param_grid <- expand.grid(A = c(0.2, 0.5, 0.9),
                          b = c(0.8, 1.2, 1.8),
                          t0 = c(0.1, 0.4),
                          v1 = c(1.5, 3),
                          v2 = c(1, 2),
                          v3 = c(0.5, 1.5))
param_grid <- param_grid[param_grid$b > param_grid$A, ]

dtimes <- c(0.05, 0.2, 0.5, 1, 2)

rtdists_call <- function(fun, rt, j, g, v) {
  fun(rt, j, A = g$A, b = g$b, t0 = g$t0, v = v, silent = TRUE)
}

test_that("the trial likelihood matches dRDM over a grid of parameters", {
  worst <- 0
  for (i in seq_len(nrow(param_grid))) {
    g <- param_grid[i, ]
    v <- c(g$v1, g$v2, g$v3)
    rt <- g$t0 + dtimes
    vmat <- matrix(rep(v, each = length(dtimes)), ncol = 3)
    for (j in 1:3) {
      mine <- rdm_loglik(rt - g$t0, rep(g$A, length(rt)), rep(g$b, length(rt)),
                         vmat, rep(j, length(rt)))
      rtd <- log(rtdists_call(rtdists::dRDM, rt, j, g, v))
      worst <- max(worst, max(abs(mine - rtd)))
    }
  }
  expect_lt(worst, 1e-6)
})

test_that("the race survival matches 1 - sum(pRDM)", {
  worst <- 0
  for (i in seq_len(nrow(param_grid))) {
    g <- param_grid[i, ]
    v <- c(g$v1, g$v2, g$v3)
    rt <- g$t0 + dtimes
    t <- rt - g$t0
    mine <- exp(Reduce(`+`, lapply(v, \(k) rdm_log_S(t, g$A, g$b, k))))
    rtd <- 1 - rowSums(sapply(1:3, \(j) rtdists_call(rtdists::pRDM, rt, j, g, v)))
    worst <- max(worst, max(abs(mine - rtd)))
  }
  expect_lt(worst, 1e-4)
})

test_that("integrating the density reproduces the race CDF", {
  ## one representative parameter set with well-spread decision times
  g <- list(A = 0.5, b = 1.2, t0 = 0.3, v = c(2.4, 1.6, 1))
  dtimes <- c(0.1, 0.4, 1)

  cdf_mine <- sapply(1:3, \(j) {
    sapply(dtimes, \(x) {
      integrate(\(s) {
        exp(rdm_loglik(s, rep(g$A, length(s)), rep(g$b, length(s)),
                       matrix(rep(g$v, each = length(s)), ncol = 3),
                       rep(j, length(s))))
      }, 0, x, rel.tol = 1e-10)$value
    })
  })
  rt <- g$t0 + dtimes
  cdf_rtd <- sapply(1:3, \(j) rtdists_call(rtdists::pRDM, rt, j, g, g$v))

  expect_lt(max(abs(cdf_mine - cdf_rtd)), 1e-4)
  ## the defective CDFs must sum to the CDF of the race
  expect_equal(rowSums(cdf_mine) + exp(Reduce(`+`, lapply(g$v, \(k) {
    rdm_log_S(dtimes, g$A, g$b, k)
  }))), rep(1, length(dtimes)), tolerance = 1e-6)
})

test_that("the fixed-start likelihood matches the Wald distribution", {
  ## with A = 0 the first-passage time of an accumulator is inverse Gaussian,
  ## which rtdists provides as dwald() and pwald()
  t <- c(0.05, 0.2, 0.5, 1, 2)
  worst_f <- worst_S <- 0
  for (g in list(c(b = 0.8, v = 1.2), c(b = 1.3, v = 2.4), c(b = 2, v = 3))) {
    mine_f <- rdm_log_f_ig(t, g[["b"]], g[["v"]])
    rtd_f <- log(rtdists::dwald(t, A = 0, b = g[["b"]], t0 = 0, v = g[["v"]]))
    worst_f <- max(worst_f, max(abs(mine_f - rtd_f)))
    ## pwald() is the CDF; 1 - it is only accurate away from the upper tail
    cdf <- rtdists::pwald(t, A = 0, b = g[["b"]], t0 = 0, v = g[["v"]])
    ok <- cdf < 1 - 1e-8
    worse <- max(abs(exp(rdm_log_S_ig(t[ok], g[["b"]], g[["v"]])) - (1 - cdf[ok])))
    worst_S <- max(worst_S, worse)
  }
  expect_lt(worst_f, 1e-8)
  expect_lt(worst_S, 1e-8)

  ## the survival is the integral of the density, also far into the tail, where
  ## 1 - pwald() has lost its digits
  g <- c(b = 1.3, v = 2.4)
  for (x in c(0.5, 3, 8)) {
    mass <- integrate(\(s) exp(rdm_log_f_ig(s, g[["b"]], g[["v"]])), 0, x,
                      rel.tol = 1e-10)$value
    expect_equal(mass + exp(rdm_log_S_ig(x, g[["b"]], g[["v"]])), 1,
                 tolerance = 1e-7)
  }
  ## and it decreases in t, so a comparison against dRDM is not needed
  expect_true(all(diff(rdm_log_S_ig(seq(0.1, 5, 0.1), 1.3, 2.4)) < 0))
})

test_that("the fixed-start likelihood is the limit of the start-averaged one", {
  t <- c(0.05, 0.2, 0.5, 1, 2)
  b <- 1.3
  v <- matrix(rep(c(2.4, 1.6, 1), each = length(t)), ncol = 3)
  fixed <- rdm_loglik(t, 0, rep(b, length(t)), v, rep(2L, length(t)),
                      fixed_start = TRUE)
  ## the difference between the two models shrinks in proportion to A
  prev <- Inf
  for (A in c(1e-2, 1e-3, 1e-4)) {
    open <- rdm_loglik(t, rep(A, length(t)), rep(b, length(t)), v, rep(2L, length(t)))
    gap <- max(abs(fixed - open))
    expect_lt(gap, prev)
    prev <- gap
  }
  expect_lt(prev, 1e-2)

  ## a whole trial: density of the winner times the survival of the others
  A0 <- 0
  mine <- exp(rdm_loglik(t, A0, rep(b, length(t)), v, rep(2L, length(t)),
                         fixed_start = TRUE))
  vv <- c(2.4, 1.6, 1)
  ref <- sapply(c(1, 2, 3), \(j) {
    exp(rdm_log_f_ig(t, b, vv[j])) *
      apply(sapply(setdiff(1:3, j), \(k) {
        exp(rdm_log_S_ig(t, b, vv[k]))
      }), 1, prod)
  })
  expect_equal(mine, ref[, 2], tolerance = 1e-10)
})

test_that("impossible response times are rejected by the density", {
  ## t <= 0 gives a log density of essentially -Inf rather than NaN
  expect_lt(rdm_log_f(1e-10, A = 0.5, b = 1, v = 2), -100)
})
