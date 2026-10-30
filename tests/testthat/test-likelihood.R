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

test_that("impossible response times are rejected by the density", {
  ## t <= 0 gives a log density of essentially -Inf rather than NaN
  expect_lt(rdm_log_f(1e-10, A = 0.5, b = 1, v = 2), -100)
})
