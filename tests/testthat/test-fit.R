test_that("fit_rdm returns a fitted object without random effects", {
  dat <- simulate_data(seed = 3)
  f <- fit_rdm(dat, rt = "rt", response = "response",
               formula = list(A = ~1, b = ~1, t0 = ~1, v = list(~1, ~1, ~1)))
  expect_s3_class(f, "rdm_fit")
  expect_true(f$converged)
  expect_named(f$beta, c("A:(Intercept)", "b:(Intercept)", "t0:(Intercept)",
                         "v[1]:(Intercept)", "v[2]:(Intercept)",
                         "v[3]:(Intercept)"))
  expect_length(f$theta, 0)
  expect_equal(f$df, 6L)
  expect_equal(f$nobs, nrow(dat))
  expect_length(f$ranef$blocks, 0)

  ## the reported likelihood is the objective without the priors
  flat <- rdm_nll(f$design, f$ranef, prior = list(sd = NULL, beta = NULL))
  expect_equal(f$loglik, -flat$fn(list(beta = f$beta, u = f$u)))
  expect_lt(f$max_grad, 1e-2)
})

test_that("fit_rdm estimates random effects", {
  f <- small_fit()
  expect_true(f$converged)
  expect_length(f$theta, 2)
  expect_named(f$theta, c("A:(1 | id).1", "t0:(1 | id).1"))
  expect_equal(f$df, 6L + 2L)
  expect_length(f$u, 10L)
  expect_true(is.logical(f$pdHess) || is.na(f$pdHess))
})

test_that("starting values can be supplied and are checked", {
  dat <- simulate_data(seed = 4)
  f <- fit_rdm(dat, rt = "rt", response = "response",
               formula = list(A = ~1, b = ~1, t0 = ~1, v = list(~1, ~1, ~1)),
               start = c("t0:(Intercept)" = log(0.3)))
  plain <- fit_rdm(dat, rt = "rt", response = "response",
                   formula = list(A = ~1, b = ~1, t0 = ~1, v = list(~1, ~1, ~1)))
  expect_true(f$converged)
  expect_equal(coef(f), coef(plain), tolerance = 1e-3)
  expect_error(
    fit_rdm(dat, rt = "rt", response = "response",
            formula = list(A = ~1, b = ~1, t0 = ~1, v = list(~1, ~1, ~1)),
            start = c(nope = 1)),
    "unknown name"
  )
})

test_that("a model with covariates on the drift rates can be fitted", {
  dat <- simulate_data(seed = 5)
  dat$cond <- factor(rep(c("a", "b", "c"), length.out = nrow(dat)))
  f <- fit_rdm(dat, rt = "rt", response = "response",
               formula = list(A = ~1, b = ~1, t0 = ~1, v = ~cond))
  expect_true(f$converged)
  expect_true(any(grepl("cond", names(f$beta))))
  expect_equal(nrow(fitted(f)), nrow(dat))
})
