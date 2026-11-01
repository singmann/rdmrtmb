test_that("fit_rdm returns a fitted object without random effects", {
  dat <- simulate_data(seed = 3)
  f <- fit_rdm(dat, rt = "rt", response = "response",
               formula = list(A = ~1, b = ~1, t0 = ~1, v = ~ acc(1)))
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

test_that("a model without start-point variability can be fitted", {
  dat <- simulate_data(seed = 5)
  f <- fit_rdm(dat, rt = "rt", response = "response",
               formula = list(t0 = ~1, b = ~1, v = ~ acc(1)))
  expect_true(f$converged)
  expect_named(f$beta, c("b:(Intercept)", "t0:(Intercept)", "v[1]:(Intercept)",
                         "v[2]:(Intercept)", "v[3]:(Intercept)"))
  expect_true(is.finite(f$loglik))
  expect_equal(unique(fitted(f)$A), 0)

  ## the likelihood of all trials, from the inverse-Gaussian density and
  ## survival of the fitted parameters
  pars <- fitted(f)
  v <- as.matrix(pars[, c("v1", "v2", "v3")])
  t <- dat$rt - pars$t0
  ref <- vapply(seq_len(nrow(dat)), function(i) {
    j <- dat$response[i]
    exp(rdm_log_f_ig(t[i], pars$b[i], v[i, j])) *
      prod(exp(rdm_log_S_ig(t[i], pars$b[i], v[i, -j])))
  }, 0)
  expect_equal(sum(log(ref)), unname(as.numeric(logLik(f))), tolerance = 1e-8)

  ## the predicted probabilities of the three responses still sum to one
  expect_equal(unname(rowSums(predict(f, type = "probs"))), rep(1, nrow(dat)),
               tolerance = 1e-6)

  ## fixing A at zero is the same model as leaving it out
  fixed_zero <- fit_rdm(dat, rt = "rt", response = "response",
                        formula = list(t0 = ~1, b = ~1, v = ~ acc(1)),
                        fixed = list(A = 0), sdreport = FALSE)
  expect_equal(coef(fixed_zero), coef(f), tolerance = 1e-6)
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
               formula = list(A = ~1, b = ~1, t0 = ~1, v = ~ acc(1)),
               start = c("t0:(Intercept)" = log(0.3)))
  plain <- fit_rdm(dat, rt = "rt", response = "response",
                   formula = list(A = ~1, b = ~1, t0 = ~1, v = ~ acc(1)))
  expect_true(f$converged)
  expect_equal(coef(f), coef(plain), tolerance = 1e-3)
  expect_error(
    fit_rdm(dat, rt = "rt", response = "response",
            formula = list(A = ~1, b = ~1, t0 = ~1, v = ~ acc(1)),
            start = c(nope = 1)),
    "unknown name"
  )
})

test_that("a fit can be restarted from a previous fit", {
  dat <- simulate_data(seed = 6)
  fml <- rdm_formula(A = ~1 + (1 | id), b = ~1, t0 = ~1, v = ~ acc(1))
  first <- fit_rdm(dat, rt = "rt", response = "response", formula = fml)
  ## restarting from a converged fit takes a couple of iterations at most
  again <- fit_rdm(dat, rt = "rt", response = "response", formula = fml,
                   start = first, sdreport = FALSE)
  expect_true(again$converged)
  expect_equal(coef(again), coef(first), tolerance = 1e-4)
  expect_equal(again$theta, first$theta, tolerance = 1e-4)
  expect_lt(again$opt$iterations, first$opt$iterations)
  ## a list with the pieces works too, and wrong lengths are caught
  expect_silent(fit_rdm(dat, rt = "rt", response = "response", formula = fml,
                        start = list(beta = first$beta, theta = first$theta),
                        sdreport = FALSE))
  expect_error(fit_rdm(dat, rt = "rt", response = "response", formula = fml,
                       start = list(theta = c(first$theta, 0)), sdreport = FALSE),
               "values but the model has")
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

test_that("thresholds can differ between accumulators", {
  set.seed(3)
  n <- 900
  pars <- list(A = rep(0.5, n), t0 = rep(0.3, n), s = rep(1, n),
               b = cbind(rep(1.4, n), rep(0.9, n)),
               v = matrix(1.5, n, 2))
  draw <- rdm_rand_race(pars)
  dat <- data.frame(rt = draw$rt, response = draw$response)

  f <- fit_rdm(dat, rt = "rt", response = "response",
               formula = list(A = ~1, t0 = ~1, b = ~ acc(1),
                              v = ~ acc(1)))
  ce <- coef(f, type = "ratio")
  ## the lower threshold wins more often; the coefficients of `b` are the
  ## distance between the threshold and the start point, so they are b - A
  expect_gt(mean(dat$response == 2), 0.5)
  expect_gt(ce[["b[1]:(Intercept)"]], ce[["b[2]:(Intercept)"]])
  expect_equal(unname(ce[["b[1]:(Intercept)"]]), 0.9, tolerance = 0.2)
  expect_equal(unname(ce[["b[2]:(Intercept)"]]), 0.4, tolerance = 0.2)

  ## the fitted() values are one row per trial and accumulator
  fit <- fitted(f)
  expect_equal(nrow(fit), 2 * n)
  expect_true(all(c("response", "A", "b", "v") %in% names(fit)))
  expect_equal(fit$b[fit$response == "1"], rep(1.4, n), tolerance = 0.2)
})

test_that("the extreme eigenvalues of the inner Hessian are reported", {
  f <- small_fit()
  expect_named(f$inner_eig, c("min", "max"))
  expect_gt(f$inner_eig[["min"]], 0)
  expect_lt(f$inner_eig[["min"]], f$inner_eig[["max"]])

  ## printed next to the log-likelihood, and in the summary
  expect_true(any(grepl("inner Hessian", capture.output(print(f)))))
  expect_true(any(grepl("Inner Hessian", capture.output(print(summary(f))))))
  expect_equal(summary(f)$inner_eig, f$inner_eig)

  ## a model without random effects has none
  expect_null(plain_fit()$inner_eig)
  expect_false(any(grepl("inner Hessian", capture.output(print(plain_fit())))))
})

test_that("the extreme eigenvalues of the outer Hessian are reported", {
  f <- small_fit()
  dat <- simulate_data(seed = 2)
  expect_named(f$outer_eig, c("min", "max"))
  expect_gt(f$outer_eig[["min"]], 0)
  expect_lt(f$outer_eig[["min"]], f$outer_eig[["max"]])
  expect_true(any(grepl("outer Hessian", capture.output(print(f)))))
  expect_true(any(grepl("Outer Hessian", capture.output(print(summary(f))))))
  expect_equal(summary(f)$outer_eig, f$outer_eig)

  ## the covariance from sdreport() inverts to the hessian
  f2 <- fit_rdm(dat, rt = "rt", response = "response",
                formula = rdm_formula(A = ~1 + (1 | id), b = ~1, t0 = ~1,
                                      v = ~ acc(1)))
  ev <- 1 / eigen(as.matrix(f2$sd$cov.fixed), symmetric = TRUE,
                  only.values = TRUE)$values
  expect_equal(f2$outer_eig[["min"]] / f2$outer_eig[["max"]],
               min(ev) / max(ev))

  ## without sdreport() there is no outer Hessian
  f3 <- fit_rdm(dat, rt = "rt", response = "response",
                formula = rdm_formula(A = ~1 + (1 | id), b = ~1, t0 = ~1,
                                      v = ~ acc(1)),
                sdreport = FALSE)
  expect_null(f3$outer_eig)
})
