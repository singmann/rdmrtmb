test_that("emmeans works on one parameter at a time", {
  skip_if_not_installed("emmeans")
  f <- small_fit()
  expect_error(emmeans::emmeans(f, ~ 1), "dpar")
  expect_error(emmeans::emmeans(f, ~ 1, dpar = "s"), "dpar")
  expect_error(emmeans::emmeans(f, ~ 1, dpar = "nope"), "dpar")
})

test_that("emmeans reports the parameters on the response scale", {
  skip_if_not_installed("emmeans")
  f <- plain_fit()
  ## without random effects the fitted values are the population parameters
  emm <- emmeans::emmeans(f, ~ .acc, dpar = "v")
  est <- as.data.frame(summary(emm, type = "response"))
  expect_equal(est$response, unname(f$pars$v[1, ]), tolerance = 1e-6)
  ## and the link scale is the log of that
  link <- as.data.frame(summary(emm))
  expect_equal(link$emmean, log(est$response), tolerance = 1e-6)

  ## the same for a parameter that is shared by the accumulators
  emm_a <- emmeans::emmeans(f, ~ 1, dpar = "A")
  est_a <- as.data.frame(summary(emm_a, type = "response"))
  expect_equal(est_a$response, f$pars$A[1, 1], tolerance = 1e-6)
  ## the threshold is reported as the distance to the start point
  emm_b <- emmeans::emmeans(f, ~ 1, dpar = "b")
  est_b <- as.data.frame(summary(emm_b, type = "response"))
  expect_equal(est_b$response, f$pars$b[1, 1] - f$pars$A[1, 1], tolerance = 1e-6)
})

test_that("emmeans uses the contrasts of the fit", {
  skip_if_not_installed("emmeans")
  dat <- simulate_data(n_id = 4, n_trials = 50)
  dat$cond <- factor(rep(c("a", "b"), length.out = nrow(dat)))
  old <- options(contrasts = c("contr.sum", "contr.poly"))
  on.exit(options(old), add = TRUE)
  f <- fit_rdm(dat, rt = "rt", response = "response",
               formula = list(A = ~1, b = ~1, t0 = ~1,
                              v = ~ cond + acc(.acc == 2)))
  ## the contrasts in force change after the fit; the reference grid must still
  ## be built with the ones the model was fitted with
  options(contrasts = c("contr.treatment", "contr.poly"))
  emm <- emmeans::emmeans(f, ~ cond + .acc, dpar = "v")
  est <- as.data.frame(summary(emm, type = "response"))
  expect_equal(nrow(est), 6L)
  ## the estimates are exp(X beta) of the cells of the fitted design
  cells <- f$design$v$X
  key <- with(f$design$long, paste(cond, .acc, sep = "|"))
  first <- match(with(est, paste(cond, .acc, sep = "|")), key)
  eta <- as.vector(cells[first, , drop = FALSE] %*%
                     f$beta[f$nll$index$beta$v])
  expect_equal(est$response, exp(eta), tolerance = 1e-6)
})

test_that("emmeans reports the random effects at the population level", {
  skip_if_not_installed("emmeans")
  f <- small_fit()
  ## v is estimated per accumulator; without `.acc` the estimate averages over
  ## the accumulators on the link scale
  emm <- emmeans::emmeans(f, ~ 1, dpar = "v")
  est <- as.data.frame(summary(emm, type = "response"))
  linked <- mean(log(unname(f$pars$v[1, ])))
  expect_equal(log(est$response), linked, tolerance = 1e-6)
})
