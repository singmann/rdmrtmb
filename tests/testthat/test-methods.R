test_that("coef() reports the coefficients on either scale", {
  f <- plain_fit()
  link <- coef(f)
  expect_equal(link, f$beta)
  expect_equal(coef(f, type = "ratio"), exp(f$beta))
  full <- coef(f, full = TRUE)
  expect_length(full, length(f$beta) + length(f$theta))
  expect_equal(fixef(f), link)
})

test_that("vcov() and confint() agree with the standard errors", {
  f <- small_fit()
  V <- vcov(f)
  expect_equal(dim(V), c(6L, 6L))
  expect_equal(V, t(V))
  expect_equal(rownames(V), names(f$beta))
  expect_equal(dim(vcov(f, full = TRUE)), c(8L, 8L))

  ci <- confint(f)
  expect_equal(dim(ci), c(6L, 2L))
  expect_equal(rownames(ci), names(f$beta))
  expect_true(all(ci[, 1] < f$beta & f$beta < ci[, 2]))
  ## the interval width is symmetric on the link scale
  z <- stats::qnorm(0.975)
  expect_equal(ci[, 2] - ci[, 1], 2 * z * sqrt(diag(V)),
               ignore_attr = TRUE)

  ci90 <- confint(f, level = 0.9, parm = "A:(Intercept)")
  expect_equal(nrow(ci90), 1L)
  expect_equal(ci90[, 1], f$beta["A:(Intercept)"] -
                 stats::qnorm(0.95) * sqrt(V["A:(Intercept)", "A:(Intercept)"]),
               ignore_attr = TRUE)

  ratio <- confint(f, type = "ratio")
  expect_equal(ratio, exp(ci), ignore_attr = TRUE)
  expect_error(confint(f, parm = "nope"), "unknown coefficient")
})

test_that("logLik(), AIC(), BIC() and nobs() are consistent", {
  f <- small_fit()
  ll <- logLik(f)
  expect_s3_class(ll, "logLik")
  expect_equal(as.numeric(ll), f$loglik)
  expect_equal(attr(ll, "df"), f$df)
  expect_equal(attr(ll, "nobs"), f$nobs)
  expect_equal(AIC(f), -2 * f$loglik + 2 * f$df)
  expect_equal(BIC(f), -2 * f$loglik + log(f$nobs) * f$df)
  expect_equal(nobs(f), nrow(f$design$data))
})

test_that("ranef() returns the conditional modes in long format", {
  f <- small_fit()
  re <- ranef(f)
  expect_s3_class(re, "data.frame")
  expect_named(re, c("term", "group", "effect", "estimate"))
  expect_equal(nrow(re), 2L * 5L)
  expect_setequal(unique(re$term), c("A:(1 | id)", "t0:(1 | id)"))
  expect_equal(ranef(f, type = "ratio")$estimate, exp(re$estimate))
  expect_null(ranef(plain_fit()))
})

test_that("fitted() and predict() return trial-level parameters", {
  f <- small_fit()
  fit <- fitted(f)
  expect_equal(nrow(fit), f$nobs)
  expect_named(fit, c("A", "b", "t0", "s", "v1", "v2", "v3"))
  expect_equal(fit$A, f$pars$A[, 1])
  expect_true(all(fit$b > fit$A))
  expect_equal(predict(f), fit)
  expect_equal(predict(f, type = "parameters"), fit)
})

test_that("predict() gives choice probabilities that sum to one", {
  f <- small_fit()
  p <- predict(f, type = "probs")
  expect_equal(dim(p), c(f$nobs, 3L))
  expect_equal(unname(rowSums(p)), rep(1, f$nobs))
  ## accumulators with a larger drift rate win more often
  expect_true(all(p[, 1] > p[, 2] & p[, 2] > p[, 3]))
})

test_that("predict() gives response-time quantiles that increase", {
  f <- small_fit()
  q <- predict(f, type = "quantiles", grid = 200)
  expect_equal(dim(q), c(f$nobs, 3L, 5L))
  expect_true(all(apply(q, c(1, 2), \(x) all(diff(x) > 0))))
  expect_true(all(q > fitted(f)$t0))
})

test_that("simulate() draws data with the same covariates", {
  f <- plain_fit()
  s <- simulate(f, seed = 11)
  expect_equal(nrow(s), f$nobs)
  expect_true(all(c("rt", "response") %in% names(s)))
  expect_true(all(s$rt > 0))
  expect_setequal(levels(s$response), f$design$response_levels)
  expect_equal(simulate(f, seed = 11), s)

  s2 <- simulate(f, nsim = 3, seed = 12)
  expect_equal(nrow(s2), 3 * f$nobs)
  expect_equal(sort(unique(s2$.sim)), 1:3)
})

test_that("the simulated data match the predicted distribution", {
  f <- plain_fit()
  n <- 20000
  mp <- list(A = rep(f$pars$A[1], n), b = rep(f$pars$b[1], n),
             t0 = rep(f$pars$t0[1], n), s = rep(1, n),
             v = matrix(rep(f$pars$v[1, ], each = n), n, 3))
  set.seed(21)
  draw <- rdm_rand_race(mp)
  p <- predict(f, type = "probs")[1, ]
  observed <- as.numeric(table(draw$response) / n)
  expect_equal(observed, unname(p), tolerance = 0.02)

  q <- predict(f, type = "quantiles")[1, 1, ]
  expect_equal(unname(quantile(draw$rt, c(0.1, 0.5, 0.9))),
               unname(q[c(1, 3, 5)]), tolerance = 0.02)
})

test_that("print() and summary() work", {
  f <- small_fit()
  expect_output(print(f), "racing diffusion model")
  s <- summary(f)
  expect_s3_class(s, "summary.rdm_fit")
  expect_named(s$coefficients, NULL)
  expect_equal(nrow(s$coefficients), length(f$beta))
  expect_equal(s$AIC, AIC(f))
  ## the effects of a block carry the name of the parameter and of the term, so
  ## that they can be matched to the fixed effects
  expect_equal(s$variance$estimate[s$variance$parameter == "sd(A:(Intercept))"],
               exp(f$theta[grep("^A:", names(f$theta))]),
               ignore_attr = TRUE)
  expect_output(print(s), "Fixed effects")
  expect_output(print(s), "Random effects")
})
