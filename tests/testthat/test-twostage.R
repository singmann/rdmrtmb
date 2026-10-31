dat <- simulate_data(n_id = 6, n_trials = 60, seed = 3)
fml <- rdm_formula(A = ~1 + (1 | id), b = ~1, t0 = ~1 + (1 | id),
                   v = ~ acc(1))

test_that("per-subject fits give a usable starting point", {
  st <- rdm_subject_start(dat, rt = "rt", response = "response", formula = fml,
                          id = "id")
  expect_s3_class(st, "rdm_start")
  expect_named(st, c("beta", "theta", "u"))
  fit <- fit_rdm(dat, rt = "rt", response = "response", formula = fml)
  expect_named(st$beta, names(fit$beta))
  expect_length(st$theta, length(fit$theta))
  expect_true(all(is.finite(c(st$beta, st$theta))))
  ## the per-subject estimates are in the region of the group estimates
  expect_equal(st$beta, fit$beta, tolerance = 0.8)
  ## and so are the random-effect standard deviations
  expect_lt(abs(exp(st$theta[1]) - exp(fit$theta[1])), 0.8)
})

test_that("a hierarchical fit starts from those values", {
  st <- rdm_subject_start(dat, rt = "rt", response = "response", formula = fml,
                          id = dat$id)
  warm <- fit_rdm(dat, rt = "rt", response = "response", formula = fml,
                  start = st, sdreport = FALSE)
  cold <- fit_rdm(dat, rt = "rt", response = "response", formula = fml,
                  sdreport = FALSE)
  expect_true(warm$converged)
  expect_equal(coef(warm), coef(cold), tolerance = 1e-3)
  expect_true(is.finite(warm$loglik))
})

test_that("subjects that cannot be fitted are reported, not fatal", {
  ## a subject that only ever chose one alternative cannot be fitted
  bad <- dat
  bad$response[bad$id == 2] <- bad$response[bad$id == 2][1]
  st <- rdm_subject_start(bad, rt = "rt", response = "response", formula = fml,
                          id = "id")
  expect_gt(attr(st, "n_failed"), 0)
  expect_identical(nrow(attr(st, "coef")), 5L)
})

test_that("the starting values print a group summary", {
  st <- rdm_subject_start(dat, rt = "rt", response = "response", formula = fml,
                          id = "id")
  out <- capture.output(print(st))
  expect_match(out[1], "6 subject fits")
  expect_true(any(grepl("A:\\(Intercept\\)", out)))
  expect_error(
    rdm_subject_start(dat, rt = "rt", response = "response", formula = fml,
                      id = rep(1, nrow(dat))),
    "fewer than two subjects"
  )
})
