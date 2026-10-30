skip_if_not_installed("rtdists")

simulate_trials <- function(n = 300, seed = 1) {
  set.seed(seed)
  d <- data.frame(rt = 0, response = 0, cond = factor(rep(c("x", "y"), length.out = n)),
                  id = rep(1:3, length.out = n))
  d$rt <- rtdists::rRDM(n, A = 0.5, b = 1, t0 = 0.4, v = c(2.4, 1.6, 1),
                        silent = TRUE)$rt
  set.seed(seed + 1)
  d$response <- rtdists::rRDM(n, A = 0.5, b = 1, t0 = 0.4, v = c(2.4, 1.6, 1),
                              silent = TRUE)$response
  d
}

## implied model parameters for a given coefficient vector
implied <- function(design, beta, u = NULL) {
  par_of <- function(nm) {
    spec <- design[[nm]]
    val <- if (spec$free) drop(spec$X %*% beta[spec$coef_names]) else spec$value
    if (!is.null(u) && !is.null(u[[nm]])) val <- val + u[[nm]]
    val
  }
  A <- exp(par_of("A"))
  list(A = A, b = A + exp(par_of("b")), t0 = exp(par_of("t0")),
       v = matrix(exp(par_of("v")), ncol = design$K))
}

test_that("the objective reproduces the rtdists likelihood without random effects", {
  dat <- simulate_trials()
  d <- rdm_design(dat, dat$rt, dat$response,
                  formula = list(A = ~1, b = ~cond, t0 = ~1,
                                 v = list(~1, ~1, ~1)))
  info <- rdm_nll(d, prior = list(sd = NULL, beta = NULL))
  beta <- c("A:(Intercept)" = log(0.6), "b:(Intercept)" = log(0.7),
            "b:condy" = 0.2, "t0:(Intercept)" = log(0.9),
            "v[1]:(Intercept)" = log(2.5), "v[2]:(Intercept)" = log(1.7),
            "v[3]:(Intercept)" = log(1.1))
  p <- info$par
  p$beta[] <- beta[names(p$beta)]
  p$beta["t0:(Intercept)"] <- log(0.4 * min(dat$rt))   # keep t feasible

  mp <- implied(d, p$beta)
  expected <- -sum(log(rtdists::dRDM(dat$rt, dat$response, A = mp$A, b = mp$b,
                                     t0 = mp$t0,
                                     v = list(mp$v[, 1], mp$v[, 2], mp$v[, 3]),
                                     silent = TRUE)))
  expect_equal(info$fn(p), expected, tolerance = 1e-6)
})

test_that("the gradients of the objective are correct", {
  dat <- simulate_trials(n = 200)
  d <- rdm_design(dat, dat$rt, dat$response,
                  formula = list(A = ~1, b = ~1, t0 = ~1, v = list(~1, ~1, ~1)))
  info <- rdm_nll(d, prior = list(sd = NULL, beta = NULL))
  p <- info$par
  p$beta["t0:(Intercept)"] <- log(0.4 * min(dat$rt))
  obj <- RTMB::MakeADFun(info$fn, p, silent = TRUE)
  expect_equal(as.vector(obj$gr(obj$par)), as.vector(numDeriv::grad(obj$fn, obj$par)),
               tolerance = 1e-4)
})

test_that("the random-effect penalty matches a multivariate normal density", {
  dat <- simulate_trials(n = 300)
  d <- rdm_design(dat, dat$rt, dat$response,
                  formula = list(A = ~1, b = ~1, t0 = ~1, v = list(~1, ~1, ~1)))
  r <- rdm_ranef(list(A = ~(1 | id)), dat, d)
  info <- rdm_nll(d, r, prior = list(sd = NULL, beta = NULL))
  p <- info$par
  p$beta["t0:(Intercept)"] <- log(0.4 * min(dat$rt))
  th <- info$index$theta[[1]]
  p$theta[th] <- log(0.3)

  L <- rdm_chol(p$theta[th], 1)
  u <- p$u[info$index$u[[1]]]
  manual_penalty <- sum(dnorm(u, 0, L[1, 1], log = TRUE) * -1)

  u_contrib <- rdm_ranef_contrib(r$blocks[[1]], u)
  mp <- implied(d, p$beta, list(A = u_contrib))
  expected <- -sum(log(rtdists::dRDM(dat$rt, dat$response, A = mp$A, b = mp$b,
                                     t0 = mp$t0,
                                     v = list(mp$v[, 1], mp$v[, 2], mp$v[, 3]),
                                     silent = TRUE))) + manual_penalty
  expect_equal(info$fn(p), expected, tolerance = 1e-5)
})

test_that("the barrier rejects response times shorter than t0", {
  dat <- simulate_trials(n = 100)
  d <- rdm_design(dat, dat$rt, dat$response,
                  formula = list(A = ~1, b = ~1, t0 = ~1, v = list(~1, ~1, ~1)))
  info <- rdm_nll(d, prior = list(sd = NULL, beta = NULL))
  p <- info$par
  ok <- p
  ok$beta["t0:(Intercept)"] <- log(0.4 * min(dat$rt))
  bad <- p
  bad$beta["t0:(Intercept)"] <- log(2 * max(dat$rt))
  expect_gt(info$fn(bad), info$fn(ok) + 1000)
})

test_that("a fixed diffusion constant rescales the model", {
  dat <- simulate_trials(n = 200)
  d <- rdm_design(dat, dat$rt, dat$response,
                  formula = list(A = ~1, b = ~1, t0 = ~1, v = list(~1, ~1, ~1)))
  info <- rdm_nll(d, prior = list(sd = NULL, beta = NULL))
  p <- info$par
  p$beta["t0:(Intercept)"] <- log(0.4 * min(dat$rt))
  ## s = 1 is the default; a model with all parameters scaled by 2 fits the
  ## same data equally well
  d2 <- rdm_design(dat, dat$rt, dat$response,
                   formula = list(A = ~1, b = ~1, t0 = ~1, v = list(~1, ~1, ~1)),
                   fixed = list(s = 2))
  info2 <- rdm_nll(d2, prior = list(sd = NULL, beta = NULL))
  ## s scales A, b and v, so their coefficients shift by log(2)
  p2 <- info2$par
  scaled <- c("A:(Intercept)", "b:(Intercept)", "v[1]:(Intercept)",
              "v[2]:(Intercept)", "v[3]:(Intercept)")
  p2$beta[scaled] <- p$beta[scaled] + log(2)
  p2$beta["t0:(Intercept)"] <- p$beta["t0:(Intercept)"]
  expect_equal(info2$fn(p2), info$fn(p), tolerance = 1e-6)
})

test_that("the model can be fitted and recovers its parameters", {
  set.seed(5)
  sim <- do.call(rbind, lapply(1:4, \(i) {
    d <- rtdists::rRDM(200, A = 0.5, b = 1, t0 = 0.4, v = c(2.4, 1.6, 1),
                       silent = TRUE)
    data.frame(id = i, rt = d$rt, response = d$response)
  }))
  des <- rdm_design(sim, sim$rt, sim$response,
                    formula = list(A = ~1, b = ~1, t0 = ~1, v = list(~1, ~1, ~1)))
  ran <- rdm_ranef(~(1 | id), sim, des)
  info <- rdm_nll(des, ran, prior = list(sd = 0.5, beta = 2))
  obj <- RTMB::MakeADFun(info$fn, info$par, random = info$random, silent = TRUE)
  ## the inner optimisation of the random effects warns about NaNs while it
  ## searches; that says nothing about the final fit
  fit <- withCallingHandlers(
    stats::nlminb(obj$par, obj$fn, obj$gr, control = list(iter.max = 300)),
    warning = function(w) {
      if (grepl("NA/NaN", conditionMessage(w), fixed = TRUE)) {
        invokeRestart("muffleWarning")
      }
    }
  )
  full <- obj$env$last.par.best
  mp <- info$pars(full[names(full) == "beta"], full[names(full) == "u"])
  expect_lt(fit$convergence, 1)
  ## A is only weakly identified with a handful of participants
  expect_equal(mean(mp$A), 0.5, tolerance = 0.35)
  expect_equal(mean(mp$b), 1, tolerance = 0.3)
  expect_equal(mean(mp$t0), 0.4, tolerance = 0.15)
  expect_equal(colMeans(mp$v), c(2.4, 1.6, 1), tolerance = 0.25)
})
