## data and a small fit shared by the tests

simulate_data <- function(n_id = 5, n_trials = 80, seed = 1) {
  set.seed(seed)
  n <- n_id * n_trials
  pars <- list(A = rep(0.6, n), b = rep(1.1, n), t0 = rep(0.35, n), s = rep(1, n),
               v = matrix(rep(c(2.4, 1.6, 1), each = n), n, 3))
  draw <- rdm_rand_race(pars)
  data.frame(id = rep(seq_len(n_id), each = n_trials),
             rt = draw$rt, response = draw$response)
}

small_fit <- local({
  cached <- NULL
  function() {
    if (is.null(cached)) {
      dat <- simulate_data()
      cached <<- rdm_fit(dat, rt = "rt", response = "response",
                         formula = rdm_formula(A = ~1 + (1 | id),
                                               b = ~1, t0 = ~1 + (1 | id),
                                               v = ~ acc(1)))
    }
    cached
  }
})

plain_fit <- local({
  cached <- NULL
  function() {
    if (is.null(cached)) {
      dat <- simulate_data(seed = 2)
      cached <<- rdm_fit(dat, rt = dat$rt, response = dat$response,
                         formula = list(A = ~1, b = ~1, t0 = ~1,
                                        v = ~ acc(1)))
    }
    cached
  }
})
