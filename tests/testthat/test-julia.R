# Tier 2: a real Julia session (PBR_TEST_JULIA=1), and compiling
# (PBR_TEST_JULIAC=1).

test_that("brm fits in Julia and post-processing works", {
  skip_without_julia()
  d <- sim_data(150)
  # what a user waits for: Julia start-up + package load, first fit (JIT),
  # the same fit again, a new model shape
  secs <- function(expr) {
    t0 <- Sys.time()
    force(expr)
    round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
  }
  t_setup <- secs(pbr_setup())
  t_first <- secs(fit <- brm(y ~ x + (1 | g), d, chains = 2, iter = 600, seed = 3, backend = "julia",
                             silent = 2))
  t_again <- secs(brm(y ~ x + (1 | g), d, chains = 2, iter = 600, seed = 4, backend = "julia", silent = 2))
  t_new <- secs(brm(k ~ x, d, family = poisson(), chains = 2, iter = 600, seed = 4, backend = "julia",
                    silent = 2))
  message(sprintf("TIMING pbr_setup %.1f s; first fit %.1f s; same model again %.1f s; new model %.1f s",
                  t_setup, t_first, t_again, t_new))
  expect_s3_class(fit, "pbrfit")
  expect_equal(fit$backend, "julia")
  expect_equal(fit$nobs, 150)
  expect_true(all(c("b_Intercept", "b_x", "sd_g__Intercept", "sigma") %in% variables(fit)))
  expect_equal(dim(fit$draws)[1:2], c(300L, 2L))
  expect_lt(abs(fixef(fit)["x", "Estimate"] - 0.5), 0.25)
  expect_output(print(fit), "Multilevel Hyperparameters")
  expect_equal(dim(posterior_predict(fit)), c(600L, 150L))
  expect_equal(dim(posterior_epred(fit, newdata = d[1:7, ])), c(600L, 7L))
  expect_equal(dim(log_lik(fit)), c(600L, 150L))
  expect_s3_class(loo(fit), "psis_loo")
  expect_equal(dim(bayes_R2(fit)), c(1L, 4L))
  # same seed, same draws
  fit2 <- brm(y ~ x + (1 | g), d, chains = 2, iter = 600, seed = 3, backend = "julia", silent = 2)
  expect_identical(fit2$draws, fit$draws)
  # the fit rebuilt from its draws (as in a later session) predicts the same
  fit$cache$session <- "stale"
  expect_equal(posterior_epred(fit), posterior_epred(fit2))
  gp <- get_prior(y ~ x + (1 | g), d)
  expect_true("sd" %in% gp$class && "Intercept" %in% gp$class)
})

test_that("priors and families reach Julia", {
  skip_without_julia()
  d <- sim_data(100)
  fit <- brm(k ~ x, d, family = poisson(), prior = prior(normal(0, 0.05), class = b), chains = 2,
             iter = 600, seed = 1, backend = "julia", silent = 2)
  expect_lt(abs(fixef(fit)["x", "Estimate"]), 0.2)   # shrunk by the tight prior
  fitd <- brm(bf(y ~ x, sigma ~ x), d, chains = 2, iter = 600, seed = 1, backend = "julia", silent = 2)
  expect_true("b_sigma_x" %in% variables(fitd))
})

test_that("compile once, refit without Julia, same draws as the Julia back end", {
  skip_without_juliac()
  store <- tempfile("store")
  old <- Sys.getenv("PBR_MODEL_STORE")
  Sys.setenv(PBR_MODEL_STORE = store)
  on.exit(Sys.setenv(PBR_MODEL_STORE = old))
  d <- sim_data(150)
  p <- prior(normal(0, 1.5), class = b)
  t0 <- Sys.time()
  m <- compile_model(y ~ x + (1 | g), d, prior = p, verbose = FALSE)
  message("compile_model: ", round(as.numeric(difftime(Sys.time(), t0, units = "secs"))), " s")
  expect_s3_class(m, "pbr_compiled_model")
  expect_equal(names(compiled_models()), m$name)
  fb <- brm(y ~ x + (1 | g), d, prior = p, chains = 2, iter = 600, seed = 5, silent = 2)   # auto -> binary
  expect_equal(fb$backend, "binary")
  fj <- brm(y ~ x + (1 | g), d, prior = p, chains = 2, iter = 600, seed = 5, backend = "julia", silent = 2)
  expect_equal(unclass(fb$draws), unclass(fj$draws), tolerance = 1e-12)
  # new data, no compilation
  d2 <- sim_data(90, 4)
  f2 <- update(fb, newdata = d2)
  expect_equal(f2$backend, "binary")
  expect_equal(f2$nobs, 90L)
  # Julia-side post-processing of a compiled fit
  expect_equal(dim(posterior_predict(f2)), c(600L, 90L))
  # compile_model again, written differently: the same stored model
  expect_equal(compile_model("y ~ x + (1|g)", d2, "gaussian", set_prior("normal(0,1.5)"))$name, m$name)
})
