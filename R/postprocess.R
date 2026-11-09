# Post-processing that needs the model, not just its draws (predictions,
# log-likelihood): done in Julia. A fit made in this session's Julia is used
# directly; any other fit (a compiled model's, or one read from disk) is
# rebuilt in Julia from its unconstrained draws, once per session.

julia_fit <- function(fit) {
  if (!is.null(fit$cache$julia) && identical(fit$cache$session, .pbr$session) && pbr_julia_ready()) {
    return(fit$cache$julia)
  }
  pbr_setup()
  JuliaCall::julia_assign("PBRR_theta", fit$theta)
  JuliaCall::julia_assign("PBRR_stats", fit$stats)
  src <- sprintf(paste0("PBRR.from_draws(%s, %s; family = %s, prior = %s, theta = PBRR_theta, ",
                        "statnames = %s, stats = PBRR_stats, times = %s, iter = %d, warmup = %d)"),
                 jl_formula(fit$formula), jl_data(fit$data, names(fit$data)),
                 jl_family(brmsfamily(fit$family$family, fit$family$link)), jl_prior(fit$prior),
                 jl_strvec(dimnames(fit$stats)$variable),
                 paste0("Float64[", paste(vapply(ifelse(is.na(fit$times), 0, fit$times), jl_num, ""),
                                          collapse = ", "), "]"),
                 fit$args$iter, fit$args$warmup)
  j <- jl_eval(src, "Julia")
  fit$cache$julia <- j
  fit$cache$session <- .pbr$session
  j
}

# call a PracticalBayesianRegression function on a fit, newdata as a keyword
jl_post <- function(fun, fit, newdata = NULL, ...) {
  j <- julia_fit(fit)
  kw <- list(...)
  if (!is.null(newdata)) {
    newdata <- as.data.frame(newdata)
    kw$newdata <- jl_eval(jl_data(newdata, names(newdata)), "Julia")
  }
  do.call(JuliaCall::julia_call, c(list(paste0("PBR.", fun), j), kw, list(need_return = "R")))
}

#' Posterior predictions
#'
#' Draws from the posterior predictive distribution (`posterior_predict`),
#' of its mean (`posterior_epred`), of the linear predictor
#' (`posterior_linpred`) and of the pointwise log-likelihood (`log_lik`), as
#' in 'brms': a matrix with one row per draw and one column per observation.
#' Computed by Julia (started if needed); a fit from a compiled model is
#' rebuilt there from its draws.
#' @param object A `pbrfit`.
#' @param newdata Data to predict for; default the data the model was fitted
#'   to.
#' @param ... Unused.
#' @return A draws x observations matrix.
#' @export
posterior_predict <- function(object, ...) UseMethod("posterior_predict")
#' @rdname posterior_predict
#' @export
posterior_predict.pbrfit <- function(object, newdata = NULL, ...) jl_post("posterior_predict", object, newdata)

#' @rdname posterior_predict
#' @export
posterior_epred <- function(object, ...) UseMethod("posterior_epred")
#' @rdname posterior_predict
#' @export
posterior_epred.pbrfit <- function(object, newdata = NULL, ...) jl_post("posterior_epred", object, newdata)

#' @rdname posterior_predict
#' @export
posterior_linpred <- function(object, ...) UseMethod("posterior_linpred")
#' @rdname posterior_predict
#' @param transform Apply the inverse link (`posterior_linpred`).
#' @export
posterior_linpred.pbrfit <- function(object, transform = FALSE, newdata = NULL, ...) {
  jl_post("posterior_linpred", object, newdata, transform = isTRUE(transform))
}

#' @rdname posterior_predict
#' @export
log_lik <- function(object, ...) UseMethod("log_lik")
#' @rdname posterior_predict
#' @export
log_lik.pbrfit <- function(object, newdata = NULL, ...) jl_post("log_lik", object, newdata)

#' Bayesian R-squared
#' @param object A `pbrfit`.
#' @param probs Quantiles.
#' @param ... Unused.
#' @return A one-row matrix, as 'brms''s `bayes_R2`.
#' @export
bayes_R2 <- function(object, ...) UseMethod("bayes_R2")
#' @rdname bayes_R2
#' @export
bayes_R2.pbrfit <- function(object, probs = c(0.025, 0.975), ...) {
  e <- posterior_epred.pbrfit(object)
  y <- JuliaCall::julia_call("PBRR.response", julia_fit(object), need_return = "R")
  r2 <- apply(e, 1, function(m) {
    ve <- stats::var(m)
    ve / (ve + stats::var(y - m))
  })
  out <- matrix(c(mean(r2), stats::sd(r2), stats::quantile(r2, probs, names = FALSE)), nrow = 1,
                dimnames = list("R2", c("Estimate", "Est.Error", paste0("Q", probs * 100))))
  out
}

#' Leave-one-out cross-validation
#'
#' PSIS-LOO and WAIC through the 'loo' package, from the pointwise
#' log-likelihood: the results work with `loo::loo_compare()` as 'brms''s do.
#' @param x A `pbrfit`.
#' @param ... Passed to [loo::loo()] / [loo::waic()].
#' @return A `psis_loo` / `waic` object.
#' @name loo.pbrfit
NULL

#' @importFrom loo loo waic loo_compare
#' @export
loo::loo
#' @export
loo::waic
#' @export
loo::loo_compare

#' @rdname loo.pbrfit
#' @exportS3Method loo::loo
loo.pbrfit <- function(x, ...) {
  ll <- log_lik.pbrfit(x)
  chain <- rep(seq_len(nchains.pbrfit(x)), each = niterations.pbrfit(x))
  r_eff <- loo::relative_eff(exp(ll - apply(ll, 2, max)[col(ll)]), chain_id = chain)
  loo::loo(ll, r_eff = r_eff, ...)
}

#' @rdname loo.pbrfit
#' @exportS3Method loo::waic
waic.pbrfit <- function(x, ...) loo::waic(log_lik.pbrfit(x), ...)

#' Default priors of a model
#'
#' The priors 'brms' would use, one row per parameter class, coefficient and
#' group (as 'brms''s `get_prior`); an empty prior is flat. Needs Julia.
#' @inheritParams brm
#' @return A `pbr_prior` data frame with a `source` column.
#' @export
get_prior <- function(formula, data, family = gaussian()) {
  formula <- as_pbr_formula(formula)
  src <- sprintf("PBRR.prior_rows(PBR.get_prior(%s, %s; family = %s))", jl_formula(formula),
                 jl_data(as.data.frame(data), formula_vars(formula)), jl_family(family))
  rows <- jl_eval(src)
  parts <- strsplit(as.character(rows), "\t", fixed = TRUE)
  p <- lapply(parts, function(r) c(r, rep("", 9 - length(r))))
  out <- set_prior(vapply(p, `[`, "", 1), class = vapply(p, `[`, "", 2), coef = vapply(p, `[`, "", 3),
                   group = vapply(p, `[`, "", 4), resp = vapply(p, `[`, "", 5), dpar = vapply(p, `[`, "", 6),
                   lb = ifelse(nzchar(vapply(p, `[`, "", 7)), vapply(p, `[`, "", 7), NA),
                   ub = ifelse(nzchar(vapply(p, `[`, "", 8)), vapply(p, `[`, "", 8), NA))
  out$source <- vapply(p, `[`, "", 9)
  out
}
