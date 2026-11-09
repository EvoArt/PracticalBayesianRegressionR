# The fit object: draws held in R (posterior's draws_array), plus what is
# needed to rebuild the fit in Julia for predictions (`theta`, the
# unconstrained draws), so a saved fit works in a later session.

#' @importFrom posterior as_draws as_draws_array as_draws_df as_draws_matrix as_draws_list
#' @importFrom posterior variables ndraws nchains niterations
#' @export
posterior::as_draws
#' @export
posterior::as_draws_array
#' @export
posterior::as_draws_df
#' @export
posterior::as_draws_matrix
#' @export
posterior::as_draws_list
#' @export
posterior::variables
#' @export
posterior::ndraws
#' @export
posterior::nchains
#' @export
posterior::niterations

#' @importFrom stats nobs update
NULL

#' Draws of a fit
#'
#' A fit converts to any of 'posterior''s draws formats, which then work with
#' 'posterior', 'bayesplot' and friends. Sampler statistics (`lp__`,
#' `divergent__`, ...) are in [nuts_params()].
#' @param x A `pbrfit`.
#' @param variable Variables to keep (names or a regular expression with
#'   `regex = TRUE`).
#' @param regex Treat `variable` as regular expressions.
#' @param ... Passed on to 'posterior'.
#' @name as_draws.pbrfit
#' @return A draws object.
NULL

#' @rdname as_draws.pbrfit
#' @exportS3Method posterior::as_draws_array
as_draws_array.pbrfit <- function(x, variable = NULL, regex = FALSE, ...) {
  d <- posterior::as_draws_array(x$draws)
  if (!is.null(variable)) d <- posterior::subset_draws(d, variable = variable, regex = regex)
  d
}
#' @rdname as_draws.pbrfit
#' @exportS3Method posterior::as_draws
as_draws.pbrfit <- function(x, variable = NULL, regex = FALSE, ...) as_draws_array.pbrfit(x, variable, regex)
#' @rdname as_draws.pbrfit
#' @exportS3Method posterior::as_draws_df
as_draws_df.pbrfit <- function(x, variable = NULL, regex = FALSE, ...) {
  posterior::as_draws_df(as_draws_array.pbrfit(x, variable, regex))
}
#' @rdname as_draws.pbrfit
#' @exportS3Method posterior::as_draws_matrix
as_draws_matrix.pbrfit <- function(x, variable = NULL, regex = FALSE, ...) {
  posterior::as_draws_matrix(as_draws_array.pbrfit(x, variable, regex))
}
#' @rdname as_draws.pbrfit
#' @exportS3Method posterior::as_draws_list
as_draws_list.pbrfit <- function(x, variable = NULL, regex = FALSE, ...) {
  posterior::as_draws_list(as_draws_array.pbrfit(x, variable, regex))
}

#' @exportS3Method posterior::variables
variables.pbrfit <- function(x, ...) dimnames(x$draws)$variable
#' @exportS3Method posterior::ndraws
ndraws.pbrfit <- function(x) dim(x$draws)[1] * dim(x$draws)[2]
#' @exportS3Method posterior::nchains
nchains.pbrfit <- function(x) dim(x$draws)[2]
#' @exportS3Method posterior::niterations
niterations.pbrfit <- function(x) dim(x$draws)[1]
#' @export
nobs.pbrfit <- function(object, ...) object$nobs

#' @export
as.data.frame.pbrfit <- function(x, row.names = NULL, optional = FALSE, ...) {
  as.data.frame(as_draws_matrix.pbrfit(x))
}
#' @export
as.matrix.pbrfit <- function(x, ...) unclass(as_draws_matrix.pbrfit(x))
#' @export
as.array.pbrfit <- function(x, ...) unclass(x$draws)

#' Sampler diagnostics
#'
#' Per-draw NUTS statistics, as Stan reports them: `accept_stat__`,
#' `stepsize__`, `treedepth__`, `n_leapfrog__`, `divergent__`, `energy__`,
#' `lp__`, in long format (as 'brms''s `nuts_params()`).
#' @param object A `pbrfit`.
#' @param ... Unused.
#' @return A data frame with columns Chain, Iteration, Parameter, Value.
#' @export
nuts_params <- function(object, ...) UseMethod("nuts_params")
#' @export
nuts_params.pbrfit <- function(object, ...) {
  s <- object$stats
  d <- dim(s)
  data.frame(Chain = rep(rep(seq_len(d[2]), each = d[1]), d[3]),
             Iteration = rep(seq_len(d[1]), d[2] * d[3]),
             Parameter = factor(rep(dimnames(s)$variable, each = d[1] * d[2])),
             Value = as.vector(s))
}

# --------------------------------------------------------------- summaries

summ_rows <- function(fit, vars, prob = 0.95) {
  a <- (1 - prob) / 2
  rows <- lapply(vars, function(v) {
    x <- fit$draws[, , v]
    if (is.null(dim(x))) x <- matrix(x, ncol = 1)
    c(Estimate = mean(x), Est.Error = stats::sd(x),
      stats::setNames(stats::quantile(x, c(a, 1 - a), names = FALSE),
                      paste0(c("l-", "u-"), round(prob * 100), "% CI")),
      Rhat = posterior::rhat(x), Bulk_ESS = posterior::ess_bulk(x), Tail_ESS = posterior::ess_tail(x))
  })
  out <- do.call(rbind, rows)
  if (is.null(out)) out <- matrix(numeric(0), 0, 7)
  rownames(out) <- vars
  out
}

#' Posterior summary
#'
#' Mean, standard deviation and quantiles of each parameter (as 'brms''s
#' `posterior_summary`).
#' @param x A `pbrfit`.
#' @param variable Variables to summarise (default: all).
#' @param probs Quantiles.
#' @param ... Unused.
#' @return A matrix.
#' @export
posterior_summary <- function(x, ...) UseMethod("posterior_summary")
#' @rdname posterior_summary
#' @export
posterior_summary.pbrfit <- function(x, variable = NULL, probs = c(0.025, 0.975), ...) {
  vars <- variable %||% variables.pbrfit(x)
  out <- t(vapply(vars, function(v) {
    s <- as.vector(x$draws[, , v])
    c(Estimate = mean(s), Est.Error = stats::sd(s), stats::quantile(s, probs, names = FALSE))
  }, numeric(2 + length(probs))))
  colnames(out) <- c("Estimate", "Est.Error", paste0("Q", probs * 100))
  out
}

#' Population-level effects
#' @param object A `pbrfit`.
#' @param probs Quantiles.
#' @param ... Unused.
#' @return A matrix, one row per coefficient (as 'brms''s `fixef`).
#' @export
fixef <- function(object, ...) UseMethod("fixef")
#' @rdname fixef
#' @export
fixef.pbrfit <- function(object, probs = c(0.025, 0.975), ...) {
  v <- grep("^b_", variables.pbrfit(object), value = TRUE)
  out <- posterior_summary.pbrfit(object, v, probs)
  rownames(out) <- sub("^b_", "", v)
  out
}

#' Group-level effects
#'
#' Per-level effects of each grouping factor, as 'brms''s `ranef`: a list
#' with one array (levels x statistics x coefficients) per grouping factor.
#' @param object A `pbrfit`.
#' @param probs Quantiles.
#' @param ... Unused.
#' @return A named list of arrays.
#' @export
ranef <- function(object, ...) UseMethod("ranef")
#' @rdname ranef
#' @export
ranef.pbrfit <- function(object, probs = c(0.025, 0.975), ...) {
  v <- grep("^r_", variables.pbrfit(object), value = TRUE)
  if (!length(v)) return(list())
  grp <- sub("^r_(.*)\\[.*$", "\\1", v)
  lev <- sub("^.*\\[(.*),[^,]*\\]$", "\\1", v)
  cf <- sub("^.*,([^,]*)\\]$", "\\1", v)
  grp_dpar <- strsplit(grp, "__", fixed = TRUE)
  g <- vapply(grp_dpar, `[`, "", 1)
  cf <- ifelse(lengths(grp_dpar) > 1, paste0(vapply(grp_dpar, function(z) z[length(z)], ""), "_", cf), cf)
  out <- list()
  for (gg in unique(g)) {
    i <- which(g == gg)
    levs <- unique(lev[i])
    cfs <- unique(cf[i])
    s <- posterior_summary.pbrfit(object, v[i], probs)
    a <- array(NA_real_, c(length(levs), ncol(s), length(cfs)), dimnames = list(levs, colnames(s), cfs))
    for (k in seq_along(i)) a[lev[i[k]], , cf[i[k]]] <- s[k, ]
    out[[gg]] <- a
  }
  out
}

#' Summarise a fit
#'
#' 'brms''s summary: population-level coefficients, group-level standard
#' deviations and correlations, further distributional parameters, with Rhat
#' and bulk / tail effective sample sizes.
#' @param object A `pbrfit`.
#' @param prob Probability of the credible intervals.
#' @param ... Unused.
#' @return An object of class `summary.pbrfit`.
#' @export
summary.pbrfit <- function(object, prob = 0.95, ...) {
  v <- variables.pbrfit(object)
  hidden <- grepl("^(r_|z_|L_|prior_|Intercept|lp__|lprior)", v)
  fixed <- grep("^b_", v, value = TRUE)
  sds <- grep("^sd_", v, value = TRUE)
  cors <- grep("^cor_", v, value = TRUE)
  rescor <- grep("^rescor__", v, value = TRUE)
  other <- v[!hidden & !v %in% c(fixed, sds, cors, rescor)]
  # sd_<group>__<coef> and cor_<group>__<coef1>__<coef2>; a group label has no "__"
  groups <- list()
  for (s in c(sds, cors)) {
    body <- sub("^(sd|cor)_", "", s)
    g <- sub("__.*$", "", body)
    rest <- substring(body, nchar(g) + 3L)
    lab <- if (startsWith(s, "sd_")) paste0("sd(", rest, ")") else
      paste0("cor(", gsub("__", ",", rest, fixed = TRUE), ")")
    groups[[g]] <- c(groups[[g]], stats::setNames(s, lab))
  }
  grp <- lapply(groups, function(vars) {
    r <- summ_rows(object, unname(vars), prob)
    rownames(r) <- names(vars)
    r
  })
  fr <- summ_rows(object, fixed, prob)
  rownames(fr) <- sub("^b_", "", fixed)
  nlev <- vapply(names(grp), function(g) {
    length(unique(sub("^.*\\[(.*),[^,]*\\]$", "\\1", grep(paste0("^r_", g, "(__[^\\[]*)?\\["), v, value = TRUE))))
  }, 1L)
  ndiv <- if ("divergent__" %in% dimnames(object$stats)$variable) sum(object$stats[, , "divergent__"]) else 0
  structure(list(fit = object, fixed = fr, random = grp, nlevels = nlev,
                 spec = summ_rows(object, other, prob), rescor = summ_rows(object, rescor, prob),
                 ndiv = ndiv, prob = prob),
            class = "summary.pbrfit")
}

fmt_table <- function(m, digits = 2) {
  m <- as.data.frame(m, check.names = FALSE)
  for (k in names(m)) {
    m[[k]] <- if (k %in% c("Bulk_ESS", "Tail_ESS")) {
      ifelse(is.na(m[[k]]), "NA", format(round(m[[k]]), big.mark = ""))
    } else {
      formatC(m[[k]], format = "f", digits = digits)
    }
  }
  print(m, right = TRUE)
}

formula_text <- function(f) {
  if (inherits(f, "pbr_mvformula")) {
    return(paste(vapply(f$parts, formula_text, ""), collapse = "\n         "))
  }
  paste(vapply(f$formulas, function(x) paste(deparse(x, width.cutoff = 500L), collapse = " "), ""),
        collapse = "\n         ")
}

#' @export
print.summary.pbrfit <- function(x, digits = 2, ...) {
  f <- x$fit
  a <- f$args
  cat(" Family:", f$family$family, "\n")
  cat("  Links: mu =", f$family$link, "\n")
  cat("Formula:", formula_text(f$formula), "\n")
  cat("   Data:", f$data_name %||% "data", paste0("(Number of observations: ", f$nobs, ")"), "\n")
  cat("  Draws: ", a$chains, " chains, each with iter = ", a$iter, "; warmup = ", a$warmup, "; thin = ",
      a$thin, ";\n         total post-warmup draws = ", ndraws.pbrfit(f), "\n", sep = "")
  cat("Backend:", if (f$backend == "binary") paste0("compiled model `", f$model$name, "`") else "Julia",
      "\n\n")
  if (length(x$random)) {
    cat("Multilevel Hyperparameters:\n")
    for (g in names(x$random)) {
      cat("~", g, " (Number of levels: ", x$nlevels[[g]], ") \n", sep = "")
      fmt_table(x$random[[g]], digits)
      cat("\n")
    }
  }
  if (nrow(x$fixed)) {
    cat("Regression Coefficients:\n")
    fmt_table(x$fixed, digits)
    cat("\n")
  }
  if (nrow(x$spec)) {
    cat("Further Distributional Parameters:\n")
    fmt_table(x$spec, digits)
    cat("\n")
  }
  if (nrow(x$rescor)) {
    cat("Residual Correlations:\n")
    fmt_table(x$rescor, digits)
    cat("\n")
  }
  cat("Draws were sampled using NUTS (PracticalBayesianRegression.jl). For each parameter, Bulk_ESS\n",
      "and Tail_ESS are effective sample size measures, and Rhat is the potential\n",
      "scale reduction factor on split chains (at convergence, Rhat = 1).\n", sep = "")
  if (x$ndiv > 0) {
    warning("There were ", x$ndiv, " divergent transitions after warmup. Increasing adapt_delta above ",
            f$args$delta, " may help.", call. = FALSE)
  }
  invisible(x)
}

#' @export
print.pbrfit <- function(x, digits = 2, ...) {
  print(summary.pbrfit(x), digits = digits)
  invisible(x)
}

#' Refit a model with new data or settings
#'
#' As 'brms''s `update`. A fit from a compiled model reruns the same
#' executable on `newdata` (no compilation); changing the formula, family or
#' priors fits the changed model as [brm()] would.
#' @param object A `pbrfit`.
#' @param formula. A new formula (`. ~ . + x` works on the main formula).
#' @param newdata New data.
#' @param ... Other [brm()] arguments to change.
#' @return A `pbrfit`.
#' @export
update.pbrfit <- function(object, formula. = NULL, newdata = NULL, ...) {
  a <- list(...)
  args <- object$args
  data <- newdata %||% object$data
  same_model <- is.null(formula.) && !any(c("family", "prior") %in% names(a))
  brm_args <- list(chains = args$chains, iter = args$iter, warmup = args$warmup, thin = args$thin,
                   cores = if (args$threads) args$chains else 1,
                   control = list(adapt_delta = args$delta, max_treedepth = args$max_depth,
                                  metric = args$metric),
                   init = args$init, sample_prior = args$sample_prior, algorithm = args$algorithm)
  brm_args[names(a)] <- a
  if (same_model && object$backend == "binary") {
    fit <- do.call(brm, c(list(formula = object$model, data = data), brm_args))
  } else {
    f <- object$formula
    if (!is.null(formula.)) {
      if (inherits(f, "pbr_mvformula")) stop("update a multivariate formula by passing a new bf() object", call. = FALSE)
      f$formulas[[1]] <- stats::update.formula(f$formulas[[1]], formula.)
    }
    fam <- a$family %||% brmsfamily(object$family$family, object$family$link)
    pri <- if ("prior" %in% names(a)) a$prior else object$prior
    brm_args$family <- NULL
    brm_args$prior <- NULL
    fit <- do.call(brm, c(list(formula = f, data = data, family = fam, prior = pri,
                               backend = if (object$backend == "julia") "julia" else "auto"), brm_args))
  }
  fit$data_name <- if (is.null(newdata)) object$data_name else "newdata"
  fit
}
