#' Fit a Bayesian regression model
#'
#' Fit a generalised (non-)linear multilevel model with the syntax of
#' 'brms'. The model is fitted by PracticalBayesianRegression.jl, either in a
#' Julia session (`backend = "julia"`) or by a compiled executable from the
#' model store (`backend = "binary"`), which needs no Julia at all once the
#' model has been compiled.
#'
#' With the default `backend = "auto"`, a model that is already in the store
#' (compiled earlier from R or from Julia, for the same formula, family,
#' priors and column types) runs from its executable; anything else runs in
#' Julia. `backend = "binary"` compiles the model first if needed (minutes,
#' once; see [compile_model()]). `brm(model, data)` with a compiled model
#' from [compiled_model()] runs it directly.
#'
#' Fits by both back ends are the same: given the same `seed`, the draws are
#' identical.
#'
#' @param formula A model formula: a [bf()] object, an R formula or a string;
#'   or a compiled model ([compiled_model()]).
#' @param data A data frame.
#' @param family The response distribution and link: see [brmsfamily()].
#' @param prior Priors: see [set_prior()]. `NULL` uses 'brms''s defaults.
#' @param chains Number of Markov chains.
#' @param iter Iterations per chain, including warm-up.
#' @param warmup Warm-up iterations per chain.
#' @param thin Thinning of the saved draws.
#' @param cores Chains run in parallel when `cores > 1`.
#' @param control A list: `adapt_delta` (0.8), `max_treedepth` (10), `metric`
#'   (`"diag"` or `"dense"`).
#' @param seed Random seed; `NA` for a random one.
#' @param init Initial values are drawn uniformly from `(-init, init)` on the
#'   unconstrained scale (default 2, as in Stan); `0` starts at zero.
#' @param backend `"auto"`, `"julia"` or `"binary"`.
#' @param file Optional file name (`.rds` added): the fit is saved there and
#'   read back instead of refitting, as in 'brms'.
#' @param file_refit When to refit a model saved in `file`: `"never"`,
#'   `"on_change"` (formula, family, priors, data or settings changed) or
#'   `"always"`.
#' @param silent Print less (1) or nothing (2).
#' @param sample_prior `"no"` or `"yes"` (also sample from the priors; Julia
#'   back end only).
#' @param algorithm `"sampling"` (NUTS) or an approximation (`"laplace"`,
#'   `"meanfield"`, `"fullrank"`, `"pathfinder"`; Julia back end only).
#' @param ... Not used; 'brms' arguments without an equivalent here are
#'   reported.
#' @return An object of class `pbrfit`.
#' @examples
#' \dontrun{
#' fit <- brm(count ~ zAge + zBase * Trt + (1 | patient), data = epilepsy,
#'            family = poisson())
#' summary(fit)
#' }
#' @export
brm <- function(formula, data, family = gaussian(), prior = NULL, chains = 4, iter = 2000,
                warmup = floor(iter / 2), thin = 1, cores = getOption("mc.cores", chains), control = NULL,
                seed = NA, init = NULL, backend = getOption("pbr.backend", "auto"), file = NULL,
                file_refit = getOption("brms.file_refit", "never"), silent = 1, sample_prior = "no",
                algorithm = "sampling", ...) {
  dots <- list(...)
  if (length(dots)) {
    warning("ignored argument(s) without an equivalent here: ", paste(names(dots), collapse = ", "),
            call. = FALSE)
  }
  # a name for the summary; never deparse a data frame passed by value (do.call)
  ds <- substitute(data)
  data_name <- if (is.symbol(ds) || (is.call(ds) && length(all.names(ds)) <= 20)) {
    deparse(ds, width.cutoff = 60L)[1]
  } else {
    "data"
  }
  if (inherits(formula, "pbr_model_fn")) formula <- attr(formula, "model")    # from load_model()
  if (inherits(formula, "pbr_compiled_model")) {
    model <- formula
    formula <- bf(model$formulas[[1]])
    if (length(model$formulas) > 1) formula <- do.call(bf, c(list(model$formulas[[1]]), as.list(model$formulas[-1])))
    family <- brmsfamily(model$family, model$link)
    prior <- model_priors(model)
    backend <- "binary"
  } else {
    model <- NULL
    formula <- as_pbr_formula(formula)
  }
  backend <- match.arg(backend, c("auto", "julia", "binary"))
  file_refit <- match.arg(file_refit, c("never", "on_change", "always"))
  stopifnot(iter > warmup, warmup >= 0, chains >= 1, thin >= 1)
  if (is.na(seed)) seed <- sample.int(.Machine$integer.max, 1)
  args <- list(chains = as.integer(chains), iter = as.integer(iter), warmup = as.integer(warmup),
               thin = as.integer(thin), seed = seed, threads = cores > 1,
               delta = control$adapt_delta %||% 0.8, max_depth = as.integer(control$max_treedepth %||% 10),
               metric = control$metric %||% "diag", init = init_value(init),
               sample_prior = sample_prior, algorithm = algorithm)
  data <- as.data.frame(data)

  if (!is.null(file)) {
    if (!grepl("\\.rds$", file)) file <- paste0(file, ".rds")
    if (file.exists(file) && file_refit != "always") {
      old <- readRDS(file)
      same <- inherits(old, "pbrfit") &&
        identical(old$key, r_model_key(formula, data, family, prior)) &&
        identical(old$data, data[intersect(names(data), formula_vars(formula))]) &&
        identical(old$args[setdiff(names(args), "seed")], args[setdiff(names(args), "seed")])
      if (file_refit == "never" || same) return(old)
    }
  }

  binary_ok <- args$sample_prior == "no" && args$algorithm == "sampling"
  if (is.null(model) && backend != "julia" && binary_ok) {
    model <- if (backend == "binary") {
      compile_model(formula, data, family, prior, verbose = silent < 2)
    } else {
      find_compiled_model(formula, data, family, prior)
    }
  }
  if (!is.null(model) && !binary_ok) {
    stop("sample_prior and approximate algorithms need the Julia back end (backend = \"julia\")", call. = FALSE)
  }
  fit <- if (!is.null(model)) {
    fit_binary(model, formula, data, family, prior, args, silent)
  } else {
    fit_julia(formula, data, family, prior, args)
  }
  fit$data_name <- data_name
  if (!is.null(file)) saveRDS(fit, file)
  fit
}

init_value <- function(init) {
  if (is.null(init) || identical(init, "random")) return(2)
  if (is.numeric(init) && length(init) == 1) return(as.numeric(init))
  if (identical(init, "0")) return(0)
  stop("`init` must be a single number (the range of random initial values) or \"random\"", call. = FALSE)
}

# model.toml's `priors`: "class|coef|group|resp|dpar|nlpar|prior|lb|ub" rows
model_priors <- function(model) {
  pl <- as.character(model$meta$priors)
  if (!length(pl)) return(NULL)
  parts <- strsplit(pl, "|", fixed = TRUE)
  do.call(c, lapply(parts, function(p) {
    p <- c(p, rep("", 9 - length(p)))
    set_prior(p[7], class = p[1], coef = p[2], group = p[3], resp = p[4], dpar = p[5], nlpar = p[6],
              lb = if (nzchar(p[8])) p[8] else NA, ub = if (nzchar(p[9])) p[9] else NA)
  }))
}

as_draws_arr <- function(a) posterior::as_draws_array(a)

fit_binary <- function(model, formula, data, family, prior, args, silent) {
  if (silent < 2) message("Running compiled model `", model$name, "` (no compilation needed)")
  r <- run_compiled(model, data, chains = args$chains, iter = args$iter, warmup = args$warmup, seed = args$seed,
                    delta = args$delta, max_depth = args$max_depth, metric = args$metric, init = args$init,
                    threads = args$threads)
  if (silent < 1) cat(r$output, sep = "\n")
  keep <- seq(1, dim(r$draws)[1], by = args$thin)
  sub3 <- function(a) a[keep, , , drop = FALSE]
  new_pbrfit(formula = formula, family = family, prior = prior, data = data, args = args,
             draws = sub3(r$draws), stats = sub3(r$stats), theta = sub3(r$theta), times = r$times,
             nobs = r$nobs, backend = "binary", model = model)
}

fit_julia <- function(formula, data, family, prior, args) {
  call <- sprintf(paste0("PBR.brm(%s, %s; family = %s, prior = %s, chains = %d, iter = %d, warmup = %d, ",
                         "thin = %d, seed = %s, delta = %s, max_depth = %d, metric = :%s, init = %s, ",
                         "threads = %s, sample_prior = :%s, algorithm = :%s)"),
                  jl_formula(formula), jl_data(data, formula_vars(formula)), jl_family(family), jl_prior(prior),
                  args$chains, args$iter, args$warmup, args$thin, jl_num(args$seed), jl_num(args$delta),
                  args$max_depth, args$metric, jl_num(args$init), tolower(args$threads), args$sample_prior,
                  args$algorithm)
  jfit <- jl_eval(call, "Julia")
  fit <- fit_from_julia(jfit, formula, data, family, prior, args)
  fit$julia_call <- call
  fit
}

fit_from_julia <- function(jfit, formula, data, family, prior, args) {
  jc <- function(f, ...) JuliaCall::julia_call(f, ..., need_return = "R")
  vars <- jc("PBRR.variable_names", jfit)
  draws <- jc("PBRR.draws_array", jfit, vars)
  stats <- jc("PBRR.stats_array", jfit)
  snames <- jc("PBRR.stat_names", jfit)
  theta <- jc("PBRR.theta_array", jfit)
  dn <- function(a, v) {
    a <- array(a, dim = c(dim(a)[1], dim(a)[2], length(v)))
    dimnames(a) <- list(iteration = NULL, chain = NULL, variable = v)
    a
  }
  fit <- new_pbrfit(formula = formula, family = family, prior = prior, data = data, args = args,
                    draws = dn(draws, vars), stats = dn(stats, snames),
                    theta = dn(theta, paste0("theta[", seq_len(dim(theta)[3]), "]")),
                    times = as.numeric(jc("PBRR.chain_times", jfit)),
                    nobs = jc("PBRR.nobs", jfit), backend = "julia", model = NULL)
  fit$cache$julia <- jfit
  fit$cache$session <- .pbr$session
  fit
}

new_pbrfit <- function(formula, family, prior, data, args, draws, stats, theta, times, nobs, backend, model) {
  data <- data[intersect(names(data), formula_vars(formula))]
  structure(list(formula = formula, family = resolve_family(family), prior = as_pbr_prior(prior),
                 data = data, args = args, draws = draws, stats = stats, theta = theta, times = times,
                 nobs = nobs, backend = backend, model = model,
                 key = r_model_key(formula, data, family, prior),
                 cache = new.env(parent = emptyenv())),
            class = "pbrfit")
}
