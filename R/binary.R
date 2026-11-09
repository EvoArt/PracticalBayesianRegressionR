# Running a compiled model: write the data as CSV, run the executable, read
# its draws back. No Julia session is involved.

model_columns <- function(model) {
  nm <- sub(" \\((numeric|factor)\\)$", "", model$columns)
  kind <- sub("^.* \\((numeric|factor)\\)$", "\\1", model$columns)
  stats::setNames(kind, nm)
}

check_model_data <- function(model, data) {
  cols <- model_columns(model)
  miss <- setdiff(names(cols), names(data))
  if (length(miss)) {
    stop("the data lack column(s) the compiled model `", model$name, "` needs: ",
         paste(miss, collapse = ", "), call. = FALSE)
  }
  for (nm in names(cols)[cols == "numeric"]) {
    if (!is.numeric(data[[nm]])) {
      stop("column `", nm, "` must be numeric for the compiled model `", model$name, "`", call. = FALSE)
    }
  }
  invisible(TRUE)
}

# set environment variables for the child process only, restoring R's own
with_env <- function(set, unset, code) {
  old <- Sys.getenv(c(names(set), unset), unset = NA)
  on.exit({
    for (nm in names(old)) {
      if (is.na(old[[nm]])) Sys.unsetenv(nm) else do.call(Sys.setenv, stats::setNames(list(old[[nm]]), nm))
    }
  })
  if (length(set)) do.call(Sys.setenv, as.list(set))
  if (length(unset)) Sys.unsetenv(unset)
  force(code)
}

# CSV of draws (chains stacked, `chain` and `iteration` first) ->
# iterations x chains x variables array
read_draws_csv <- function(file, chains) {
  d <- utils::read.csv(file, check.names = FALSE)
  vars <- setdiff(names(d), c("chain", "iteration"))
  niter <- nrow(d) %/% chains
  m <- as.matrix(d[vars])
  storage.mode(m) <- "double"
  array(m, c(niter, chains, length(vars)),
        dimnames = list(iteration = NULL, chain = NULL, variable = vars))
}

#' Run a compiled model's executable
#'
#' The low-level runner behind `brm(model, data)`: returns the draws, sampler
#' statistics and unconstrained draws as arrays (iterations x chains x
#' variables), the chain times and what the executable printed.
#' @param model A compiled model ([compiled_model()]).
#' @param data A data frame with the columns the model needs.
#' @param chains,iter,warmup,seed Sampling settings, as in [brm()].
#' @param delta Target acceptance rate (`adapt_delta`).
#' @param max_depth Maximum tree depth (`max_treedepth`).
#' @param metric `"diag"` or `"dense"` mass matrix.
#' @param init Initial values uniform(-init, init) on the unconstrained scale.
#' @param threads Run chains in parallel threads.
#' @return A list.
#' @export
run_compiled <- function(model, data, chains = 4, iter = 2000, warmup = floor(iter / 2), seed = NA,
                         delta = 0.8, max_depth = 10, metric = "diag", init = 2, threads = TRUE) {
  if (is.character(model)) model <- compiled_model(model)
  data <- as.data.frame(data)
  check_model_data(model, data)
  if (!file.exists(model$exe)) stop("the executable of `", model$name, "` is missing: ", model$exe, call. = FALSE)
  if (is.na(seed)) seed <- sample.int(.Machine$integer.max, 1)
  dir <- tempfile("pbr_run_")
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  datafile <- file.path(dir, "data.csv")
  drawsfile <- file.path(dir, "draws.csv")
  ufile <- file.path(dir, "theta.csv")
  write_data_csv(data, datafile, names(model_columns(model)))
  args <- c("--data", datafile, "--output", drawsfile, "--unconstrained", ufile,
            "--chains", chains, "--iter", iter, "--warmup", warmup, "--seed", format(seed, scientific = FALSE),
            "--delta", delta, "--max-depth", max_depth, "--metric", metric, "--init", init, "--no-summary")
  if (!threads) args <- c(args, "--no-threads")
  nthreads <- if (threads) min(chains, max(1L, parallel_cores())) else 1L
  # a bundle carries its own libraries; an inherited LD_LIBRARY_PATH (a
  # cluster's compiler module) can outrank them and break start-up
  unset <- if (store_os() == "linux") "LD_LIBRARY_PATH" else character(0)
  t0 <- Sys.time()
  out <- with_env(c(JULIA_NUM_THREADS = as.character(nthreads)), unset,
                  suppressWarnings(system2(model$exe, shQuote(args), stdout = TRUE, stderr = TRUE)))
  wall <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  status <- attr(out, "status")
  if (!is.null(status) && status != 0) {
    stop("compiled model `", model$name, "` failed (exit code ", status, "):\n", paste(out, collapse = "\n"),
         call. = FALSE)
  }
  times <- rep(NA_real_, chains)
  tl <- grep("^chain times \\(s\\): ", out, value = TRUE)
  if (length(tl)) times <- as.numeric(strsplit(sub("^chain times \\(s\\): ", "", tl[1]), ",\\s*")[[1]])
  nl <- grep("^fitting [0-9]+ observations", out, value = TRUE)
  nobs <- if (length(nl)) as.integer(sub("^fitting ([0-9]+) observations.*", "\\1", nl[1])) else NA_integer_
  all <- read_draws_csv(drawsfile, chains)
  vars <- dimnames(all)$variable
  st <- endsWith(vars, "__")
  list(draws = all[, , !st, drop = FALSE], stats = all[, , st, drop = FALSE],
       theta = read_draws_csv(ufile, chains), times = times, wall = wall, nobs = nobs, seed = seed,
       output = out)
}
