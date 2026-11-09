# The model store, shared with PracticalBayesianRegression.jl: one directory
# per compiled model, each with a model.toml describing it (see the Julia
# package's src/store.jl). Reading the store, finding a model and running it
# need no Julia; only compiling does.

#' The model store
#'
#' Where compiled models are kept: the environment variable `PBR_MODEL_STORE`
#' if set, else `~/.pbr/models`. The Julia package uses the same default, so
#' a model compiled from either language is found by both.
#' @return A path.
#' @export
model_store <- function() {
  env <- Sys.getenv("PBR_MODEL_STORE")
  if (nzchar(env)) return(env)
  home <- if (.Platform$OS.type == "windows") Sys.getenv("USERPROFILE") else Sys.getenv("HOME")
  if (!nzchar(home)) home <- path.expand("~")
  file.path(home, ".pbr", "models")
}

store_os <- function() {
  s <- tolower(Sys.info()[["sysname"]])
  switch(s, darwin = "macos", s)
}
store_arch <- function() {
  a <- R.version$arch
  switch(a, arm64 = "aarch64", a)
}

# ------------------------------------------------------- model.toml (subset)

# The top-level `key = value` lines of a TOML file (strings, numbers,
# booleans, arrays of them), up to the first [table]: the part of model.toml
# written for R. Tables (Julia's prior rows) are skipped.
read_toml_top <- function(file) {
  txt <- paste(readLines(file, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  p <- new.env()
  p$s <- txt
  p$i <- 1L
  n <- nchar(txt)
  ch <- function() substr(p$s, p$i, p$i)
  skip_ws <- function(newlines = TRUE) {
    repeat {
      c <- ch()
      if (c %in% c(" ", "\t", "\r") || (newlines && c == "\n")) {
        p$i <- p$i + 1L
      } else if (c == "#") {
        while (p$i <= n && ch() != "\n") p$i <- p$i + 1L
      } else {
        break
      }
    }
  }
  parse_value <- function() {
    c <- ch()
    if (c == "\"") return(parse_string())
    if (c == "'") {
      j <- regexpr("'", substr(p$s, p$i + 1L, n), fixed = TRUE)
      v <- substr(p$s, p$i + 1L, p$i + j - 1L)
      p$i <- p$i + j + 1L
      return(v)
    }
    if (c == "[") {
      p$i <- p$i + 1L
      out <- list()
      repeat {
        skip_ws()
        if (ch() == "]") {
          p$i <- p$i + 1L
          break
        }
        out[[length(out) + 1L]] <- parse_value()
        skip_ws()
        if (ch() == ",") p$i <- p$i + 1L
      }
      return(if (length(out) && all(vapply(out, function(x) is.atomic(x) && length(x) == 1, TRUE)))
        unlist(out) else out)
    }
    rest <- substr(p$s, p$i, n)
    # (a `]` first in a bracket expression is literal in POSIX regex)
    m <- regmatches(rest, regexpr("^[^],\n#]+", rest))
    p$i <- p$i + nchar(m)
    m <- trimws(m)
    if (m == "true") return(TRUE)
    if (m == "false") return(FALSE)
    as.numeric(m)
  }
  parse_string <- function() {
    p$i <- p$i + 1L
    buf <- character(0)
    repeat {
      c <- ch()
      if (c == "\"") {
        p$i <- p$i + 1L
        break
      }
      if (c == "\\") {
        e <- substr(p$s, p$i + 1L, p$i + 1L)
        if (e %in% c("u", "U")) {
          k <- if (e == "u") 4L else 8L
          buf <- c(buf, intToUtf8(strtoi(substr(p$s, p$i + 2L, p$i + 1L + k), 16L)))
          p$i <- p$i + 2L + k
        } else {
          buf <- c(buf, switch(e, n = "\n", t = "\t", r = "\r", b = "\b", f = "\f", e))
          p$i <- p$i + 2L
        }
      } else {
        buf <- c(buf, c)
        p$i <- p$i + 1L
      }
    }
    paste(buf, collapse = "")
  }
  out <- list()
  repeat {
    skip_ws()
    if (p$i > n || ch() == "[") break
    eq <- regexpr("=", substr(p$s, p$i, n), fixed = TRUE)
    key <- trimws(substr(p$s, p$i, p$i + eq - 2L))
    key <- gsub("^\"|\"$", "", key)
    p$i <- p$i + eq
    skip_ws(newlines = FALSE)
    out[[key]] <- parse_value()
  }
  out
}

# ------------------------------------------------------------ compiled models

load_compiled <- function(dir) {
  m <- read_toml_top(file.path(dir, "model.toml"))
  structure(list(name = m$name, dir = normalizePath(dir, winslash = "/"),
                 exe = file.path(normalizePath(dir, winslash = "/"), m$exe),
                 formulas = as.character(m$formulas), family = m$family, link = m$link,
                 columns = as.character(m$columns), keys = as.character(m$keys), meta = m),
            class = "pbr_compiled_model")
}

runs_here <- function(m) identical(m$meta$os, store_os()) && identical(m$meta$arch, store_arch())

#' Compiled models
#'
#' `compiled_models()` lists the compiled models in the store that run on this
#' machine, newest first; index the result by name. `compiled_model()` gets
#' one by name, by the path of its directory, or by the path of its
#' executable. Fit one with `brm(model, data)`, or load it as a function with
#' [load_model()]; delete one with `remove_compiled_model()`. None of these
#' start Julia.
#'
#' @param name Name of a compiled model, the path of its directory, or the
#'   path of its executable.
#' @param store The model store directory.
#' @return `compiled_models()`: a named list of class `pbr_model_store`;
#'   `compiled_model()`: an object of class `pbr_compiled_model`.
#' @seealso [compile_model()], [model_store()]
#' @export
compiled_models <- function(store = model_store()) {
  dirs <- if (dir.exists(store)) list.dirs(store, recursive = FALSE) else character(0)
  dirs <- dirs[!startsWith(basename(dirs), ".") & file.exists(file.path(dirs, "model.toml"))]
  ms <- lapply(dirs, function(d) tryCatch(load_compiled(d), error = function(e) {
    warning("skipping unreadable compiled model ", d, ": ", conditionMessage(e), call. = FALSE)
    NULL
  }))
  ms <- Filter(function(m) !is.null(m) && runs_here(m), ms)
  if (length(ms)) {
    created <- vapply(ms, function(m) as.character(m$meta$created %||% ""), "")
    ms <- ms[order(created, decreasing = TRUE)]
  }
  names(ms) <- vapply(ms, `[[`, "", "name")
  structure(ms, class = "pbr_model_store", store = store)
}

#' @rdname compiled_models
#' @export
compiled_model <- function(name, store = model_store()) {
  dir <- model_dir(name, store)
  if (is.null(dir)) {
    stop("no compiled model `", name, "` in ", store, " (see compiled_models())", call. = FALSE)
  }
  m <- load_compiled(dir)
  if (!runs_here(m)) {
    stop("`", name, "` was compiled for ", m$meta$os, " ", m$meta$arch, "; compile it again here",
         call. = FALSE)
  }
  m
}

# A model's directory from a store name, the directory itself, or any path
# inside it (the executable bin/<name>[.exe]): the nearest ancestor with a
# model.toml, at most three levels up. NULL if none.
model_dir <- function(name, store) {
  if (!grepl("[/\\\\]", name) && !file.exists(name)) {
    d <- file.path(store, name)
    return(if (file.exists(file.path(d, "model.toml"))) d else NULL)
  }
  d <- normalizePath(name, winslash = "/", mustWork = FALSE)
  if (!file.exists(d)) return(NULL)
  if (!dir.exists(d)) d <- dirname(d)
  for (i in 0:3) {
    if (file.exists(file.path(d, "model.toml"))) return(d)
    d <- dirname(d)
  }
  NULL
}

#' Load a compiled model as a function
#'
#' `mymod <- load_model("mymod")` returns a function that fits the compiled
#' model to data: `mymod(df)` is `brm(compiled_model("mymod"), df)`, and any
#' [brm()] sampling argument passes through (`mymod(df, chains = 2, seed =
#' 1)`). Printing it describes the model. No Julia is started, neither to
#' load nor to fit.
#'
#' @param name A model name in the store, the model's directory, or the path
#'   of its executable (`.../mymod/bin/mymod.exe`). The executable needs the
#'   rest of its directory (libraries, `model.toml`), so copy the whole
#'   directory when moving a model.
#' @param store The model store directory.
#' @return A function `function(data, ...)` of class `pbr_model_fn`.
#' @examples
#' \dontrun{
#' compile_model(y ~ x + (1 | g), d, family = bernoulli(), name = "mymod")
#' mymod <- load_model("mymod")
#' fit <- mymod(new_data, seed = 1)
#' }
#' @export
load_model <- function(name, store = model_store()) {
  model <- compiled_model(name, store)
  f <- function(data, ...) brm(model, data, ...)
  structure(f, class = c("pbr_model_fn", "function"), model = model)
}

#' @export
print.pbr_model_fn <- function(x, ...) {
  print(attr(x, "model"))
  cat("Call it with data: fit <- f(data, chains = 4, seed = 1)\n")
  invisible(x)
}

#' @rdname compiled_models
#' @export
remove_compiled_model <- function(name, store = model_store()) {
  m <- compiled_model(name, store)
  unlink(m$dir, recursive = TRUE)
  invisible(NULL)
}

#' @export
print.pbr_compiled_model <- function(x, ...) {
  cat("Compiled model `", x$name, "`\n", sep = "")
  cat("  formula:", paste(x$formulas, collapse = ", "), "\n")
  cat("  family: ", x$family, "(link = ", x$link, ")\n", sep = "")
  cat("  columns:", paste(x$columns, collapse = ", "), "\n")
  cat("  built:  ", x$meta$created %||% "?", " (PracticalBayesianRegression ",
      x$meta$package_version %||% "?", ", Julia ", x$meta$julia_version %||% "?", ")\n", sep = "")
  cat("  exe:    ", x$exe, "\n")
  invisible(x)
}

#' @export
print.pbr_model_store <- function(x, ...) {
  n <- length(x)
  cat(n, " compiled model", if (n != 1) "s", " in ", attr(x, "store"), "\n", sep = "")
  if (n == 0) return(invisible(x))
  df <- data.frame(name = names(x),
                   family = vapply(x, `[[`, "", "family"),
                   formula = vapply(x, function(m) paste(m$formulas, collapse = ", "), ""),
                   built = vapply(x, function(m) as.character(m$meta$created %||% ""), ""),
                   stringsAsFactors = FALSE)
  print(df, row.names = FALSE, right = FALSE)
  invisible(x)
}

#' @export
`[.pbr_model_store` <- function(x, i) structure(unclass(x)[i], class = "pbr_model_store", store = attr(x, "store"))

`%||%` <- function(a, b) if (is.null(a)) b else a

# --------------------------------------------------------------- lookup keys

# This package's lookup key: the Julia text of the formula (spaces removed),
# family and link, prior rows (sorted), and the kind of every data column the
# formula reads. Julia keys a model by its own canonical form; a model
# compiled from R carries both, and Julia adds R's key when it finds a match
# under its own, so the store answers either language without a rebuild.
r_model_key <- function(formula, data, family, prior) {
  f <- as_pbr_formula(formula)
  fam <- resolve_family(family)
  p <- as_pbr_prior(prior)
  nows <- function(s) gsub("[[:space:]]", "", s)
  lines <- c("R", paste0("formula: ", nows(jl_formula(f))),
             paste0("family: ", fam$family, "(", fam$link %||% "", ")"))
  if (nrow(p)) {
    pl <- vapply(seq_len(nrow(p)), function(i) {
      r <- p[i, ]
      paste(r$class, r$coef, r$group, r$resp, r$dpar, r$nlpar, nows(r$prior),
            if (is.na(r$lb)) "" else r$lb, if (is.na(r$ub)) "" else r$ub, sep = "|")
    }, "")
    lines <- c(lines, paste0("prior: ", sort(pl, method = "radix")))
  }
  vars <- intersect(formula_vars(f), names(data))
  if (length(vars)) {
    kinds <- vapply(as.data.frame(data)[vars], column_kind, "")
    lines <- c(lines, paste0("column: ", sort(paste0(vars, " (", kinds, ")"), method = "radix")))
  }
  paste0(paste(lines, collapse = "\n"), "\n")
}

# a model binary's view of a column
column_kind <- function(x) if (is.numeric(x)) "numeric" else "factor"

#' Find a compiled model
#'
#' The compiled model for this formula, family, priors and data column types,
#' if one is in the store. Does not start Julia.
#' @inheritParams brm
#' @param store The model store directory.
#' @return A `pbr_compiled_model`, or `NULL`.
#' @export
find_compiled_model <- function(formula, data, family = gaussian(), prior = NULL, store = model_store()) {
  key <- r_model_key(formula, data, family, prior)
  for (m in compiled_models(store)) if (key %in% m$keys) return(m)
  NULL
}

#' Compile a model
#'
#' Compile a model into a standalone executable and keep it in the model
#' store, or return the stored one if this formula, family, priors and column
#' types have been compiled before (by this package or the Julia package).
#' Fit it to any data with the right columns with `brm(model, data)`, which
#' needs no Julia. `brm(..., backend = "binary")` calls this for you.
#'
#' Compiling needs Julia (see [pbr_setup()]), a C compiler and JuliaC: by
#' default the `juliac` app on the PATH, or a command given by
#' `options(pbr.juliac = ...)` or the environment variable `PBR_JULIAC`, e.g.
#' `julia --project=/path/env -e "using JuliaC; JuliaC.main(ARGS)" --`. It
#' takes a few minutes. Not every model can be compiled yet: unsupported
#' features (smooths, ordinal families, censoring, ...) are named in the error.
#'
#' @inheritParams brm
#' @param name Name for the model in the store; default
#'   `<response>_<family>_<hash>`.
#' @param force Compile again even if the model is in the store.
#' @param juliac Command line for JuliaC.
#' @param store The model store directory.
#' @param verbose Show Julia's build output.
#' @return A `pbr_compiled_model`.
#' @export
compile_model <- function(formula, data, family = gaussian(), prior = NULL, name = NULL, force = FALSE,
                          juliac = getOption("pbr.juliac", Sys.getenv("PBR_JULIAC")),
                          store = model_store(), verbose = TRUE) {
  key <- r_model_key(formula, data, family, prior)
  if (!force) {
    for (m in compiled_models(store)) {
      if (key %in% m$keys && (is.null(name) || m$name == name)) return(m)
    }
  }
  src <- sprintf(paste0("PBRR.compile(%s, %s; family = %s, prior = %s, name = %s, store = %s, ",
                        "force = %s, juliac = %s, key = %s, verbose = %s)"),
                 jl_formula(formula), jl_data(data, formula_vars(formula)), jl_family(family),
                 jl_prior(prior), jl_str(name %||% ""), jl_str(julia_path(store)), tolower(isTRUE(force)),
                 jl_str(juliac %||% ""), jl_str(key), tolower(isTRUE(verbose)))
  nm <- jl_eval(src)
  compiled_model(nm, store)
}
