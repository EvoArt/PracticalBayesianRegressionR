# The Julia session: one per R process, started on first need. Fitting a
# compiled model never needs it (see store.R / binary.R).

.pbr <- new.env(parent = emptyenv())
.pbr$ready <- FALSE
.pbr$session <- NULL

#' Start Julia and load PracticalBayesianRegression.jl
#'
#' Called automatically the first time something needs Julia (fitting with
#' `backend = "julia"`, compiling a model, predictions). Call it yourself to
#' choose the Julia installation, the Julia project or the number of threads.
#'
#' The Julia project defaults to the environment variable
#' `PBR_JULIA_PROJECT`, else a project in this package's user data directory
#' that is created on first use. It installs the registered
#' PracticalBayesianRegression.jl (a version this package supports, see
#' [pbr_julia_compat()]); with `dev = path` it uses your local checkout instead
#' (`Pkg.develop`). The first set-up downloads and precompiles the Julia
#' packages, which takes a few minutes.
#'
#' @param julia_home Directory containing the `julia` executable; `NULL`
#'   lets JuliaCall find it (it respects `JULIA_HOME`).
#' @param project Julia project directory to activate.
#' @param dev Path of a local PracticalBayesianRegression.jl checkout to use.
#' @param threads Julia threads, so that chains run in parallel. Only takes
#'   effect if Julia has not started yet in this R session; default
#'   `JULIA_NUM_THREADS` if set, else up to 4.
#' @param update Update PracticalBayesianRegression.jl to the newest supported
#'   version.
#' @param verbose Show Julia's set-up output.
#' @param force Run the set-up again even if Julia is already up.
#' @return Invisibly `TRUE`.
#' @export
pbr_setup <- function(julia_home = NULL, project = NULL, dev = getOption("pbr.julia_dev"),
                      threads = NULL, update = FALSE, verbose = FALSE, force = FALSE) {
  if (.pbr$ready && !force && !update) return(invisible(TRUE))
  if (isFALSE(getOption("pbr.julia"))) {
    stop("this needs Julia, which is switched off (options(pbr.julia = FALSE))", call. = FALSE)
  }
  if (!requireNamespace("JuliaCall", quietly = TRUE)) {
    stop("this needs Julia through the JuliaCall package: install.packages(\"JuliaCall\"). ",
         "(Compiled models run without it: see compiled_models().)", call. = FALSE)
  }
  if (!nzchar(Sys.getenv("JULIA_NUM_THREADS"))) {
    if (is.null(threads)) threads <- min(4L, max(1L, parallel_cores()))
    Sys.setenv(JULIA_NUM_THREADS = as.character(threads))
  }
  if (is.null(project)) project <- pbr_julia_project()
  JuliaCall::julia_setup(JULIA_HOME = julia_home, verbose = verbose)
  jc <- function(src) JuliaCall::julia_command(src)
  jc(sprintf("import Pkg; Pkg.activate(%s; io = devnull)", jl_str(julia_path(project))))
  # Pkg operations are slow (Pkg.instantiate alone was 38 s a session on a
  # cluster file system), so run them only when something must change: the
  # package missing, or a dev checkout other than the one in use
  where <- JuliaCall::julia_eval(paste0(
    "let id = Base.identify_package(\"PracticalBayesianRegression\"); ",
    "p = id === nothing ? nothing : Base.locate_package(id); ",
    "p === nothing ? \"\" : dirname(dirname(p)) end"))
  changed <- FALSE
  if (!is.null(dev)) {
    dev <- normalizePath(dev, winslash = "/", mustWork = TRUE)
    if (!identical(normalizePath(where, winslash = "/", mustWork = FALSE), dev)) {
      message("Julia: using PracticalBayesianRegression.jl from ", dev)
      jc(sprintf("Pkg.develop(path = %s)", jl_str(julia_path(dev))))
      changed <- TRUE
    }
  } else if (!nzchar(where)) {
    message("Julia: installing PracticalBayesianRegression.jl ", pbr_julia_compat(), " into ", project,
            " (once; downloading and precompiling takes a few minutes)")
    url <- getOption("pbr.julia_url")
    if (is.null(url)) {
      # the registered package, held to the versions this R package supports
      jc(sprintf("Pkg.add(Pkg.PackageSpec(name = \"PracticalBayesianRegression\", version = %s))",
                 jl_str(pbr_julia_compat())))
      jc(sprintf("Pkg.compat(\"PracticalBayesianRegression\", %s)", jl_str(pbr_julia_compat())))
    } else {
      jc(sprintf("Pkg.add(url = %s)", jl_str(url)))
    }
    changed <- TRUE
  } else if (update) {
    jc("Pkg.update(\"PracticalBayesianRegression\")")
    changed <- TRUE
  }
  if (changed) jc("Pkg.instantiate(; io = devnull)")
  loaded <- tryCatch({
    jc("import PracticalBayesianRegression")
    TRUE
  }, error = function(e) FALSE)
  if (!loaded) {
    # e.g. a manifest whose packages are not installed yet
    jc("Pkg.instantiate(; io = devnull)")
    jc("import PracticalBayesianRegression")
  }
  jc("const PBR = PracticalBayesianRegression")
  # an older install (e.g. from before this R version) lacks what the glue
  # code uses: say how to update rather than fail somewhere obscure
  ver <- JuliaCall::julia_eval("string(pkgversion(PracticalBayesianRegression))")
  if (package_version(ver) < package_version(.pbr_julia_min)) {
    stop("PracticalBayesianRegression.jl ", ver, " is too old for this R package (needs >= ", .pbr_julia_min,
         "): run pbr_setup(update = TRUE) in a new R session", call. = FALSE)
  }
  .pbr$julia_version <- ver
  jc(sprintf("include(%s)", jl_str(julia_path(glue_file()))))
  .pbr$session <- paste0(Sys.getpid(), "-", format(Sys.time(), "%H%M%OS6"))
  .pbr$ready <- TRUE
  invisible(TRUE)
}

# the PracticalBayesianRegression.jl versions this package supports: the
# model store and compiled models arrived in 0.2.0
.pbr_julia_min <- "0.2.0"

#' The PracticalBayesianRegression.jl versions this package supports
#'
#' A Julia compat specification (`"0.2"` means 0.2.0 up to, not including,
#' 0.3.0), used when [pbr_setup()] installs the Julia package.
#' @return A string.
#' @export
pbr_julia_compat <- function() "0.2"

#' The Julia project PracticalBayesianRegressionR activates
#'
#' `PBR_JULIA_PROJECT` if set; otherwise a directory in the package's user
#' data directory.
#' @return A path.
#' @export
pbr_julia_project <- function() {
  env <- Sys.getenv("PBR_JULIA_PROJECT")
  if (nzchar(env)) return(env)
  p <- file.path(tools::R_user_dir("PracticalBayesianRegressionR", "data"), "julia")
  dir.create(p, recursive = TRUE, showWarnings = FALSE)
  p
}

glue_file <- function() {
  p <- system.file("julia", "glue.jl", package = "PracticalBayesianRegressionR")
  if (!nzchar(p)) p <- file.path(getwd(), "inst", "julia", "glue.jl")
  p
}

parallel_cores <- function() {
  n <- tryCatch(parallel::detectCores(logical = FALSE), error = function(e) NA)
  if (is.na(n)) 1L else as.integer(n)
}

#' Is Julia running?
#' @return `TRUE` or `FALSE`.
#' @export
pbr_julia_ready <- function() isTRUE(.pbr$ready)

# Julia string literals want forward slashes on Windows
julia_path <- function(p) gsub("\\\\", "/", p)

jl_eval <- function(src, need_return = c("R", "Julia", "None")) {
  need_return <- match.arg(need_return)
  pbr_setup()
  JuliaCall::julia_eval(src, need_return = need_return)
}

# ---------------------------------------------------------------- data

# R column type -> what the Julia side rebuilds (glue.jl read_table)
col_kind <- function(x) {
  if (is.factor(x) || is.character(x)) return("string")
  if (is.logical(x)) return("logical")
  if (is.integer(x)) return("integer")
  if (is.numeric(x)) return("double")
  stop("data columns must be numeric, integer, logical, character or factor; got ", class(x)[1],
       call. = FALSE)
}

# The data as a Julia table, through a CSV file written exactly (%.17g), so
# what Julia fits is bit for bit what R holds. Only the columns the formula
# reads are sent. Returns the Julia source of an expression for the table.
jl_data <- function(data, vars) {
  data <- as.data.frame(data)
  cols <- intersect(names(data), vars)
  warn_factor_levels(data[cols])
  file <- tempfile("pbr_data_", fileext = ".csv")
  write_data_csv(data, file, cols)
  kinds <- vapply(data[cols], col_kind, "")
  sprintf("PBRR.read_table(%s, %s, %s)", jl_str(julia_path(file)), jl_strvec(cols), jl_strvec(kinds))
}

# Julia (like R's factor()) sorts the levels of a text column; a factor with
# another level order would be coded with a different reference level
warn_factor_levels <- function(data) {
  for (nm in names(data)) {
    x <- data[[nm]]
    if (is.factor(x)) {
      lv <- levels(droplevels(x))
      if (!identical(lv, sort(lv, method = "radix"))) {
        warning("factor `", nm, "` has levels in a non-alphabetical order; the Julia back end sorts ",
                "levels, so the reference level is `", sort(lv, method = "radix")[1],
                "` (rename the levels to choose another)", call. = FALSE)
      }
    }
  }
}

# One data column as CSV cells: numbers exactly, text quoted where needed
csv_cells <- function(x) {
  if (is.factor(x)) x <- as.character(x)
  out <- if (is.logical(x)) {
    ifelse(x, "TRUE", "FALSE")
  } else if (is.integer(x)) {
    as.character(x)
  } else if (is.numeric(x)) {
    sprintf("%.17g", x)
  } else {
    csv_quote(as.character(x))
  }
  out[is.na(x)] <- "NA"
  out
}

csv_quote <- function(s) {
  q <- grepl("[\",\r\n]", s) | (!is.na(s) & s == "")
  s[q] <- paste0("\"", gsub("\"", "\"\"", s[q], fixed = TRUE), "\"")
  s
}

write_data_csv <- function(data, file, columns = names(data)) {
  cells <- lapply(columns, function(nm) csv_cells(data[[nm]]))
  body <- if (length(columns) == 0) character(0) else do.call(paste, c(cells, sep = ","))
  con <- file(file, "wb")
  on.exit(close(con))
  writeLines(c(paste(csv_quote(columns), collapse = ","), body), con, sep = "\n", useBytes = TRUE)
  invisible(file)
}
