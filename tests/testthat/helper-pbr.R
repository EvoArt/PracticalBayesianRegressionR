# Two tiers. Tier 1 (most files) needs no Julia. Tier 2 needs a Julia
# session; it is skipped unless PBR_TEST_JULIA=1, and the compile tests also
# need PBR_TEST_JULIAC=1 (they build an executable: minutes, several GB).

pbr <- asNamespace("PracticalBayesianRegressionR")

# tier 1 must never start Julia: a store lookup that misses would otherwise
# fall back to the Julia back end (and install it)
if (!identical(Sys.getenv("PBR_TEST_JULIA"), "1")) options(pbr.julia = FALSE)

skip_without_julia <- function() {
  testthat::skip_if_not(identical(Sys.getenv("PBR_TEST_JULIA"), "1"), "PBR_TEST_JULIA not set")
  testthat::skip_if_not_installed("JuliaCall")
}

skip_without_juliac <- function() {
  skip_without_julia()
  testthat::skip_if_not(identical(Sys.getenv("PBR_TEST_JULIAC"), "1"), "PBR_TEST_JULIAC not set")
}

sim_data <- function(n = 120, seed = 1) {
  set.seed(seed)
  g <- sprintf("g%02d", sample(1:8, n, replace = TRUE))
  x <- rnorm(n)
  u <- rnorm(8, 0, 0.5)[as.integer(factor(g))]
  data.frame(y = 1 + 0.5 * x + u + rnorm(n, 0, 0.7), x = x, g = g,
             k = rpois(n, exp(0.3 + 0.4 * x)), stringsAsFactors = FALSE)
}

# a store with one model directory per spec, whose "executable" is a shell
# script that writes canned output: the binary code path without Julia
fake_model <- function(store, name, keys = character(0), os = pbr$store_os(), arch = pbr$store_arch(),
                       columns = c("g (factor)", "x (numeric)", "y (numeric)"), exe = TRUE,
                       priors = character(0), created = "2026-11-09 12:00:00") {
  dir <- file.path(store, name)
  dir.create(file.path(dir, "bin"), recursive = TRUE, showWarnings = FALSE)
  # TOML basic strings: escape " and newlines (literal replacements: fixed)
  q <- function(s) {
    s <- gsub("\"", paste0("\\", "\""), s, fixed = TRUE)
    paste0("\"", gsub("\n", paste0("\\", "n"), s, fixed = TRUE), "\"")
  }
  arr <- function(v) paste0("[", paste(vapply(v, q, ""), collapse = ", "), "]")
  toml <- c("# Written by PracticalBayesianRegression.compile_model.",
            paste0("arch = ", q(arch)), paste0("columns = ", arr(columns)), paste0("created = ", q(created)),
            paste0("exe = ", q(paste0("bin/", name))), "family = \"gaussian\"",
            "formula = \"y~x+(1|g)\"", "formulas = [\"y~x+(1|g)\"]", "julia_version = \"1.12.7\"",
            paste0("keys = ", arr(keys)), "link = \"identity\"", paste0("name = ", q(name)),
            paste0("os = ", q(os)), "package_version = \"0.1.0\"", paste0("priors = ", arr(priors)), "",
            "[[prior]]", "class = \"b\"", "prior = \"normal(0, 2)\"")
  writeLines(toml, file.path(dir, "model.toml"))
  if (exe) {
    f <- file.path(dir, "bin", name)
    writeLines(c(
      "#!/bin/sh",
      "# fake model binary: 2 chains x 3 draws, sigma and one coefficient",
      "while [ $# -gt 0 ]; do",
      "  case \"$1\" in",
      "    --output) out=\"$2\"; shift ;;",
      "    --unconstrained) unc=\"$2\"; shift ;;",
      "    --data) data=\"$2\"; shift ;;",
      "  esac; shift",
      "done",
      "n=$(($(wc -l < \"$data\") - 1))",
      "echo \"fitting $n observations, 2 parameters, 2 chains x 6 iterations (3 warm-up), seed 1\"",
      "printf 'chain,iteration,lp__,accept_stat__,stepsize__,treedepth__,n_leapfrog__,divergent__,energy__,b_x,sigma\n' > \"$out\"",
      "for c in 1 2; do for i in 1 2 3; do",
      "  echo \"$c,$i,-1.5,0.9,0.5,2,3,0,2.0,0.$c$i,1.$i\" >> \"$out\"",
      "done; done",
      "printf 'chain,iteration,theta[1],theta[2]\n' > \"$unc\"",
      "for c in 1 2; do for i in 1 2 3; do echo \"$c,$i,0.$c$i,0.0$i\" >> \"$unc\"; done; done",
      "echo 'chain times (s): 0.25, 0.5'",
      "echo \"draws written to $out\""), f)
    Sys.chmod(f, "755")
  }
  dir
}
