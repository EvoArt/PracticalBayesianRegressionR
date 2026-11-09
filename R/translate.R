# R -> Julia source text.
#
# Models reach Julia as source code (`PBR.bf("y ~ x + (1 | g)")`, ...), not as
# JuliaCall-converted objects: the text is exactly what a Julia user would
# type, it can be shown when something goes wrong (`fit$julia_call`), and no
# conversion rule of JuliaCall's (length-one vectors become scalars, integers
# Int32, ...) can change what the model means.

# A Julia string literal: `$` would interpolate, `\` and `"` need escapes.
jl_str <- function(x) {
  x <- gsub("\\", "\\\\", x, fixed = TRUE)
  x <- gsub("\"", "\\\"", x, fixed = TRUE)
  x <- gsub("$", "\\$", x, fixed = TRUE)
  x <- gsub("\n", "\\n", x, fixed = TRUE)
  paste0("\"", x, "\"")
}

jl_strvec <- function(x) paste0("String[", paste(vapply(x, jl_str, ""), collapse = ", "), "]")

jl_num <- function(x) {
  if (is.na(x)) stop("NA is not a number Julia can use here", call. = FALSE)
  if (is.infinite(x)) return(if (x > 0) "Inf" else "-Inf")
  if (x == round(x) && abs(x) < 1e15) return(format(x, scientific = FALSE))
  format(x, digits = 17)
}

# operators written infix; R's parse tree keeps every `(`, so the text has the
# grouping R gave it, and the Julia front end (which accepts R formulas as
# strings) restores R's precedence for `:`
.infix <- c("+", "-", "*", "/", "^", ":", "|", "||", "~", "==", "!=", "<", ">", "<=", ">=", "&", "&&")

#' Julia source for one R expression from a formula.
#' @noRd
jl_expr <- function(e) {
  if (is.symbol(e)) {
    nm <- as.character(e)
    if (!grepl("^[A-Za-z_][A-Za-z0-9_]*$", nm)) {
      stop("variable name `", nm, "` cannot be used in a formula here: Julia names may contain ",
           "only letters, digits and _ (rename the column, e.g. ", gsub("[^A-Za-z0-9_]", "_", nm), ")",
           call. = FALSE)
    }
    return(nm)
  }
  if (is.character(e)) return(jl_str(e))
  if (is.logical(e)) {
    if (is.na(e)) stop("NA in a formula", call. = FALSE)
    return(if (e) "true" else "false")
  }
  if (is.numeric(e)) return(jl_num(e))
  if (!is.call(e)) stop("cannot translate `", deparse(e), "` to Julia", call. = FALSE)
  fn <- e[[1]]
  if (!is.symbol(fn)) stop("cannot translate `", deparse(e), "` to Julia", call. = FALSE)
  fn <- as.character(fn)
  args <- as.list(e)[-1]
  if (fn == "(") return(paste0("(", jl_expr(args[[1]]), ")"))
  if (fn %in% .infix && length(args) == 2) {
    return(paste(jl_expr(args[[1]]), fn, jl_expr(args[[2]])))
  }
  if (fn %in% c("-", "+", "~") && length(args) == 1) return(paste0(fn, jl_expr(args[[1]])))
  if (grepl("^%.*%$", fn)) {
    stop("operator ", fn, " is not supported in formulas (for nesting write a + a:b)", call. = FALSE)
  }
  if (fn %in% c("{", "if", "for", "while", "function", "<-", "=", "[", "[[", "$", "@")) {
    stop("`", deparse(e)[1], "` is not supported in a formula", call. = FALSE)
  }
  nms <- names(args)
  if (is.null(nms)) nms <- rep("", length(args))
  parts <- vapply(seq_along(args), function(i) {
    a <- jl_expr(args[[i]])
    if (nzchar(nms[i])) paste0(nms[i], " = ", a) else a
  }, "")
  paste0(fn, "(", paste(parts, collapse = ", "), ")")
}

# ------------------------------------------------------------------ formulas

as_formula_obj <- function(f) {
  if (inherits(f, "formula")) return(f)
  if (is.character(f) && length(f) == 1) return(stats::as.formula(f, env = globalenv()))
  stop("expected a formula, got ", class(f)[1], call. = FALSE)
}

#' Model formulas
#'
#' Set up a model formula in the syntax of 'brms': a main formula
#' `response ~ terms`, optionally with formulas for distributional parameters
#' (`sigma ~ x`, `shape ~ ...`). Formulas can be added (`bf(y1 ~ x) + bf(y2 ~ x)`)
#' for a multivariate model, with [set_rescor()] to model residual
#' correlations.
#'
#' Formula syntax is that of 'brms': addition terms on the left
#' (`y | trials(n)`, `y | weights(w)`, `y | cens(c)`), group-level terms
#' `(1 + x | g)`, `(1 | g1:g2)`, `(x || g)`, and special terms such as
#' `s()`, `gp()`, `mo()`, `ar()`. Functions inside the formula (`log(x)`,
#' `I(x^2)`) are evaluated by Julia, so use functions Julia has under the same
#' name (`log`, `exp`, `sqrt`, `abs`, ...). Variable names may contain only
#' letters, digits and `_`.
#'
#' @param formula The main formula, `response ~ terms`, or a string.
#' @param ... Formulas for distributional or non-linear parameters.
#' @param nl Is the main formula non-linear? Then its right-hand side is an
#'   ordinary expression and each parameter gets its own formula.
#' @param center Centre the population-level design matrix (the default, as in
#'   'brms').
#' @return An object of class `pbr_formula`.
#' @examples
#' bf(y ~ x + (1 | g), sigma ~ x)
#' bf(count | trials(n) ~ treatment + (1 | site))
#' @export
bf <- function(formula, ..., nl = FALSE, center = TRUE) {
  if (inherits(formula, "pbr_formula")) return(formula)
  fs <- c(list(as_formula_obj(formula)), lapply(list(...), as_formula_obj))
  structure(list(formulas = fs, nl = isTRUE(nl), center = !isFALSE(center)), class = "pbr_formula")
}

#' @rdname bf
#' @export
brmsformula <- bf

#' Linear formulas for distributional parameters
#'
#' `bf(y ~ x) + lf(sigma ~ z)` is `bf(y ~ x, sigma ~ z)`.
#' @param ... Formulas.
#' @return An object of class `pbr_lf`.
#' @export
lf <- function(...) structure(list(formulas = lapply(list(...), as_formula_obj)), class = "pbr_lf")

#' Residual correlations in multivariate models
#'
#' Add to a multivariate formula: `bf(y1 ~ x) + bf(y2 ~ x) + set_rescor(TRUE)`.
#' @param rescor Model residual correlations between the (gaussian or
#'   student) responses?
#' @return An object of class `pbr_rescor`.
#' @export
set_rescor <- function(rescor = TRUE) structure(list(rescor = isTRUE(rescor)), class = "pbr_rescor")

#' @export
`+.pbr_formula` <- function(e1, e2) {
  if (inherits(e2, "pbr_lf")) {
    e1$formulas <- c(e1$formulas, e2$formulas)
    return(e1)
  }
  if (inherits(e2, "formula")) e2 <- bf(e2)
  if (inherits(e2, "pbr_formula")) {
    return(structure(list(parts = list(e1, e2), rescor = NULL), class = "pbr_mvformula"))
  }
  if (inherits(e2, "pbr_rescor")) {
    return(structure(list(parts = list(e1), rescor = e2$rescor), class = "pbr_mvformula"))
  }
  stop("cannot add a ", class(e2)[1], " to a formula", call. = FALSE)
}

#' @export
`+.pbr_mvformula` <- function(e1, e2) {
  if (inherits(e2, "formula")) e2 <- bf(e2)
  if (inherits(e2, "pbr_formula")) {
    e1$parts <- c(e1$parts, list(e2))
  } else if (inherits(e2, "pbr_rescor")) {
    e1$rescor <- e2$rescor
  } else {
    stop("cannot add a ", class(e2)[1], " to a multivariate formula", call. = FALSE)
  }
  e1
}

#' @export
print.pbr_formula <- function(x, ...) {
  for (f in x$formulas) cat(paste(deparse(f, width.cutoff = 500L), collapse = " "), "\n")
  if (x$nl) cat("(non-linear)\n")
  invisible(x)
}

#' @export
print.pbr_mvformula <- function(x, ...) {
  for (p in x$parts) print(p)
  if (!is.null(x$rescor)) cat("rescor:", x$rescor, "\n")
  invisible(x)
}

as_pbr_formula <- function(formula) {
  if (inherits(formula, c("pbr_formula", "pbr_mvformula"))) return(formula)
  if (inherits(formula, "brmsformula")) {
    stop("a 'brms' bf() object: use this package's bf() (brms::bf objects carry R-only parts)",
         call. = FALSE)
  }
  bf(formula)
}

formula_strings <- function(f) {
  vapply(f$formulas, function(x) {
    if (length(x) != 3) stop("formulas need a left-hand side: `", deparse(x), "`", call. = FALSE)
    jl_expr(x)
  }, "")
}

#' Julia source for a formula.
#' @noRd
jl_formula <- function(f) {
  f <- as_pbr_formula(f)
  if (inherits(f, "pbr_mvformula")) {
    s <- paste(vapply(f$parts, jl_formula, ""), collapse = " + ")
    if (!is.null(f$rescor)) s <- paste0(s, " + PBR.set_rescor(", tolower(f$rescor), ")")
    return(s)
  }
  paste0("PBR.bf(", paste(vapply(formula_strings(f), jl_str, ""), collapse = ", "),
         "; nl = ", tolower(f$nl), ", center = ", tolower(f$center), ")")
}

# the variables a formula reads (for the data columns to send)
formula_vars <- function(f) {
  f <- as_pbr_formula(f)
  fs <- if (inherits(f, "pbr_mvformula")) unlist(lapply(f$parts, `[[`, "formulas")) else f$formulas
  unique(unlist(lapply(fs, all.vars)))
}

# ------------------------------------------------------------------ families

# name, default link (PracticalBayesianRegression's family table)
.family_table <- c(
  gaussian = "identity", student = "identity", lognormal = "identity", bernoulli = "logit",
  binomial = "logit", beta_binomial = "logit", poisson = "log", negbinomial = "log", geometric = "log",
  gamma = "log", weibull = "log", exponential = "log", beta = "logit",
  zero_inflated_poisson = "log", zero_inflated_negbinomial = "log", zero_inflated_binomial = "logit",
  hurdle_poisson = "log", hurdle_negbinomial = "log", hurdle_gamma = "log", hurdle_lognormal = "identity",
  cumulative = "logit", sratio = "logit", cratio = "logit", acat = "logit", categorical = "logit",
  skew_normal = "identity", exgaussian = "identity", inverse.gaussian = "inverse_squared",
  frechet = "log", shifted_lognormal = "identity", asym_laplace = "identity", von_mises = "tan_half",
  zero_inflated_beta = "logit", zero_one_inflated_beta = "logit", discrete_weibull = "logit",
  negbinomial2 = "log", com_poisson = "log", hurdle_cumulative = "logit", multinomial = "logit",
  cox = "log", logistic_normal = "identity", gen_extreme_value = "identity", xbeta = "logit",
  dirichlet = "logit"
)

#' Model families
#'
#' Response distributions and link functions, named as in 'brms'.
#' `brmsfamily("negbinomial", link = "log")` names any supported family; the
#' constructors below are shorthands. Families from the stats package
#' (`gaussian()`, `binomial()`, `poisson()`, `Gamma()`, `inverse.gaussian()`)
#' work as well, as do 'brms' family objects.
#'
#' One difference from 'brms': `brmsfamily("gamma")` without a link uses the
#' log link. `Gamma()` (stats) carries its inverse link explicitly, as in
#' 'brms'.
#'
#' @param family Name of the family.
#' @param link Link function of the mean; `NULL` for the family's default.
#' @param ... Families (for `mixture()`), as names or family objects.
#' @return An object of class `brmsfamily`.
#' @export
brmsfamily <- function(family, link = NULL) {
  structure(list(family = family, link = link), class = c("brmsfamily", "family"))
}

#' @rdname brmsfamily
#' @export
student <- function(link = "identity") brmsfamily("student", link)
#' @rdname brmsfamily
#' @export
bernoulli <- function(link = "logit") brmsfamily("bernoulli", link)
#' @rdname brmsfamily
#' @export
negbinomial <- function(link = "log") brmsfamily("negbinomial", link)
#' @rdname brmsfamily
#' @export
geometric <- function(link = "log") brmsfamily("geometric", link)
#' @rdname brmsfamily
#' @export
lognormal <- function(link = "identity") brmsfamily("lognormal", link)
#' @rdname brmsfamily
#' @export
shifted_lognormal <- function(link = "identity") brmsfamily("shifted_lognormal", link)
#' @rdname brmsfamily
#' @export
skew_normal <- function(link = "identity") brmsfamily("skew_normal", link)
#' @rdname brmsfamily
#' @export
exgaussian <- function(link = "identity") brmsfamily("exgaussian", link)
#' @rdname brmsfamily
#' @export
weibull <- function(link = "log") brmsfamily("weibull", link)
#' @rdname brmsfamily
#' @export
exponential <- function(link = "log") brmsfamily("exponential", link)
#' @rdname brmsfamily
#' @export
frechet <- function(link = "log") brmsfamily("frechet", link)
#' @rdname brmsfamily
#' @export
Beta <- function(link = "logit") brmsfamily("beta", link)
#' @rdname brmsfamily
#' @export
beta_binomial <- function(link = "logit") brmsfamily("beta_binomial", link)
#' @rdname brmsfamily
#' @export
von_mises <- function(link = "tan_half") brmsfamily("von_mises", link)
#' @rdname brmsfamily
#' @export
asym_laplace <- function(link = "identity") brmsfamily("asym_laplace", link)
#' @rdname brmsfamily
#' @export
gen_extreme_value <- function(link = "identity") brmsfamily("gen_extreme_value", link)
#' @rdname brmsfamily
#' @export
zero_inflated_poisson <- function(link = "log") brmsfamily("zero_inflated_poisson", link)
#' @rdname brmsfamily
#' @export
zero_inflated_negbinomial <- function(link = "log") brmsfamily("zero_inflated_negbinomial", link)
#' @rdname brmsfamily
#' @export
zero_inflated_binomial <- function(link = "logit") brmsfamily("zero_inflated_binomial", link)
#' @rdname brmsfamily
#' @export
zero_inflated_beta <- function(link = "logit") brmsfamily("zero_inflated_beta", link)
#' @rdname brmsfamily
#' @export
zero_one_inflated_beta <- function(link = "logit") brmsfamily("zero_one_inflated_beta", link)
#' @rdname brmsfamily
#' @export
hurdle_poisson <- function(link = "log") brmsfamily("hurdle_poisson", link)
#' @rdname brmsfamily
#' @export
hurdle_negbinomial <- function(link = "log") brmsfamily("hurdle_negbinomial", link)
#' @rdname brmsfamily
#' @export
hurdle_gamma <- function(link = "log") brmsfamily("hurdle_gamma", link)
#' @rdname brmsfamily
#' @export
hurdle_lognormal <- function(link = "identity") brmsfamily("hurdle_lognormal", link)
#' @rdname brmsfamily
#' @export
cumulative <- function(link = "logit") brmsfamily("cumulative", link)
#' @rdname brmsfamily
#' @export
sratio <- function(link = "logit") brmsfamily("sratio", link)
#' @rdname brmsfamily
#' @export
cratio <- function(link = "logit") brmsfamily("cratio", link)
#' @rdname brmsfamily
#' @export
acat <- function(link = "logit") brmsfamily("acat", link)
#' @rdname brmsfamily
#' @export
categorical <- function(link = "logit") brmsfamily("categorical", link)
#' @rdname brmsfamily
#' @export
mixture <- function(...) {
  fams <- lapply(list(...), function(f) resolve_family(f)$family)
  brmsfamily(paste0("mixture(", paste(fams, collapse = ", "), ")"))
}

#' @export
print.brmsfamily <- function(x, ...) {
  f <- resolve_family(x)
  cat("Family:", f$family, "\nLink function:", f$link, "\n")
  invisible(x)
}

# -> list(family = Julia name, link = Julia link name)
resolve_family <- function(family) {
  if (is.function(family)) family <- family()
  if (is.character(family) && length(family) == 1) family <- list(family = family, link = NULL)
  if (!is.list(family) || is.null(family$family)) {
    stop("`family` must be a family object such as poisson() or a family name", call. = FALSE)
  }
  name <- family$family
  if (grepl("^quasi", name)) stop("quasi families are not likelihoods; not supported", call. = FALSE)
  if (!startsWith(name, "mixture(")) name <- tolower(name)
  link <- family$link
  if (!is.null(link)) link <- switch(link, "1/mu^2" = "inverse_squared", link)
  if (is.null(link) && name %in% names(.family_table)) link <- .family_table[[name]]
  if (!startsWith(name, "mixture(") && !name %in% names(.family_table)) {
    stop("family `", name, "` is not supported", call. = FALSE)
  }
  list(family = name, link = link)
}

jl_family <- function(family) {
  f <- resolve_family(family)
  if (is.null(f$link)) return(paste0("PBR.brmsfamily(", jl_str(f$family), ")"))
  paste0("PBR.brmsfamily(", jl_str(f$family), "; link = ", jl_str(f$link), ")")
}
