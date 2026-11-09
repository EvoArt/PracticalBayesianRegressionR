# Priors, as in 'brms': a data frame with one row per prior statement.

.prior_cols <- c("prior", "class", "coef", "group", "resp", "dpar", "nlpar", "lb", "ub")

#' Prior distributions
#'
#' Specify priors as 'brms' does. `set_prior()` (and its alias
#' `prior_string()`) takes the distribution as a string; `prior()` takes it
#' unquoted, along with unquoted class, coefficient and group names. Combine
#' several with `c()` or `+`. Defaults are 'brms''s default priors; see
#' [get_prior()].
#'
#' Supported distributions: `normal`, `student_t`, `cauchy`, `exponential`,
#' `gamma`, `inv_gamma`, `lognormal`, `beta`, `uniform`, `logistic`,
#' `double_exponential`, `skew_normal`, `constant`, `lkj` (class `cor`),
#' `dirichlet`, `horseshoe`, `R2D2`.
#'
#' @param prior The prior distribution, e.g. `"normal(0, 1)"`.
#' @param class Parameter class: `"b"`, `"Intercept"`, `"sd"`, `"cor"`,
#'   `"sigma"`, `"shape"`, ...
#' @param coef Coefficient name (within class `b` or `sd`).
#' @param group Grouping factor (classes `sd`, `cor`).
#' @param resp Response (multivariate models).
#' @param dpar Distributional parameter the prior applies to.
#' @param nlpar Non-linear parameter the prior applies to.
#' @param lb,ub Bounds on population-level coefficients (class `b`).
#' @param check Ignored (for compatibility with 'brms').
#' @return An object of class `pbr_prior`: a data frame.
#' @examples
#' set_prior("normal(0, 1)", class = "b")
#' prior(student_t(3, 0, 2.5), class = sd, group = g) + prior(lkj(2), class = cor)
#' @export
set_prior <- function(prior, class = "b", coef = "", group = "", resp = "", dpar = "", nlpar = "",
                      lb = NA, ub = NA, check = TRUE) {
  if (length(prior) == 0) {
    df <- as.data.frame(stats::setNames(rep(list(character(0)), length(.prior_cols)), .prior_cols),
                        stringsAsFactors = FALSE)
    return(structure(df, class = c("pbr_prior", "data.frame")))
  }
  df <- data.frame(prior = as.character(prior), class = as.character(class), coef = as.character(coef),
                   group = as.character(group), resp = as.character(resp), dpar = as.character(dpar),
                   nlpar = as.character(nlpar), lb = as.character(lb), ub = as.character(ub),
                   stringsAsFactors = FALSE)
  structure(df, class = c("pbr_prior", "data.frame"))
}

#' @rdname set_prior
#' @export
prior_string <- set_prior

#' @rdname set_prior
#' @export
prior <- function(prior, class = "b", coef = "", group = "", resp = "", dpar = "", nlpar = "",
                  lb = NA, ub = NA, check = TRUE) {
  call <- as.list(match.call())[-1]
  txt <- function(nm, default) {
    e <- call[[nm]]
    if (is.null(e)) return(default)
    if (is.character(e) || is.numeric(e) || is.logical(e)) return(as.character(e))
    paste(deparse(e, width.cutoff = 500L), collapse = " ")
  }
  set_prior(txt("prior", ""), class = txt("class", "b"), coef = txt("coef", ""), group = txt("group", ""),
            resp = txt("resp", ""), dpar = txt("dpar", ""), nlpar = txt("nlpar", ""),
            lb = txt("lb", NA), ub = txt("ub", NA))
}

#' @export
c.pbr_prior <- function(x, ...) {
  parts <- c(list(x), list(...))
  df <- do.call(rbind, lapply(parts, function(p) as.data.frame(as_pbr_prior(p))))
  rownames(df) <- NULL
  structure(df, class = c("pbr_prior", "data.frame"))
}

#' @export
`+.pbr_prior` <- function(e1, e2) c(e1, e2)

#' @export
print.pbr_prior <- function(x, ...) {
  if (nrow(x) == 0) {
    cat("(no priors: brms defaults)\n")
    return(invisible(x))
  }
  y <- as.data.frame(x)
  y$lb[is.na(y$lb)] <- ""
  y$ub[is.na(y$ub)] <- ""
  print(y, row.names = FALSE, right = FALSE)
  invisible(x)
}

# NULL, a pbr_prior or a 'brms' brmsprior (same columns) -> pbr_prior
as_pbr_prior <- function(p) {
  if (is.null(p)) return(set_prior(character(0)))
  if (is.character(p)) return(set_prior(p))
  if (!is.data.frame(p) || !"prior" %in% names(p)) {
    stop("`prior` must come from prior(), set_prior() or prior_string()", call. = FALSE)
  }
  p <- as.data.frame(p, stringsAsFactors = FALSE)
  for (nm in .prior_cols) {
    if (!nm %in% names(p)) p[[nm]] <- if (nm %in% c("lb", "ub")) NA_character_ else ""
    p[[nm]] <- as.character(p[[nm]])
    if (!nm %in% c("lb", "ub")) p[[nm]][is.na(p[[nm]])] <- ""
  }
  p$lb[p$lb %in% ""] <- NA
  p$ub[p$ub %in% ""] <- NA
  structure(p[, .prior_cols], class = c("pbr_prior", "data.frame"))
}

jl_prior <- function(p) {
  p <- as_pbr_prior(p)
  if (nrow(p) == 0) return("PBR.PriorSpec[]")
  rows <- vapply(seq_len(nrow(p)), function(i) {
    r <- p[i, ]
    kw <- c(paste0("class = ", jl_str(r$class)))
    for (nm in c("coef", "group", "resp", "dpar", "nlpar")) {
      if (nzchar(r[[nm]])) kw <- c(kw, paste0(nm, " = ", jl_str(r[[nm]])))
    }
    for (nm in c("lb", "ub")) {
      if (!is.na(r[[nm]])) kw <- c(kw, paste0(nm, " = ", jl_num(as.numeric(r[[nm]]))))
    }
    paste0("PBR.prior(", jl_str(r$prior), "; ", paste(kw, collapse = ", "), ")")
  }, "")
  paste0("PBR.PriorSpec[", paste(rows, collapse = ", "), "]")
}
