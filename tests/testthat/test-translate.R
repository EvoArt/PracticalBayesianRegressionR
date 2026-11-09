test_that("formulas translate to Julia text", {
  jl <- function(e) pbr$jl_expr(e)
  expect_equal(jl(quote(y ~ x + (1 | g))), "y ~ x + (1 | g)")
  expect_equal(jl(quote(y | trials(n) ~ x * z + (1 + x || g))), "y | trials(n) ~ x * z + (1 + x || g)")
  expect_equal(jl(quote(y ~ s(x, k = 5L, bs = "cr"))), "y ~ s(x, k = 5, bs = \"cr\")")
  expect_equal(jl(quote(y ~ gr(g, cor = FALSE))), "y ~ gr(g, cor = false)")
  expect_equal(jl(quote(y ~ I(x^2) + log(z) - 1)), "y ~ I(x ^ 2) + log(z) - 1")
  expect_equal(jl(quote(y ~ x + offset(0.25))), "y ~ x + offset(0.25)")
  expect_error(jl(quote(y ~ a %in% b)), "not supported")
  expect_error(jl(quote(y ~ log.dose)), "Julia names")
  expect_error(jl(quote(y ~ x[1])), "not supported")
})

test_that("bf builds single, distributional and multivariate formulas", {
  f <- bf(y ~ x + (1 | g), sigma ~ x)
  expect_s3_class(f, "pbr_formula")
  expect_equal(pbr$jl_formula(f),
               "PBR.bf(\"y ~ x + (1 | g)\", \"sigma ~ x\"; nl = false, center = true)")
  expect_equal(pbr$jl_formula(bf(y ~ x) + lf(sigma ~ x)), pbr$jl_formula(bf(y ~ x, sigma ~ x)))
  expect_equal(pbr$jl_formula("y ~ x"), "PBR.bf(\"y ~ x\"; nl = false, center = true)")
  mv <- bf(y1 ~ x) + bf(y2 ~ x) + set_rescor(TRUE)
  expect_s3_class(mv, "pbr_mvformula")
  expect_match(pbr$jl_formula(mv),
               "PBR.bf\\(\"y1 ~ x\".*\\) \\+ PBR.bf\\(\"y2 ~ x\".*\\) \\+ PBR.set_rescor\\(true\\)")
  expect_setequal(pbr$formula_vars(bf(y | trials(n) ~ x + (1 | g))), c("y", "n", "x", "g"))
  expect_equal(pbr$jl_str("a\"b$c\\d"), "\"a\\\"b\\$c\\\\d\"")
})

test_that("families resolve to Julia names and links", {
  r <- pbr$resolve_family
  expect_equal(r(poisson()), list(family = "poisson", link = "log"))
  expect_equal(r(Gamma()), list(family = "gamma", link = "inverse"))       # brms's Gamma default
  expect_equal(r("gamma"), list(family = "gamma", link = "log"))           # the Julia default
  expect_equal(r(inverse.gaussian())$link, "inverse_squared")
  expect_equal(r(negbinomial()), list(family = "negbinomial", link = "log"))
  expect_equal(r(poisson), list(family = "poisson", link = "log"))         # a bare function
  expect_equal(r(mixture(gaussian, "student"))$family, "mixture(gaussian, student)")
  expect_error(r(quasipoisson()), "quasi")
  expect_error(r("nonsense"), "not supported")
  expect_equal(pbr$jl_family(binomial(link = "probit")), "PBR.brmsfamily(\"binomial\"; link = \"probit\")")
})

test_that("priors: set_prior, prior() and combining", {
  p <- prior(normal(0, 1), class = b) + prior(student_t(3, 0, 2.5), class = sd, group = g) +
    set_prior("lkj(2)", class = "cor")
  expect_s3_class(p, "pbr_prior")
  expect_equal(p$prior, c("normal(0, 1)", "student_t(3, 0, 2.5)", "lkj(2)"))
  expect_equal(p$class, c("b", "sd", "cor"))
  expect_equal(p$group, c("", "g", ""))
  expect_equal(c(prior(normal(0, 1)), prior_string("normal(0, 2)", coef = "x"))$coef, c("", "x"))
  expect_equal(pbr$jl_prior(prior(normal(0, 1), class = b, lb = 0)),
               "PBR.PriorSpec[PBR.prior(\"normal(0, 1)\"; class = \"b\", lb = 0)]")
  expect_equal(pbr$jl_prior(NULL), "PBR.PriorSpec[]")
  # a brms prior data frame (brmsprior) has the same columns
  bp <- data.frame(prior = "normal(0, 5)", class = "b", coef = "", group = "", resp = "", dpar = "",
                   nlpar = "", lb = NA, ub = NA, source = "user")
  expect_equal(pbr$as_pbr_prior(bp)$prior, "normal(0, 5)")
})
