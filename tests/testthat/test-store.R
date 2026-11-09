test_that("model.toml's top-level keys are read", {
  f <- tempfile()
  writeLines(c("# comment", "arch = \"x86_64\"", "columns = [\"g (factor)\", \"x (numeric)\"]",
               "keys = [\"julia\\nformula: y~x\\n\", \"R\\nform \\\"q\\\" \\u00e9\"]", "n = 3", "ok = true",
               "empty = []", "multi = [", "  \"a\",", "  \"b\",", "]", "", "[[prior]]", "class = \"b\""), f)
  m <- pbr$read_toml_top(f)
  expect_equal(m$arch, "x86_64")
  expect_equal(m$columns, c("g (factor)", "x (numeric)"))
  expect_equal(m$keys, c("julia\nformula: y~x\n", "R\nform \"q\" \u00e9"))
  expect_equal(m$n, 3)
  expect_true(m$ok)
  expect_length(m$empty, 0)
  expect_equal(m$multi, c("a", "b"))
  expect_null(m$class)                                # tables are not read
})

test_that("the store lists, finds by key, and skips foreign and half-built models", {
  store <- tempfile("store")
  d <- sim_data()
  key <- pbr$r_model_key(bf(y ~ x + (1 | g)), d, gaussian(), prior(normal(0, 2), class = b))
  fake_model(store, "m_new", keys = c("julia\nwhatever\n", key), created = "2026-11-09 13:00:00", exe = FALSE)
  fake_model(store, "m_old", created = "2026-11-01 13:00:00", exe = FALSE)
  fake_model(store, "m_foreign", keys = key, os = "plan9", exe = FALSE)
  dir.create(file.path(store, ".build_m_x_123"))
  s <- compiled_models(store)
  expect_equal(names(s), c("m_new", "m_old"))         # newest first; plan9 and the build dir skipped
  expect_output(print(s), "2 compiled models")
  expect_equal(compiled_model("m_old", store)$family, "gaussian")
  expect_error(compiled_model("m_foreign", store), "compiled for plan9")
  expect_error(compiled_model("nope", store), "no compiled model")
  # the same model however it is written; another prior, family or column type is another model
  expect_equal(find_compiled_model("y ~ x + (1|g)", d, "gaussian", set_prior("normal(0,2)"), store)$name,
               "m_new")
  expect_null(find_compiled_model(y ~ x + (1 | g), d, gaussian(), NULL, store))
  expect_null(find_compiled_model(y ~ x + (1 | g), d, student(), prior(normal(0, 2)), store))
  d2 <- d
  d2$g <- as.numeric(factor(d2$g))
  expect_null(find_compiled_model(y ~ x + (1 | g), d2, gaussian(), prior(normal(0, 2)), store))
  # compile_model returns a stored model without starting Julia
  expect_equal(compile_model(y ~ x + (1 | g), d, gaussian(), prior(normal(0, 2), class = b),
                             store = store)$name, "m_new")
  remove_compiled_model("m_old", store)
  expect_equal(names(compiled_models(store)), "m_new")
})

test_that("data are checked against a compiled model's columns", {
  store <- tempfile("store")
  fake_model(store, "m", exe = FALSE)
  m <- compiled_model("m", store)
  d <- sim_data()
  expect_true(pbr$check_model_data(m, d))
  expect_error(pbr$check_model_data(m, d[c("y", "x")]), "lack column\\(s\\).*g")
  d$x <- as.character(d$x)
  expect_error(pbr$check_model_data(m, d), "`x` must be numeric")
})

test_that("priors stored with a model come back as a prior table", {
  store <- tempfile("store")
  fake_model(store, "m", exe = FALSE, priors = c("b||||||normal(0, 2)||", "sd||g||||student_t(3, 0, 1)|0|"))
  p <- pbr$model_priors(compiled_model("m", store))
  expect_equal(p$prior, c("normal(0, 2)", "student_t(3, 0, 1)"))
  expect_equal(p$group, c("", "g"))
  expect_equal(p$lb, c(NA, "0"))
})
