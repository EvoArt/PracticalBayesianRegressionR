test_that("data CSV round-trips exactly and quotes text", {
  d <- data.frame(x = c(0.1, 1 / 3, -2e-300, NA, 1e15 + 0.5), i = c(1L, NA, 3L, 4L, 5L),
                  s = c("a", "b,c", "say \"hi\"", NA, "z"), b = c(TRUE, FALSE, NA, TRUE, TRUE),
                  f = factor(c("u", "v", "u", "v", "u")), stringsAsFactors = FALSE)
  file <- tempfile(fileext = ".csv")
  pbr$write_data_csv(d, file)
  back <- utils::read.csv(file, stringsAsFactors = FALSE, na.strings = "NA")
  expect_identical(back$x, d$x)                      # every bit
  expect_identical(back$i, d$i)
  expect_equal(back$s, d$s)
  expect_equal(back$b, d$b)
  expect_equal(back$f, as.character(d$f))
  expect_equal(readLines(file)[1], "x,i,s,b,f")
  expect_equal(pbr$col_kind(d$f), "string")
  expect_equal(pbr$col_kind(d$i), "integer")
  expect_equal(pbr$col_kind(d$b), "logical")
})

test_that("factor level order other than sorted is warned about", {
  d <- data.frame(f = factor(c("b", "a"), levels = c("b", "a")))
  expect_warning(pbr$warn_factor_levels(d), "reference level is `a`")
  expect_silent(pbr$warn_factor_levels(data.frame(f = factor(c("b", "a")))))
})
