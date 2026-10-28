test_that("inline R syntax distinguishes code from mathematical text", {
  bad <- tempfile(fileext = ".Rmd")
  writeLines("The signal ratio is `r = A / sigma`.", bad)
  expect_length(inline_r_parse_errors(bad), 1L)
  good <- tempfile(fileext = ".Rmd")
  writeLines(c("The signal ratio is `A / sigma`.", "The value is `r 1 + 1`."), good)
  expect_length(inline_r_parse_errors(good), 0L)
})

test_that("vignettes contain parseable inline R expressions", {
  source <- testthat::test_path("..", "..", "vignettes")
  paths <- list.files(source, pattern = "[.]Rmd$", full.names = TRUE)
  if (!length(paths)) {
    paths <- list.files(system.file("doc", package = "SCP"),
                         pattern = "[.]Rmd$", full.names = TRUE)
  }
  skip_if(!length(paths), "Vignette sources are checked in the source CI and release build.")
  expect_length(inline_r_parse_errors(paths), 0L)
})
