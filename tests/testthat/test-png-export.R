read_png_opener <- function() {
  root <- Sys.getenv("SCP_TEST_SOURCE_ROOT")
  path <- if (nzchar(root)) file.path(root, "inst", "shiny", "app.R") else
    system.file("shiny", "app.R", package = "SCP")
  assignments <- Filter(function(x) is.call(x) &&
    identical(x[[1]], as.name("<-")) && identical(x[[2]], as.name("server")),
    as.list(parse(path)))
  statements <- as.list(assignments[[1]][[3]][[3]])
  opener <- Filter(function(x) is.call(x) &&
    identical(x[[1]], as.name("<-")) && identical(x[[2]], as.name(".open_png")),
    statements)
  stopifnot(length(opener) == 1L)
  environment <- new.env(parent = baseenv())
  eval(opener[[1]], envir = environment)
  environment
}

test_that("the app writes a PNG with the requested dimensions", {
  app <- read_png_opener()
  file <- tempfile(fileext = ".png")
  on.exit(unlink(file), add = TRUE)
  previous <- grDevices::dev.cur()
  app$.open_png(file, width = 640, height = 480, res = 96)
  expect_false(grDevices::dev.cur() == previous)
  graphics::plot(1:3, 1:3)
  grDevices::dev.off()
  header <- readBin(file, what = "raw", n = 24L)
  expect_identical(header[1:8], as.raw(c(137, 80, 78, 71, 13, 10, 26, 10)))
  width <- sum(as.integer(header[17:20]) * 256^(3:0))
  height <- sum(as.integer(header[21:24]) * 256^(3:0))
  expect_equal(c(width, height), c(640, 480))
})

test_that("an unavailable PNG device fails clearly", {
  app <- read_png_opener()
  app$requireNamespace <- function(...) FALSE
  app$capabilities <- function(...) TRUE
  local_mocked_bindings(png = function(...) warning("PNG library unavailable"),
                         .package = "grDevices")
  expect_error(app$.open_png(tempfile(fileext = ".png"), 640, 480),
               "Unable to open a PNG device")
})
