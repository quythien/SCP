read_app_functions <- function() {
  root <- Sys.getenv("SCP_TEST_SOURCE_ROOT")
  path <- if (nzchar(root)) file.path(root, "inst", "shiny", "app.R") else
    system.file("shiny", "app.R", package = "SCP")
  stopifnot(file.exists(path))
  environment <- new.env(parent = asNamespace("SCP"))
  wanted <- c("%||%", ".app_explorer_table", ".app_apply_threshold")
  for (expression in parse(path)) {
    if (is.call(expression) && identical(expression[[1]], as.name("<-")) &&
        is.symbol(expression[[2]]) && as.character(expression[[2]]) %in% wanted) {
      eval(expression, envir = environment)
    }
  }
  stopifnot(all(vapply(wanted, exists, logical(1), envir = environment,
                        inherits = FALSE)))
  environment
}
