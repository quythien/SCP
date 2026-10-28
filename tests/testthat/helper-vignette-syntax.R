inline_r_parse_errors <- function(paths) {
  errors <- character()
  for (path in paths) {
    lines <- readLines(path, warn = FALSE)
    for (i in seq_along(lines)) {
      matches <- gregexpr("`r[ \t]+[^`\r\n]*`", lines[i], perl = TRUE)
      tokens <- regmatches(lines[i], matches)[[1]]
      for (token in tokens) {
        code <- sub("^`r[ \t]+", "", substr(token, 1, nchar(token) - 1L))
        error <- tryCatch({ parse(text = code); NULL }, error = identity)
        if (!is.null(error)) {
          errors <- c(errors, sprintf("%s:%d: %s", basename(path), i,
                                      conditionMessage(error)))
        }
      }
    }
  }
  errors
}
