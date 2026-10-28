# Run the demo with the same app code as the installed package.
version_file <- file.path(getwd(), "SCP_VERSION")
expected_version <- Sys.getenv("SCP_SHINY_PACKAGE_VERSION",
  if (file.exists(version_file)) readLines(version_file, n = 1L) else "0.5.1")
if (!requireNamespace("SCP", quietly = TRUE) ||
    as.character(utils::packageVersion("SCP")) != expected_version) {
  stop("Install SCP ", expected_version, " before launching this demo.")
}
app <- system.file("shiny", "app.R", package = "SCP")
if (!nzchar(app)) stop("The installed SCP app is missing.")
checksum_file <- file.path(getwd(), "SCP_APP_MD5")
if (file.exists(checksum_file) &&
    !identical(readLines(checksum_file, n = 1L), unname(tools::md5sum(app))))
  stop("The installed SCP app does not match this demo bundle. Install its recorded source revision.")
bundle <- file.path(getwd(), "pilots")
old_demo <- Sys.getenv("SCP_SHINY_DEMO", unset = NA_character_)
old_bundle <- Sys.getenv("SCP_SHINY_DEMO_PILOTS", unset = NA_character_)
Sys.setenv(SCP_SHINY_DEMO = if (file.exists(file.path(bundle, "manifest.csv"))) "1" else "",
           SCP_SHINY_DEMO_PILOTS = bundle)
tryCatch(source(app, local = TRUE)$value, finally = {
  if (is.na(old_demo)) Sys.unsetenv("SCP_SHINY_DEMO") else
    Sys.setenv(SCP_SHINY_DEMO = old_demo)
  if (is.na(old_bundle)) Sys.unsetenv("SCP_SHINY_DEMO_PILOTS") else
    Sys.setenv(SCP_SHINY_DEMO_PILOTS = old_bundle)
})
