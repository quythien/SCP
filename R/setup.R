# setup.R -- Legacy source-and-run entry point (retained for backward compatibility).
#
# When using the installed SCP package, this file is NOT needed. All dependencies
# are declared in DESCRIPTION (Imports/Suggests) and loaded automatically by R.
#
# For scripts that still use source("code/setup.R"), the original file in code/
# remains intact. This copy inside R/ is inert when the package is installed and
# is included only so that R CMD check does not flag a missing file reference.
#
# To load the package in a fresh session: library(SCP)
