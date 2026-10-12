# =============================================================================
# Paths used by every script. Each script finds the repository root (the folder
# containing this file) and sources it, so no script contains a machine-specific path.
#
#   RAW       public input data (data/raw; filled by scripts/00_download_data.R)
#   PCCF_DIR  licensed PCCF+ Version 8D files (data/licensed/PCCF8D; not distributed)
#   DERIVED   intermediate and output files written by the pipeline (data/derived)
#
# To keep data elsewhere (for example on an external drive), copy
# config_local_example.R to config_local.R and edit it. config_local.R is not
# tracked by git.
# =============================================================================

RAW      <- file.path(ROOT, "data", "raw")
PCCF_DIR <- file.path(ROOT, "data", "licensed", "PCCF8D")
DERIVED  <- file.path(ROOT, "data", "derived")

if (file.exists(file.path(ROOT, "config_local.R"))) source(file.path(ROOT, "config_local.R"))
RAW <- path.expand(RAW); PCCF_DIR <- path.expand(PCCF_DIR); DERIVED <- path.expand(DERIVED)
dir.create(DERIVED, recursive = TRUE, showWarnings = FALSE)
